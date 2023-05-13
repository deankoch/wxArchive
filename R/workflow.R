#' List all times and variables available locally in the RAP/RUC/GFS archive
#'
#' This checks the sub-directories of `project_dir` (both 'rap' and 'gfs') for
#' NetCDF data on the variables named in `.nm_output_var`, and reports some stats
#' to the console. All file paths are returned in a list but only the existing files
#' are reported on.
#'
#' With `quiet=TRUE`, this is essentially just a wrapper for `nc_list`
#'
#' @param project_dir character path to the project root directory
#' @param quiet logical suppresses console messages
#'
#' @return a list of character vectors, paths to individual NetCDF files
#' @export
workflow_list = function(project_dir, quiet=FALSE) {

  # make a list of all datasets with preference for RAP archive over GFS
  p_all = project_dir |> nc_list()
  if( quiet ) return(p_all)

  # fixed width names for printout
  nm_all = names(p_all)
  nm_len = nm_all |> sapply(nchar) |> max()
  nm_fix_wid = lapply(nm_all, \(nm) paste0(nm, paste(rep(' ', nm_len - nchar(nm)), collapse='')))

  # check available times in a loop, printing as we go
  message('\nchecking available times for ', length(p_all), ' variable(s)')
  for( i in seq_along(p_all) ) {

    # check existence of source files
    p = p_all[[i]]
    p = p[ file.exists(p) ]

    # check number of source files
    n_file = p |> length()
    if( n_file == 0 ) {

      paste0('\n', nm_fix_wid[[i]], ' | 0 files') |> cat()
      next
    }

    # count times and find their range
    time_all = p |> time_wx()
    time_min = time_all[['time_obs']] |> min()
    time_max = time_all[['time_obs']] |> max()
    n = time_all[['time_obs']] |> length()

    # check if time series is complete
    ts_df = data.frame(posix_pred=time_all[['time_obs']]) |> archive_pad(quiet=TRUE)
    n_miss = ts_df[['ts_hours']] |> is.na() |> sum()
    msg_miss = ifelse(n_miss==0, ' (complete)', paste0(' (', n_miss, ' missing)'))

    # print stats to console
    paste0('\n', nm_fix_wid[[i]], ' | ',
           n_file, ' file(s) | ',
           n, ' time(s) | ',
           time_min, ' to ', time_max,
           msg_miss) |> cat()

    # print file paths to console
    c('\n', rep('-', nm_len-1)) |> paste(collapse='') |> cat()
    paste0('\n > ', p) |> cat()
    cat('\n')
  }

  cat('\n')
  return( invisible(p_all) )
}


#' Update the RUC/RAP archive
#'
#' Download the latest files, convert to NetCDF and replace some of the
#' missing layers. This is a wrapper for a sequence of calls to `archive_update`,
#' `nc_update`, `pcp_update`, then `nc_resample`, in that order.
#'
#' A geometry object must be supplied to define the bounding box for NetCDF data.
#' This function expects to find it in the file 'aoi.geojson' in `project_dir`.
#'
#' This creates/modifies several sub-directories of `project_dir` (named in
#' `.nm_resample_rap`), writing NetCDF versions of the variables named in
#' `.rap_regex` and `.var_pcp`.
#'
#' @param project_dir character path to the project root directory
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_update_rap = function(project_dir, from=NULL, to=NULL) {

  # path to area of interest polygon
  aoi_path = project_dir |> file.path('aoi.geojson')
  aoi = sf::st_read(aoi_path)

  # all the output files go here
  base_dir_rap = project_dir |> file.path('rap')

  ## RAP/RUC

  # download RUC/RAP files
  cat('\n')
  message('updating RAP/RUC GRIB archive')
  archive_update(base_dir = base_dir_rap,
                 hour_rel = .hour_rel_rap,
                 from = from,
                 to = to,
                 model = 'rap_archive') |> invisible()

  # export to nc
  cat('\n')
  message('updating NetCDF files')
  nc_update(aoi = aoi,
            base_dir = base_dir_rap,
            output_nm = .nm_src_rap,
            regex = .rap_regex) |> invisible()

  # copy precip from components (applies to early years)
  cat('\n')
  message('processing precipitation layers')
  pcp_update(base_dir = base_dir_rap,
             pcp_nm = .var_pcp,
             input_nm = .nm_src_rap,
             output_nm = .nm_rap) |> invisible()

  # resampled coarse to fine
  cat('\n')
  message('resampling')
  nc_resample(var_nm = .nm_output_var,
              base_dir = base_dir_rap,
              input_nm = .nm_src_rap,
              output_nm = .nm_resample) |> invisible()
}

#' Fit a temporal model to the RUC/RAP times series
#'
#' A wrapper for `time_fit`. This reads its input from the sub-directories
#' named in `.nm_resample_rap` and writes its output to the first of these
#' sub-directories (for each variable). See `?time_fit`
#'
#' Call this function at least once after running `workflow_update_rap` for
#' the first time (and before running `workflow_impute_rap`)
#'
#' @param project_dir character path to the project root directory
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_fit_temporal = function(project_dir) {

  # all the output files go here
  base_dir_rap = project_dir |> file.path('rap')

  # part 6: fit temporal model to fine grid (include all layers)
  time_fit(var_nm = .nm_output_var,
           base_dir = base_dir_rap,
           input_nm = .nm_resample_rap) |> invisible()
}


#' Impute missing time layers
#'
#' This creates the `.nm_complete` sub-directory of `project_dir` and
#' fills it with NetCDF files containing imputed layers for times missing
#' from the time series named in `.nm_resample_rap`.
#'
#' Users must first run `workflow_update_rap` and `workflow_fit_temporal` at
#' least once before calling this function.
#'
#' After running `workflow_update_rap` and `workflow_impute_rap`, users should
#' be able to access an up-to-date completed RAP time series (with no NA layers) by
#' passing `.nm_complete_rap` to `file_wx` then `nc_layers`.
#'
#' @param project_dir character path to the project root directory
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_impute_rap = function(project_dir) {

  # all the output files go here
  base_dir_rap = project_dir |> file.path('rap')

  # fill gaps
  cat('\n')
  message('imputing missing times')
  time_impute(var_nm = .nm_output_var,
              base_dir = base_dir_rap,
              input_nm = .nm_resample_rap,
              output_nm = .nm_complete) |> invisible()
}

#' Create wind speed layer
#'
#' This creates the `wnd` sub-directory of `project_dir` and writes a NetCDF file
#' for variable "wnd" after computing it from components "wnd_u" and "wnd_v".
#'
#' Users should either call this function after running the workflow up to
#' (and including) `workflow_impute_rap` to get a complete wind speed time series
#' from the gap-filled components.
#'
#' @param project_dir character path to the project root directory
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_wnd_rap = function(project_dir) {

  # all the output files go here
  base_dir_rap = project_dir |> file.path('rap')

  # create/update the file
  cat('\n')
  message('computing wind speed from u/v components')
  base_dir_rap |> wnd_update(wnd_nm = .var_wnd,
                             uv_nm = .var_wnd_uv,
                             input_nm = .nm_complete_rap) |> invisible()
}


#' Update the GFS archive
#'
#' Download the latest files, convert the times needed to NetCDF. This is a wrapper
#' for a sequence of calls to `archive_update`, `nc_update` in that order.
#'
#' Users should first run `workflow_update_rap` and `workflow_impute_rap`, in that
#' order, to update the RAP archive. The function checks the output from those functions
#' (sub-directories named in `.nm_complete_rap`) and identifies, for each variable. the
#' latest time available. Only times occurring after this point are extracted from GFS.
#'
#' A geometry object must be supplied to define the bounding box to fetch from GFS.
#' The function expects to find it in the file 'aoi.geojson' in `project_dir`.
#'
#' This will create/modify the 'gfs' sub-directory of `project_dir`, writing NetCDF
#' versions of the variables named in `.gfs_regex`.
#'
#' @param project_dir character path to the project root directory
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_update_gfs = function(project_dir) {

  # path to area of interest polygon
  aoi_path = project_dir |> file.path('aoi.geojson')
  aoi = sf::st_read(aoi_path, quiet=TRUE)

  # all the output files go here
  base_dir_gfs = project_dir |> file.path('gfs')

  ## GFS

  # RAP archive end time determines start time of the GFS output
  base_dir_rap = project_dir |> file.path('rap')
  rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
  rap_time = lapply(rap_nc_path, \(p) time_wx(p))
  from = do.call(c, lapply(rap_time, \(x) max(x[['time_obs']]) ))

  # download GFS files
  cat('\n')
  message('updating GFS GRIB archive')
  gfs_result = archive_update(base_dir = base_dir_gfs,
                              hour_pred = .hour_pred_gfs,
                              hour_rel = .hour_rel_gfs,
                              aoi = aoi,
                              model = 'gfs_0p25',
                              alternate = FALSE)

  # delete the old GFS NetCDF directories
  base_dir_gfs |> file.path('coarse') |> unlink(recursive=TRUE)
  base_dir_gfs |> file.path(.nm_resample) |> unlink(recursive=TRUE)
  base_dir_gfs |> file.path(.var_wnd) |> unlink(recursive=TRUE)

  # export latest GFS data to nc (creates "coarse" subdirectory)
  cat('\n')
  nc_update(aoi = aoi,
            base_dir = base_dir_gfs,
            output_nm = list(coarse=.nm_gfs),
            regex = .gfs_regex,
            from = from) |> invisible()

  # load an example grid at fine resolution (second rast call drops cell values)
  r_fine = file_wx('nc', base_dir_rap, .nm_spatial, names(.rap_regex)[[1]]) |>
    terra::rast() |> terra::rast()

  # resample to match RAP grid
  cat('\n')
  message('resampling')
  nc_resample(var_nm = .nm_gfs_var,
              base_dir = base_dir_gfs,
              input_nm = list(coarse=.nm_gfs),
              output_nm = .nm_resample,
              r_fine = r_fine) |> invisible()

  # create wind speed layers
  cat('\n')
  message('computing wind speed from u/v components')
  base_dir_gfs |> wnd_update(wnd_nm = 'wnd',
                             uv_nm = c('wnd_u', 'wnd_v'),
                             input_nm = .nm_resample) |> invisible()

}

#' Exported completed time series to daily aggregate values and write to disk
#'
#'
#' @param project_dir character path to the project root directory
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_export = function(project_dir, write_csv=TRUE) {

  cat('\n')
  message('resampling')

  # chosen to match SWAT inputs: each of these has a specifically chosen aggregation function
  tmp_max_path = project_dir |> nc_export('tmp', write_csv=write_csv, fun='max', tz='MST')
  tmp_min_path = project_dir |> nc_export('tmp', write_csv=write_csv, fun='min', tz='MST')
  hum_mean_path = project_dir |> nc_export('hum', write_csv=write_csv, fun='mean', tz='MST')
  pcp_mean_path = project_dir |> nc_export('pcp', write_csv=write_csv, fun='mean', tz='MST')
  wnd_mean_path = project_dir |> nc_export('wnd', write_csv=write_csv, fun='mean', tz='MST')
}

