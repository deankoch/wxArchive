#' List all times and variables available locally in the RAP/RUC/GFS archive
#'
#' This checks the sub-directories of `project_dir` (both 'rap' and 'gfs') for
#' NetCDF data on the variables named in `.nm_output_var`, and reports some stats
#' to the console
#'
#' @param project_dir character path to the project root directory
#'
#' @return a list of character vectors, paths to individual NetCDF files
#' @export
workflow_list = function(project_dir) {

  # base directories for all NetCDF and GRIB files from RAP/RUC and GFS
  base_dir_rap = project_dir |> file.path('rap')
  base_dir_gfs = project_dir |> file.path('gfs')

  # make a list of all datasets with preference for RAP archive over GFS
  rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
  gfs_nc_path = file_wx('nc', base_dir_gfs, .nm_resample, as.list(.nm_gfs_var))
  p_all = Map(\(rap, gfs) c(rap, gfs), rap = rap_nc_path, gfs = gfs_nc_path)

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

    # print stats to console
    paste0('\n', nm_fix_wid[[i]], ' | ',
           n_file, ' files | ',
           n, ' times | ',
           time_min, ' to ', time_max) |> cat()

    # print file paths to console
    c('\n', rep('-', nm_len-1)) |> paste(collapse='') |> cat()
    paste0('\n > ', p) |> cat()
    cat('\n')
  }

  cat('\n')
  return(invisible(p_all))
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
#' @return Returns nothing but possible writes to `project_dir`
#' @export
workflow_update_rap = function(project_dir, from=NULL, to=NULL) {

  # path to area of interest polygon
  aoi_path = project_dir |> file.path('aoi.geojson')
  aoi = sf::st_read(aoi_path)

  # all the output files go here
  base_dir_rap = project_dir |> file.path('rap')

  ## RAP/RUC

  # download RUC/RAP files
  message('\nupdating RAP/RUC GRIB archive')
  archive_update(base_dir = base_dir_rap,
                 hour_rel = .hour_rel_rap,
                 from = from,
                 to = to,
                 model = 'rap_archive') |> invisible()

  # export to nc
  message('\nupdating NetCDF files')
  nc_update(aoi = aoi,
            base_dir = base_dir_rap,
            output_nm = .nm_src_rap,
            regex = .rap_regex) |> invisible()

  # copy precip from components (applies to early years)
  message('\nprocessing precipitation layers')
  pcp_update(base_dir = base_dir_rap,
             pcp_nm = .var_pcp,
             input_nm = .nm_src_rap,
             output_nm = .nm_rap) |> invisible()

  # resampled coarse to fine
  message('\nresampling')
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
#' @return Returns nothing but possible writes to `project_dir`
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
#' from the time series named in in `.nm_resample_rap`.
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
#' @return Returns nothing but possible writes to `project_dir`
#' @export
workflow_impute_rap = function(project_dir) {

  # all the output files go here
  base_dir_rap = project_dir |> file.path('rap')

  # fill gaps
  message('\nimputing missing times')
  time_impute(var_nm = .nm_output_var,
              base_dir = base_dir_rap,
              input_nm = .nm_resample_rap,
              output_nm = .nm_complete) |> invisible()
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
#' @return Returns nothing but possible writes to `project_dir`
#' @export
workflow_update_gfs = function(project_dir) {

  # path to area of interest polygon
  aoi_path = project_dir |> file.path('aoi.geojson')
  aoi = sf::st_read(aoi_path)

  # all the output files go here
  base_dir_gfs = project_dir |> file.path('gfs')

  ## GFS

  # RAP archive end time determines start time of the GFS output
  base_dir_rap = project_dir |> file.path('rap')
  rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
  rap_time = lapply(rap_nc_path, \(p) time_wx(p))
  from = do.call(c, lapply(rap_time, \(x) max(x[['time_obs']]) ))

  # download GFS files
  message('\nupdating GFS GRIB archive')
  gfs_result = archive_update(base_dir = base_dir_gfs,
                              hour_pred = .hour_pred_gfs,
                              hour_rel = .hour_rel_gfs,
                              aoi = aoi,
                              model = 'gfs_0p25',
                              alternate = FALSE)

  # delete the old GFS NetCDF directories
  base_dir_gfs |> file.path('coarse') |> unlink(recursive=TRUE)
  base_dir_gfs |> file.path(.nm_resample) |> unlink(recursive=TRUE)

  # part 8: export latest GFS data to nc (creates "coarse" subdirectory)
  nc_update(aoi = aoi,
            base_dir = base_dir_gfs,
            output_nm = list(coarse=.nm_gfs),
            regex = .gfs_regex,
            from = from) |> invisible()

  # load an example grid at fine resolution (second rast call drops cell values)
  r_fine = file_wx('nc', base_dir_rap, .nm_spatial, names(.rap_regex)[[1]]) |>
    terra::rast() |> terra::rast()

  # part 9: resample to match RAP grid
  nc_resample(var_nm = .nm_gfs_var,
              base_dir = base_dir_gfs,
              input_nm = list(coarse=.nm_gfs),
              output_nm = .nm_resample,
              r_fine = r_fine) |> invisible()

}


