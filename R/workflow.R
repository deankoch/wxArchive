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
#' @param daily logical shows daily data first
#' @param tz character time zone string for displaying times
#'
#' @return a list of character vectors, paths to individual NetCDF files
#' @export
workflow_list = function(project_dir, daily=TRUE, quiet=FALSE, tz='MST') {

  # make a list of all datasets with preference for RAP archive over GFS
  p_all = project_dir |> nc_list()
  if( quiet ) return(p_all)

  # AGGREGATE DATA

  # report on daily variables
  if( daily ) {

    # base directories for all daily files
    daily_path = file_wx('nc', project_dir, .nm_daily, .var_daily)

    # fixed width names for printout
    nm_all = .var_daily_pairs |> sapply(\(x) paste(x['var'], '|', x['fun']))
    nm_len = nm_all |> sapply(nchar) |> max()
    nm_fix_wid = lapply(nm_all, \(nm) paste0(nm, paste(rep(' ', nm_len - nchar(nm)), collapse='')))

    # same for paths
    p_len = daily_path |> sapply(nchar) |> max()
    p_fix_wid = lapply(daily_path, \(nm) paste0(nm, paste(rep(' ', p_len - nchar(nm)), collapse='')))

    # loop over daily variables
    n_daily = length(daily_path)
    message('\n', n_daily, ' daily variable(s)')
    for(i in seq(n_daily)) {

      # this returns empty list when the file is missing
      p = daily_path[i]
      time_all = p |> time_wx()
      if( length(time_all) == 0 ) {

        paste0('\n', p_fix_wid[i], ' | ', nm_fix_wid[i], ' | not found') |> cat()
        next
      }

      # count times and find their range
      time_min = time_all[['time_obs']] |> min() |> as.character(tz=tz)
      time_max = time_all[['time_obs']] |> max() |> as.character(tz=tz)
      n = time_all[['time_obs']] |> length()

      # check if time series is complete
      ts_df = data.frame(posix_pred=time_all[['time_obs']]) |> archive_pad(quiet=TRUE)
      n_miss = ts_df[['ts_hours']] |> is.na() |> sum()
      msg_miss = ifelse(n_miss==0, ' (complete)', paste0(' (', n_miss, ' missing)'))

      # print stats and path to console
      paste0('\n', p_fix_wid[i], ' | ',
             nm_fix_wid[i], ' | ',
             n, ' day(s) | ',
             time_min, ' to ', time_max,
             msg_miss) |> cat()
    }
    cat('\n')
  }


  # HOURLY DATA

  # fixed width names for printout
  nm_all = names(p_all)
  nm_len = nm_all |> sapply(nchar) |> max()
  nm_fix_wid = lapply(nm_all, \(nm) paste0(nm, paste(rep(' ', nm_len - nchar(nm)), collapse='')))

  # check available times in a loop, printing as we go
  message('\n', length(p_all), ' two-hourly variable(s)')
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
            from = from,
            to = to,
            output_nm = .nm_src_rap,
            regex = .rap_regex) |> invisible()

  # copy precip from components (applies mostly to early years)
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

#' Fit a temporal model to the RUC/RAP times series or load it from a file
#'
#' A wrapper for `time_fit`. This reads its input from the sub-directories
#' named in `.nm_resample_rap` and writes its output to `.nm_model` (see
#' `?time_fit`).
#'
#' Alternatively, load a previously fitted model from a zip archive by providing
#' the file name in `from_file`. This zip should contain the contents of
#' `.nm_model` (including the "temporal" subdirectory) generated from a `time_fit`
#' call on a similar dataset (in terms of time period and AOI)
#'
#' Call this function at least once after running `workflow_update_rap` for
#' the first time (and before running `workflow_impute_rap`)
#'
#' @param project_dir character path to the project root directory
#' @param from_file character file name of zip containing previously fitted model
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_fit_temporal = function(project_dir, from_file='model.zip') {

  # all the output files go here, in sub-directory .nm_model
  base_dir_rap = project_dir |> file.path('rap')

  # fit the models from scratch
  if( is.null(from_file) ) {

    # part 6: fit temporal model to fine grid (include all layers)
    time_fit(var_nm = .nm_output_var,
             base_dir = base_dir_rap,
             model_nm = .nm_model,
             input_nm = .nm_resample_rap) |> invisible()

  } else {

    # load from archive (overwrites existing model)
    zip_path = project_dir |> file.path(from_file)
    dest_path = base_dir_rap |> file.path(.nm_model)
    if( !file.exists(zip_path) ) stop('file ', zip_path, ' not found')
    zip_path |> unzip(exdir=dest_path)
  }
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
              model_nm = .nm_model,
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
  wnd_update(base_dir = base_dir_rap,
             wnd_nm = .var_wnd,
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
#' @param n_ahead integer, predictions run ahead of GFS data by this many days
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_update_gfs = function(project_dir, n_ahead=3) {

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
  base_dir_gfs |> file.path(.nm_complete) |> unlink(recursive=TRUE)

  # export latest GFS data to nc (creates "coarse" subdirectory)
  cat('\n')
  message('updating NetCDF files')
  nc_update(aoi = aoi,
            base_dir = base_dir_gfs,
            output_nm = list(coarse=.nm_gfs),
            regex = .gfs_regex,
            from = from) |> invisible()

  # load an example grid at fine resolution (second rast call drops cell values)
  r_fine = nc_chunk(file_wx('nc', base_dir_rap, 'fine', names(.rap_regex)[[1]]))[1] |>
    terra::rast() |> terra::rast()

  # resample to match RAP grid
  cat('\n')
  message('resampling')
  nc_resample(var_nm = .nm_gfs_var,
              base_dir = base_dir_gfs,
              input_nm = list(coarse=.nm_gfs),
              output_nm = .nm_resample,
              r_fine = r_fine) |> invisible()

  # predictions will run n_ahead days past the end of the GFS data
  gfs_nc_path = file_wx('nc', base_dir_gfs, .nm_resample, .nm_gfs_var)
  gfs_time = lapply(gfs_nc_path, \(p) time_wx(p))
  gfs_end = do.call(c, lapply(gfs_time, \(x) max(x[['time_obs']]) )) |> max()
  until = gfs_end + ( n_ahead * 24 * 60 * 60 )

  # prediction using model results from RAP analysis
  cat('\n')
  message('extending forecasts by ', n_ahead, ' day(s)')
  time_impute(var_nm = .nm_gfs_var,
              base_dir = base_dir_gfs,
              to = until,
              model_dir = base_dir_rap,
              model_nm = .nm_model,
              input_nm = .nm_resample,
              output_nm = .nm_complete) |> invisible()

  # create wind speed layers
  cat('\n')
  message('computing wind speed from u/v components')
  base_dir_gfs |> wnd_update(wnd_nm = .var_wnd,
                             uv_nm = .var_wnd_uv,
                             input_nm = .nm_complete_gfs) |> invisible()
}

#' Exported completed time series to daily aggregate values and write to disk
#'
#' Wrapper for `nc_aggregate_time`
#'
#' This prepares five output variables aggregated to daily average, or maximum or minimum.
#' The pairing of variable names and aggregation functions is set up in the global
#' constant `.var_daily_pairs` and the name of the output folder is `.nm_daily`.
#'
#' By default the function overwrites layers starting from 10 days before the latest
#' date in the existing output files. This ensures that any/all GFS forecasts are
#' updated with fresh estimates, or else replaced by the more precise RAP archive
#' estimates. Specify a different time range to process using `from` and `to`.
#'
#' The time zone string `tz` controls the alignment of the output frames, so that each
#' "day" starts (and ends) on the hour 12AM in the time zone `tz`. Set this to the local
#' time zone for your AOI.
#'
#' @param project_dir character path to the project root directory
#' @param from POSIXct start of time range (default is earliest available)
#' @param to POSIXct end of time range (default is latest available)
#'
#' @return returns nothing but possible writes to `project_dir`
#' @export
workflow_daily = function(project_dir, from=NULL, to=NULL, tz='MST') {

  cat('\n')
  message('aggregating data by day')

  # export all in a loop
  export_paths = .var_daily_pairs |> lapply(\(x) nc_aggregate_time(base_dir = project_dir,
                                                                   var_nm = x['var'],
                                                                   output_nm = .nm_daily,
                                                                   fun = x['fun'],
                                                                   tz = tz,
                                                                   from = from,
                                                                   to = to))
}


#' Construct down-scaled (fine resolution) estimates of daily spatial grid data
#'
#' Wrapper for `nc_downscale` to produced down-scaled versions of the daily outputs
#' produced by `workflow_daily`.
#'
#' This requires a digital elevation model (DEM) with data in meters. This is expected
#' in the GeoTIFF file "elev_m.tif" in `project_dir`, but you can specify an alternate
#' name or location with `dem_path` (passed to `terra::rast`).
#'
#' The DEM is warped to the target coordinate reference system by bilinear averaging
#' Output is written to the sub-directory `.nm_export` of `project_dir`. By default the
#' function overwrites layers starting from 10 days before the latest date in the existing
#' output files (see `?workflow_daily`). Change this by setting start/end dates
#' in `from`/`to`.
#'
#'
#'
workflow_downscale = function(project_dir, down=100,
                              dem_path=NULL, poly_path=NULL,
                              from=NULL, to=NULL) {

  cat('\n')
  message('downscaling')

  # check for polygons to crop output
  if( is.null(poly_path) ) poly_path = project_dir |> file.path('export.geojson')
  if( !file.exists(poly_path) ) {

    warning('"export.geojson" not found. Defaulting to bounding box of AOI polygon')
    poly_path = project_dir |> file.path('aoi.geojson')
    if( !file.exists(poly_path) ) stop('"aoi.geojson" not found in ', project_dir)
  }

  # load the polygons and DEM
  poly_in = poly_path |> sf::st_read()
  if( is.null(dem_path) ) dem_path = project_dir |> file.path('elev_m.tif')
  dem = terra::rast(dem_path)

  # run the workflow
  nc_downscale(base_dir = project_dir,
               dem = dem,
               down = down,
               input_nm = .nm_daily,
               model_nm = .nm_model,
               output_nm = .nm_down,
               var_nm = .var_daily,
               poly_out = poly_in,
               edge_buffer = NULL,
               from = from,
               to = to,
               write_nc = TRUE)
}
