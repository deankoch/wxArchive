#' Aggregate the completed 2-hourly time series to daily
#'
#' This creates a (single) NetCDF file in `output_nm` for each of the variables named in
#' `var_nm` by applying the function named in `fun` to aggregate by day. `var_nm` should
#' be a subset of the names returned by `nc_list(base_dir)`.
#'
#' If `fun=NULL`, the function does no aggregation and writes all times to the output file.
#'
#' If `fun` is not `NULL` then the last three arguments are passed to `.nc_aggregate_time` to
#' control the alignment of the aggregation window. `fun` specifies the function to use for
#' combining times within the window (see `?.nc_aggregate_time`).  Set `tz` to the desired
#' output time zone, and leave `origin_hour=0` to have each day begin at 12AM (in time zone
#' `tz`).
#'
#' Use `from` and `to` to specify a date range to update (inclusive), or leave them `NULL`
#' to use a default range. The default range is meant to all layers originating from GFS
#' (allowing them to be can be replaced by newly added RAP layers, or more recently released
#' GFS forecasts).
#'
#' The default for `to` is always the latest available date in the input. The default for
#' `from` is 10 days before the latest date found in the existing output files. If there
#' are no existing outputs, the default is set to the earliest available date in the input.
#'
#' File names for aggregate data are given the suffix `_daily_<fun>` - eg with `fun='mean'`,
#' the output "tmp.nc" becomes "tmp_daily_mean.nc".
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param var_nm character vector, the variable(s) to export (NULL to export all)
#' @param output_nm character, the sub-directory name for output
#' @param fun function, a function name like "mean", "min", or "max" (see `?.nc_aggregate_time`)
#' @param tz character time zone for `origin_hour`
#' @param origin_hour integer, steps are aligned to start at this hour of the day
#' @param from POSIXct, the start of the time range to process
#' @param to POSIXct, the end of the time range to process
#'
#' @return vector of file paths written to the `output_nm` directory
#' @export
nc_aggregate_time = function(base_dir,
                             var_nm = NULL,
                             output_nm = .nm_daily,
                             fun = 'mean',
                             tz = 'UTC',
                             origin_hour = 0L,
                             from = NULL,
                             to = NULL) {

  # get list of all relevant input netCDF files
  p_all = nc_list(base_dir)
  if( is.null(var_nm) ) var_nm = names(p_all)

  # validity check for requested variable names
  is_valid = var_nm %in% names(p_all)
  msg_invalid = paste('variable name(s)', paste(var_nm[!is_valid], collapse=', '), 'not recognized')
  if( any(!is_valid) ) stop(msg_invalid)
  p_fetch = p_all[var_nm]

  # define output directory
  output_nc = file_wx('nc', base_dir, output_nm, as.list(var_nm), make_dir=TRUE)
  is_agg = !is.null(fun)
  if( is_agg ) output_nc = gsub('.nc$',  paste0('_', fun, '.nc'), output_nc)

  # copy existing output dates
  date_existing = time_wx( nc_chunk(output_nc) )[['time_obs']]
  is_initial = length(date_existing) == 0

  # loop over files
  for( i in seq_along(p_fetch) ) {

    cat('\nprocessing', names(p_fetch)[[i]], '...')
    t1 = proc.time()

    # find available times from RAP/GFS combined
    time_i = time_wx(p_fetch[[i]])[['time_obs']]
    if( length(time_i) == 0 ) {

      cat('nothing to write\n')
      next
    }

    # set default starting/ending times
    if( is.null(to) ) to = max(time_i)
    if( is.null(from) ) {

      # on first call this writes everything
      from = min(time_i)

      # subsequently the default start time is 10 days before latest time
      if(!is_initial) from = as.POSIXct(max(date_existing)) - ( 10 * (60*60*24) )
    }

    # silently fix invalid start/end times
    from = as.POSIXct(from)
    to = as.POSIXct(to)
    if( from < min(time_i) ) from = min(time_i)
    if( to > max(time_i) ) to = max(time_i)

    # filter to requested range
    time_i = time_i[ ( time_i >= from ) & ( time_i <= to ) ]
    if( length(time_i) == 0 ) {

      cat('nothing to write\n')
      next
    }

    # split by year
    yr_i = time_i |> format('%Y', tz=tz)
    yr_unique = yr_i |> unique()
    cat('\nprocessing', length(yr_unique), 'year(s)...\n')
    for(yr in yr_unique) {

      cat('\nyear', yr, paste0('(', tz, ')'))
      time_i_yr = time_i[yr_i == yr]

      # make a SpatRaster of output layers
      r_i = if( is_agg ) {

        # aggregate to daily (NULL if failed)
        aggregate_result = p_fetch[[i]] |> .nc_aggregate_time(times = time_i_yr,
                                                              fun = fun,
                                                              tz = tz,
                                                              origin_hour = origin_hour)

        # skip when there are not enough layers to make a single day
        if( is.null(aggregate_result) ) {

          cat('\nskipped (not enough layers)\n')
          next
        }

        # else set times and return SpatRaster
        r_out = aggregate_result |> terra::rast()
        terra::time(r_out) = as.Date(names(r_out))
        r_out

      } else {

        # copy all layers unchanged
        p_fetch[[i]] |> nc_layers(times=time_i_yr, na_rm=TRUE)
      }

      # create the output nc file and write to it
      r_i |> nc_write_chunk(output_nc[[i]], insert=TRUE)

      # remove the large source raster from memory
      rm(r_i)
      gc()
    }

    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  }

  # report the file paths of the times series that was modified
  return( nc_chunk(output_nc) )
}



#' Aggregate sub-daily data to daily
#'
#' Helper function for `nc_aggregate_time`
#'
#' Returns a list of SpatRasters, one for each day in the time series `p`. Times
#' within a given day are aggregated using the function `fun` should be a string
#' naming the function, like `"mean"` (and not `mean` itself) and this function
#' must be compatible with `terra::app`.
#'
#' The function clips the data to the nearest `origin_hour` in the time zone `tz`,
#' omitting times from incomplete days that appear at the beginning and/or end of
#' `times`. If this results in a time series of length-0, the function returns NULL
#'
#' @param p character vector path to the nc file(s)
#' @param times vector of POSIXct times to include (default is all)
#' @param fun character naming a function, such as "mean", "min", or "max"
#' @param tz character time zone for `origin_hour`
#' @param origin_hour integer, steps are aligned to start at this hour of the day
#'
#' @return SpatRaster, the aggregated data or NULL
#' @export
.nc_aggregate_time = function(p, times=NULL, fun='mean', tz='UTC', origin_hour=0L) {

  # check all available times
  cat('\nchecking times in', length(p), 'file(s)')
  p_time = time_wx(p)
  if( is.null(times) ) times = p_time[['time_obs']]
  times = times[times %in% p_time[['time_obs']]]
  if(length(times) < 2) return(NULL)

  # lead all requested data
  cat('\ncopying', length(times), 'layers to RAM')
  r = nc_layers(p, times=times, na_rm=TRUE)

  # check for gaps and detect step size (hours)
  ts_df = data.frame(posix_pred=terra::time(r)) |> archive_pad()
  step_data = time_step(ts_df)
  n_per = 24L / step_data
  if(length(times) < n_per) return(NULL)

  # starting hour of the data in target time zone
  start_hour = ts_df[['posix_pred']] |> min() |> format('%H', tz=tz) |> as.integer()
  start_idx = ( step_data + (origin_hour - start_hour) %% 24 ) / step_data
  start_date = ts_df[['posix_pred']][start_idx] |> as.Date(tz=tz)

  # list of indices to aggregate in each step and the corresponding dates
  n_out = floor( ( length(ts_df[['posix_pred']]) - start_idx + 1 ) / n_per )
  list_idx = seq(n_out) |> lapply(\(i) seq(n_per) + (i-1)*n_per )
  date_out = seq.Date(start_date, by='day', length.out=n_out)

  # compute stats in loop over days then remove the large source raster from memory
  cat('\ncomputing', fun, 'of', n_per, 'steps on', n_out, 'day(s)')
  r_result = lapply(list_idx, \(j) terra::app(r[[j]], get(fun)) )
  rm(r)
  gc()

  names(r_result) = date_out
  return(r_result)
}
