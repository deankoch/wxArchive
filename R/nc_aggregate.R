#' Export completed time series, optionally aggregating to daily
#'
#' @param var_nm character vector, the variable(s) to export
#' @param agg list with entries "fun", "tz", "origin_hour" to pass to `nc_aggregate`
#'
#' @return returns nothing but writes to the `output_nm` directory
#' @export
nc_export = function(base_dir, p_all=NULL, output_nm=fun, agg=NULL) {

  # get list of all relevant input netCDF files
  if(is.null(p_all)) p_all = workflow_list(base_dir, quiet=TRUE)

  # define output files and create the directory then loop over files
  output_nc = file_wx('nc', base_dir, output_nm, as.list(names(p_all)), make_dir=TRUE)
  for( i in seq_along(p_all) ) {

    cat('\n\nprocessing', names(p_all)[[i]], '...')
    t1 = proc.time()
    r_i = if( is.null(agg) ) {

      # copy all layers unchanged
      time_i = time_wx(p_all[[i]])[['time_obs']]
      p_all[[i]] |> nc_layers(times=time_i)

    } else {

      # aggregate to daily
      p_all[[i]] |> nc_aggregate(fun = agg[['fun']],
                                 tz = agg[['tz']],
                                 origin_hour = agg[['origin_hour']])
    }

    # append results to existing data file (or create the file and write to it)
    r_i |> nc_write(output_nc[[nm]])

    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  }
}






#' Aggregate sub-daily data to daily
#'
#' @param p character vector path to the nc file(s)
#' @param fun function, one of "mean", "min", or "max"
#' @param tz character time zone for `origin_hour`
#' @param origin_hour integer, steps are aligned to start at this hour of the day
#'
#' @return SpatRaster, the aggregated data
#' @export
nc_aggregate = function(p, fun='mean', tz='UTC', origin_hour=0L) {

  # load everything into RAM (~5GB)
  p_time = time_wx(p)
  cat('\nloading', length(p_time[['time_obs']]), 'layers into RAM')
  r = nc_layers(p, times=p_time[['time_obs']], na_rm=TRUE)

  # check for gaps and detect step size (hours)
  ts_df = data.frame(posix_pred=terra::time(r)) |> archive_pad()
  step_data = time_step(ts_df)

  # starting hour of the data in target time zone
  start_hour = ts_df[['posix_pred']] |> min() |> format('%H', tz=tz) |> as.integer()
  start_idx = 1 + ( ( (origin_hour - start_hour) %% 24 ) / step_data )
  start_date = ts_df[['posix_pred']][start_idx] |> as.Date(tz=tz)

  # list of indices to aggregate in each step and the corresponding dates
  n_per = 24L / step_data
  n_out = floor( ( length(ts_df[['posix_pred']]) - start_idx + 1 ) / n_per )
  list_idx = seq(n_out) |> lapply(\(i) seq(n_per) + (i-1)*n_per )
  date_out = seq.Date(start_date, by='day', length.out=n_out)

  # compute stats in loop over days
  cat('\ncomputing', fun, 'of', n_per, 'steps for', n_out, 'day(s)')
  if( fun == 'mean' ) r_result = do.call(c, lapply(list_idx, \(j) terra::mean(r[[j]])) )
  if( fun == 'min' ) r_result = do.call(c, lapply(list_idx, \(j) terra::min(r[[j]])) )
  if( fun == 'max' ) r_result = do.call(c, lapply(list_idx, \(j) terra::max(r[[j]])) )
  terra::time(r_result) = date_out
  return(r_result)
}
