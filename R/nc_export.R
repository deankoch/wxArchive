#' Export completed time series and (optionally) aggregate to daily
#'
#' The last two arguments are passed to `nc_aggregate` to control the alignment of
#' the aggregation window. `fun` specifies the function (one of "mean", "min", or "max")
#' to use for combining times within the window. If `fun=NA`, the function does no
#' aggregation and writes all time points.
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param var_nm character vector, the variable(s) to export (NULL to export all)
#' @param output_nm character, the sub-directory name for output
#' @param fun function, one of "mean", "min", or "max"
#' @param tz character time zone for `origin_hour`
#' @param origin_hour integer, steps are aligned to start at this hour of the day
#' @param write_csv logical, Writes a CSV copy of the data, and point locations to geojson
#'
#' @return vector of file paths written to the `output_nm` directory
#' @export
nc_export = function(base_dir,
                     var_nm = NULL,
                     output_nm = 'export',
                     write_csv = FALSE,
                     fun = 'mean',
                     tz = 'UTC',
                     origin_hour = 0L) {

  # get list of all relevant input netCDF files
  p_all = nc_list(base_dir)
  if( is.null(var_nm) ) var_nm = names(p_all)

  # validity check for requested variable names
  is_valid = var_nm %in% names(p_all)
  msg_invalid = paste('variable name(s)', paste(var_nm[!is_valid], collapse=', '), 'not recognized')
  if( any(!is_valid) ) stop(msg_invalid)
  p_fetch = p_all[var_nm]

  # define output files
  is_agg = !is.na(fun)
  output_nc = file_wx('nc', base_dir, output_nm, as.list(var_nm), make_dir=TRUE)
  if( is_agg ) output_nc = gsub('.nc$',  paste0('_daily_', fun, '.nc'), output_nc)
  output_csv = gsub('.nc$', '.csv', output_nc)
  output_geojson = base_dir |> file.path(output_nm, 'grid_points.geojson')

  # define output files and create the directory then loop over files
  for( i in seq_along(p_fetch) ) {

    cat('\n\nprocessing', names(p_fetch)[[i]], '...')
    t1 = proc.time()

    # make a SpatRaster of output layers
    r_i = if( is_agg ) {

      # aggregate to daily
      p_fetch[[i]] |> nc_aggregate(fun=fun, tz=tz, origin_hour=origin_hour)

    } else {

      # copy all layers unchanged
      time_i = time_wx(p_fetch[[i]])[['time_obs']]
      p_fetch[[i]] |> nc_layers(times=time_i, na_rm=TRUE)
    }

    # create the output nc file and write to it
    r_i |> nc_write(output_nc[[i]])
    if( write_csv ) {

      # export to data frame
      cat('\nwriting data matrix to', output_csv[i])
      cell_id = paste0('grid_', seq( terra::ncell(r_i)))
      df_i = terra::time(r_i) |>
        as.character(tz=tz) |>
        data.frame() |>
        cbind( t(r_i[]) ) |>
        stats::setNames( nm=c('time', cell_id) )

      # write the CSV to disk
      df_i |> write.csv(output_csv[i], row.names=FALSE)
    }

    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  }

  # grid point locations as st POINT collection in WGS84 coordinates
  pts = terra::crds(r_i) |>
    sf::st_multipoint() |>
    sf::st_sfc(crs=terra::crs(r_i)) |>
    sf::st_cast('POINT') |>
    sf::st_transform(4326L)

  # add key and write to disk as geojson
  unlink(output_geojson)
  sf::st_sf(data.frame(id=cell_id), geometry=pts) |> sf::st_write(output_geojson)

  # report all files written
  output_paths = output_nc
  if(write_csv) output_paths = c(output_paths, output_csv, output_geojson)
  return(output_paths)
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
  cat('\nchecking times in', length(p), 'file(s)')
  p_time = time_wx(p)
  cat('\ncopying', length(p_time[['time_obs']]), 'layers to RAM')
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
  cat('\ncomputing', fun, 'of', n_per, 'steps on', n_out, 'day(s)')
  if( fun == 'mean' ) r_result = do.call(c, lapply(list_idx, \(j) terra::mean(r[[j]])) )
  if( fun == 'min' ) r_result = do.call(c, lapply(list_idx, \(j) terra::min(r[[j]])) )
  if( fun == 'max' ) r_result = do.call(c, lapply(list_idx, \(j) terra::max(r[[j]])) )
  terra::time(r_result) = date_out
  return(r_result)
}
