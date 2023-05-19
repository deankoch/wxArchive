#' Export completed time series and (optionally) aggregate to daily
#'
#' This exports all data each of the variables named in `var_nm` to their own (single)
#' NetCDF file in `output_nm`. `var_nm` should be a subset of the names returned by
#' `nc_list(base_dir)`.
#'
#' If `write_csv=TRUE`, the function additionally writes, for each variable, a CSV file
#' containing the matrix of data associated with each variable, and a time column with a
#' character string giving the date/time in UTC. These files have the same file-names as
#' the NetCDF versions (up to file extension).
#'
#' If `fun=NA`, the function does no aggregation and writes all time points.
#'
#' If `fun` is not `NA` then the last three arguments are passed to `nc_aggregate` to control
#' the alignment of the aggregation window. `fun` specifies the function to use for combining
#' times within the window (see `?nc_aggregate`).  Set `tz` to the desired output time
#' zone, and leave `origin_hour=0` to have each day begin at 12AM (in time zone `tz`).
#'
#' File names for aggregate data are given the suffix `_daily_<fun>` - eg with `fun='mean'`,
#' the output "tmp.nc" becomes "tmp_daily_mean.nc".
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
                     output_nm = .nm_export,
                     write_csv = FALSE,
                     fun = mean,
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

  # get the name of the supplied function
  is_agg = !is.null(fun)
  fun_nm = fun |> substitute() |> deparse()

  # define output files
  output_nc = file_wx('nc', base_dir, output_nm, as.list(var_nm), make_dir=TRUE)
  if( is_agg ) output_nc = gsub('.nc$',  paste0('_daily_', fun_nm, '.nc'), output_nc)
  output_csv = gsub('.nc$', '.csv', output_nc)
  output_geojson = base_dir |> file.path(output_nm, 'grid_points.geojson')

  # loop over files
  for( i in seq_along(p_fetch) ) {

    cat('\nprocessing', names(p_fetch)[[i]], '...')
    t1 = proc.time()

    # make a SpatRaster of output layers
    r_i = if( is_agg ) {

      # aggregate to daily
      r_out = p_fetch[[i]] |> nc_aggregate(fun=fun, tz=tz, origin_hour=origin_hour) |> terra::rast()
      terra::time(r_out) = as.Date(names(r_out))
      r_out

    } else {

      # copy all layers unchanged
      time_i = time_wx(p_fetch[[i]])[['time_obs']]
      p_fetch[[i]] |> nc_layers(times=time_i, na_rm=TRUE)
    }

    # create the output nc file and write to it
    r_i |> nc_write(output_nc[[i]])

    if( write_csv ) {

      # export to data frame
      cat('writing data matrix to', output_csv[i])
      cell_id = paste0('grid_', seq( terra::ncell(r_i)))
      df_i = terra::time(r_i) |>
        as.character(tz=tz) |>
        data.frame() |>
        cbind( t(r_i[]) ) |>
        stats::setNames( nm=c('time', cell_id) )

      # write the CSV to disk
      df_i |> write.csv(output_csv[i], row.names=FALSE)

      # grid point locations as st POINT collection in WGS84 coordinates
      pts = terra::crds(r_i) |>
        sf::st_multipoint() |>
        sf::st_sfc(crs=terra::crs(r_i)) |>
        sf::st_cast('POINT') |>
        sf::st_transform(4326)

      # add key and write to disk as geojson
      unlink(output_geojson)
      sf::st_sf(data.frame(id=cell_id), geometry=pts) |>
        sf::st_write(output_geojson) |>
        suppressWarnings() # GDAL generates spurious warnings on linux
    }

    # remove the large source raster from memory
    rm(r_i)
    gc()
    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  }

  # report all files written
  output_paths = output_nc
  if(write_csv) output_paths = c(output_paths, output_csv, output_geojson)
  return(output_paths)
}



#' Aggregate sub-daily data to daily
#'
#' Returns a list of SpatRasters, one for each day in the time series
#' `p`. Times within a given day are aggregated using the function
#' `fun` should be a function (not a string naming it), such as `mean`
#' or `max` (see `?terra::app`).
#'
#' This has not been tested on incomplete data
#'
#' @param p character vector path to the nc file(s)
#' @param fun function, one of "mean", "min", or "max"
#' @param tz character time zone for `origin_hour`
#' @param origin_hour integer, steps are aligned to start at this hour of the day
#'
#' @return SpatRaster, the aggregated data
#' @export
nc_aggregate = function(p, fun=mean, tz='UTC', origin_hour=0L) {

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

  # compute stats in loop over days then remove the large source raster from memory
  cat('\ncomputing', fun, 'of', n_per, 'steps on', n_out, 'day(s)')
  r_result = lapply(list_idx, \(j) terra::app(r[[j]], fun) )
  rm(r)
  gc()

  names(r_result) = date_out
  return(r_result)
}
