#' Downscale and export completed time series by aggregating over areas
#'
#' This exports all data each of the variables named in `var_nm` to their own (single)
#' NetCDF file in `output_nm`. `var_nm` should be a subset of the names returned by
#' `nc_list(base_dir)`.
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param var_nm character vector, the variable(s) to export (NULL to export all)
#' @param output_nm character, the sub-directory name for output
#' @param fun character, a function name like "mean", "min", or "max"
#' @param write_csv logical, Writes a CSV copy of the data, and point locations to geojson
#'
#' @return vector of file paths written to the `output_nm` directory
#' @export
nc_export = function(base_dir,
                     var_nm = NULL,
                     output_nm = .nm_daily,
                     write_csv = FALSE,
                     fun = 'mean') {

  stop('not implemented yet')

  # # get list of all relevant input netCDF files
  # p_all = nc_list(base_dir)
  # if( is.null(var_nm) ) var_nm = names(p_all)
  #
  # # validity check for requested variable names
  # is_valid = var_nm %in% names(p_all)
  # msg_invalid = paste('variable name(s)', paste(var_nm[!is_valid], collapse=', '), 'not recognized')
  # if( any(!is_valid) ) stop(msg_invalid)
  # p_fetch = p_all[var_nm]
  #
  # # define output files
  # output_nc = file_wx('nc', base_dir, output_nm, as.list(var_nm), make_dir=TRUE)
  # is_agg = !is.null(fun)
  # if( is_agg ) output_nc = gsub('.nc$',  paste0('_', fun, '.nc'), output_nc)
  # output_csv = gsub('.nc$', '.csv', output_nc)
  # output_geojson = base_dir |> file.path(output_nm, 'grid_points.geojson')
  #
  # # loop over files
  # for( i in seq_along(p_fetch) ) {
  #
  #   cat('\nprocessing', names(p_fetch)[[i]], '...')
  #   t1 = proc.time()
  #
  #   # make a SpatRaster of output layers
  #   r_i = if( is_agg ) {
  #
  #     # aggregate to daily
  #     r_out = p_fetch[[i]] |> .nc_aggregate_time(fun=fun, tz=tz, origin_hour=origin_hour) |> terra::rast()
  #     terra::time(r_out) = as.Date(names(r_out))
  #     r_out
  #
  #   } else {
  #
  #     # copy all layers unchanged
  #     time_i = time_wx(p_fetch[[i]])[['time_obs']]
  #     p_fetch[[i]] |> nc_layers(times=time_i, na_rm=TRUE)
  #   }
  #
  #   # create the output nc file and write to it
  #   r_i |> nc_write(output_nc[[i]])
  #
  #   if( write_csv ) {
  #
  #     # export to data frame
  #     cat('writing data matrix to', output_csv[i])
  #     cell_id = paste0('grid_', seq( terra::ncell(r_i)))
  #     df_i = terra::time(r_i) |>
  #       as.character(tz=tz) |>
  #       data.frame() |>
  #       cbind( t(r_i[]) ) |>
  #       stats::setNames( nm=c('time', cell_id) )
  #
  #     # write the CSV to disk
  #     df_i |> write.csv(output_csv[i], row.names=FALSE)
  #
  #     # grid point locations as st POINT collection in WGS84 coordinates
  #     pts = terra::crds(r_i) |>
  #       sf::st_multipoint() |>
  #       sf::st_sfc(crs=terra::crs(r_i)) |>
  #       sf::st_cast('POINT') |>
  #       sf::st_transform(4326)
  #
  #     # add key and write to disk as geojson
  #     unlink(output_geojson)
  #     sf::st_sf(data.frame(id=cell_id), geometry=pts) |>
  #       sf::st_write(output_geojson) |>
  #       suppressWarnings() # GDAL generates spurious warnings on linux
  #   }
  #
  #   # remove the large source raster from memory
  #   rm(r_i)
  #   gc()
  #   t2 = proc.time()
  #   cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  # }
  #
  # # report all files written
  # output_paths = output_nc
  # if(write_csv) output_paths = c(output_paths, output_csv, output_geojson)
  # return(output_paths)
}


