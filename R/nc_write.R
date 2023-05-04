#' Create or modify a NetCDF file by adding new time layers from SpatRaster
#'
#' This creates or updates an existing NetCDF file by adding the layers (times)
#' in the supplied SpatRaster `r`.
#'
#' The function does little checking to see if the request is valid. It assumes
#' that `p` is a valid path to a NetCDF file (.nc), and that raster grid in `r`
#' is compatible with the file at `p` in terms of dimensions etc.
#'
#' `r` must have a time for each layer (check `terra::time(r)`). Only new times
#' are copied - ie those not found in the file at `p`.
#'
#' If `overwrite=FALSE` (the default), the function overwrites nothing but
#' returns the times that would be added to the file.
#'
#' When overwriting an existing file, the function initially writes to a temporary
#' file in the same directory. At the end of the function call the old file is
#' deleted and the temporary one is renamed to take its place. If something goes
#' wrong and the function halts, this temporary file can be safely deleted.
#'
#' After writing to the .nc file, the function creates/updates the JSON attributes
#' file in the 'time' subdirectory.
#'
#' @param r SpatRaster with POSIXct vector `terra::time(r)`, the data to write
#' @param p character path the (.nc) time series data file to write
#' @param overwrite logical indicating to write changes to the nc file on disk
#'
#' @return vector of POSIXct times, the layers added to the file
#' @export
nc_write = function(r, p, overwrite=FALSE) {

  # create/load JSON for nc file at p and copy times
  is_update = file.exists(p)
  p_time = if( is_update ) time_wx(p)[['time']] else NULL

  # collect times from input raster
  r_time = terra::time(r)
  if( is.null(r_time) ) stop('terra::time(r) must return a time for each layer in r')

  # append mode omits duplicates in r, otherwise omit duplicates in nc file
  is_new = !( r_time %in% p_time )
  p_time_fetch = as.POSIXct(p_time[ !(p_time %in% r_time[is_new]) ])

  # if there's nothing to add we are done
  if( !any(is_new) ) {

    # compute index JSON if it's not there
    cat('\nup to date')
    return( as.POSIXct(integer(0)) )
  }

  # +0 forces terra to load data into RAM
  cat('\nloading and sorting', sum(is_new), 'input SpatRaster layer(s)')
  r_add = r[[ which(is_new) ]] + 0

  # merge new and old layers in memory
  if( !is_update | ( length(p_time_fetch) == 0 ) ) { r_out = r_add } else {

    cat('\nmerging with', length(p_time_fetch), 'existing nc layer(s)')
    r_existing = p |> nc_layers(times=p_time_fetch, preload=TRUE)
    r_out = c(r_existing, r_add)
  }

  # sort and name output layers
  r_out = r_out[[ order(terra::time(r_out)) ]]
  names(r_out) = paste0('lyr_', seq(terra::nlyr(r_out)))

  # update the nc and JSON files
  if( overwrite ) {

    # a name for the dataset pulled from the file path
    nm = basename(p) |> tools::file_path_sans_ext()
    cat('\nwriting to', p)

    # safer overwrite via tempfile
    if( is_update ) {

      # a temporary file name for the existing file
      suffix_temp = paste0('_', basename(tempfile()), '.nc')
      p_dest = file.path(dirname(p), paste0(nm, suffix_temp))

    } else { p_dest = p }

    # write result
    terra::writeCDF(r_out, p_dest, varname=nm)
    if( is_update ) {

      # remove old file and rename new one to replace it
      unlink(p)
      file.rename(from=p_dest, to=p)
      p_dest = p
    }

    # if all new times appear at the end use the (faster) append mode
    cat(' \U2713')
    is_appended = all( p_time < min(terra::time(r_add)) )

    # update attributes
    if( is_appended ) r_add = r_out
    p_dest |> write_time_json(r=r_add, append=is_appended)
  }

  return( r_time[is_new] )
}
