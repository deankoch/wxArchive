#' Create or modify a NetCDF file by appending new time layers from SpatRaster
#'
#' This creates or updates an existing NetCDF file with new times (layers) in the
#' supplied SpatRaster `r`. With `in_mem=TRUE` the file is not modified and the
#' function instead returns the SpatRaster that would have been written.
#'
#' `r` must have a time for each layer (check `terra::time(r)`). When `insert=FALSE`
#' (the default), only new times in `r` are written, ie those not found already in
#' the file at `p`. Set `insert=TRUE` to allow existing times to be modified in `p`.
#' Existing times that do not appear in `r` are never modified. This means that once
#' a time has been added to a file, it cannot be removed except by deleting the file
#' and starting over.
#'
#' The function assumes that `p` is a valid path to a NetCDF file with extension
#' .nc (and not a directory of them, as in calls to `nc_write_chunk`). It also assumes
#' that the raster grid in `r` is compatible with the file at `p` in terms of
#' dimensions, projection, etc.
#'
#' When writing to an existing file, the function initially writes to a temporary
#' file in the same directory. At the end of the function call the old file is
#' deleted and the temporary one is renamed to take its place. If something goes
#' wrong and the function halts, this temporary file can be safely deleted.
#'
#' After writing to the .nc file, the function creates/updates the JSON attributes
#' file in the 'time' subdirectory. If something goes wrong with this file, it can
#' be safely deleted then rebuilt with `time_wx(p)`.
#'
#' I have encountered some unexpected memory management issues when attempting to
#' write NetCDF files containing grids with many spatial points. In my case, I was
#' calling `terra::writeCDF` on a grid with 365 layers, and total size 3GB in RAM
#' (`terra::mem_info`). I would expect this to require about 6-9 GB of free RAM
#' (`ncdf4` is making its own copy in RAM). However this function call consistently
#' led to peak memory usage exceeding 28GB, and subsequent slowdowns as `terra` is
#' forced to move grids from RAM to disk.
#'
#' I believe this is a garbage collection malfunction, with intermediate (temporary)
#' grid objects not getting removed from RAM, but I don't know enough about `ncdf4`
#' or NetCDF drivers to understand where the issue is happening. The new argument
#' `in_mem=TRUE` should allow users to write other file formats (like GeoTIFF)
#' manually as a workaround. `file_ext` may be changed to `.tif` in this case (do not
#' change it when `in_mem=FALSE`). This will cause the function to scan for existing
#' GeoTIFF output data instead of NetCDF.
#'
#' @param r SpatRaster with POSIXct vector `terra::time(r)`, the data to write
#' @param p character path the (.nc) time series data file/directory to write
#' @param quiet logical suppresses console messages
#' @param in_mem logical if `TRUE` returns a SpatRaster instead of modifying the NetCDF
#' @param file_ext character, either 'tif' or 'nc'
#'
#' @return vector of POSIXct times, the layers added to the file
#' @export
nc_write = function(r, p, quiet=FALSE, insert=FALSE, in_mem=FALSE, file_ext='.nc') {

  # create/load JSON for nc file at p and copy times
  is_update = file.exists(p)
  p_time = if( is_update ) time_wx(p, file_ext=file_ext)[['time']] else NULL

  # collect times from input raster
  r_time = terra::time(r)
  if( is.null(r_time) ) stop('terra::time(r) must return a time for each layer in r')

  # insert mode ignores existing times in `p` that also appear in `r`
  if( insert )  p_time = p_time[ !( p_time %in% r_time ) ]

  # resolve duplicates (in non-insert mode) by omitting layers from input raster
  is_new = !( r_time %in% p_time )
  p_time_fetch = as.POSIXct(p_time[ !(p_time %in% r_time[is_new]) ])

  # if there's nothing to add we are done
  if( !any(is_new) ) {

    if( !quiet ) cat('\nup to date\n')

    # return an empty object
    out_val = NULL
    if(!in_mem) as.POSIXct(integer(0), tz='UTC')
    return(out_val)
  }

  # explicitly remove source raster from RAM after copying relevant layers
  if( !quiet ) cat('\nloading and sorting', sum(is_new), 'input SpatRaster layer(s)')
  r_out = r[[ which(is_new) ]]
  rm(r)
  gc()

  # merge new and old layers in memory
  if( is_update & ( length(p_time_fetch) != 0 ) ) {

    if( !quiet ) cat('\nmerging with', length(p_time_fetch), 'existing nc layer(s)')
    r_out = nc_layers(p, times=p_time_fetch, preload=TRUE, file_ext=file_ext) |> c(r_out)
    gc()
  }

  # sort and name output layers
  r_out = r_out[[ order(terra::time(r_out)) ]]
  names(r_out) = paste0('lyr_', seq(terra::nlyr(r_out)))
  gc()

  # return from memory mode
  if(in_mem) return(r_out)

  # a name for the dataset pulled from the file path
  nm = basename(p) |> tools::file_path_sans_ext()
  if( !quiet ) cat('\nwriting to', p)

  # prepare for safer overwrite via temporary file/directory
  if( is_update ) {

    # a temporary file name for the existing file
    suffix_temp = paste0('_', basename(tempfile()), '.nc')
    p_dest = file.path(dirname(p), paste0(nm, suffix_temp))

  } else { p_dest = p }

  # write result to file
  terra::writeCDF(r_out, p_dest, varname=nm)

  # rename the tempfile if needed
  if( is_update ) {

    # remove old file/directory and rename new one to replace it
    unlink(p, recursive=TRUE)
    file.rename(from=p_dest, to=p)
    p_dest = p
  }

  # update attributes on disk
  p_dest |> write_time_json(r=r_out)
  if( !quiet ) cat('\n')

  # remove remaining SpatRaster object from memory
  rm(r_out)
  gc()

  # return the times added
  return( r_time[is_new] )
}


#' Create or update the JSON file associated with a NetCDF time series
#'
#' This manages the JSON file read by `time_json`. To get the path to this
#' file, do `nc_path |> time_json() |> names()`. This function will create both
#' the file and its parent directory (if they don't exist already). It only modifies
#' the JSON and never the NetCDF file.
#'
#' Note that the JSON contains only the fields "na" and "time" and not the derived
#' fields "time_na" and "time_obs" returned by `time_nc` and `time_wx`.
#'
#' If `r` is supplied, the function indexes the SpatRaster. If `r` is `NULL`, the
#' function instead indexes all layers in the file at `nc_path`. This can be very slow
#' if there are many time points.
#'
#' If `nc_path` doesn't exist, the function returns `FALSE` and writes nothing
#'
#' @param nc_path character path to a single NetCDF file
#' @param r SpatRaster or NULL, containing times to add
#'
#' @return logical indicating if the file was modified
#' @export
write_time_json = function(nc_path, r=NULL) {

  # nc_path must point to a single, existing file
  if( length(nc_path) > 1 ) stop('nc_path had length > 1')
  if( !file.exists(nc_path) ) return(FALSE)

  # output is a like-named JSON in subdirectory "time"
  json_dir = file.path(dirname(nc_path), 'time')
  if( !dir.exists(json_dir) ) dir.create(json_dir)
  json_nm = nc_path |> tools::file_path_sans_ext() |> basename() |> paste0('.json')
  json_path = json_dir |> file.path(json_nm)
  cat('\nwriting time index to', json_path)

  # read times from file if r not supplied
  if( is.null(r) ) r = nc_path
  json_data = time_nc(r)

  # report problems with `time_nc` or write results to disk
  if( anyNA(json_data) ) stop('there was a problem reading times from the SpatRaster')
  json_data[c('na', 'time')] |> jsonlite::toJSON(pretty=TRUE) |> writeLines(json_path)
  return(TRUE)
}
