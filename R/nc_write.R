#' Create or modify a NetCDF file by appending new time layers from SpatRaster
#'
#' This creates or updates an existing NetCDF file with new times (layers) in the
#' supplied SpatRaster `r`.
#'
#' `r` must have a time for each layer (check `terra::time(r)`). Only new times
#' are copied - ie those occurring after the latest time in the file at `p`. To
#' add/modify earlier times you will need to delete the file and start over.
#'
#' The function assumes that `p` is a valid path to a NetCDF file (.nc), and that
#' the raster grid in `r` is compatible with the file at `p` in terms of dimensions,
#' projection, etc.
#'
#' When overwriting an existing file, the function initially writes to a temporary
#' file in the same directory. At the end of the function call the old file is
#' deleted and the temporary one is renamed to take its place. If something goes
#' wrong and the function halts, this temporary file can be safely deleted.
#'
#' After writing to the .nc file, the function creates/updates the JSON attributes
#' file in the 'time' subdirectory. If something goes wrong with this file, it can
#' be safely deleted then rebuilt with `time_wx(p)`.
#'
#' @param r SpatRaster with POSIXct vector `terra::time(r)`, the data to write
#' @param p character path the (.nc) time series data file to write
#'
#' @return vector of POSIXct times, the layers added to the file
#' @export
nc_write = function(r, p) {

  # create/load JSON for nc file at p and copy times
  is_update = file.exists(p)
  p_time = if( is_update ) time_wx(p)[['time']] else NULL

  # collect times from input raster
  r_time = terra::time(r)
  if( is.null(r_time) ) stop('terra::time(r) must return a time for each layer in r')

  # resolve duplicates by omitting layers from input raster
  is_new = !( r_time %in% p_time )
  p_time_fetch = as.POSIXct(p_time[ !(p_time %in% r_time[is_new]) ])

  # if there's nothing to add we are done
  if( !any(is_new) ) {

    # compute index JSON if it's not there
    cat('\nup to date\n')
    return( as.POSIXct(integer(0), tz='UTC') )
  }

  # +0 forces data into RAM
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

  # a name for the dataset pulled from the file path
  nm = basename(p) |> tools::file_path_sans_ext()
  cat('\nwriting to', p)

  # prepare for safer overwrite via tempfile
  if( is_update ) {

    # a temporary file name for the existing file
    suffix_temp = paste0('_', basename(tempfile()), '.nc')
    p_dest = file.path(dirname(p), paste0(nm, suffix_temp))

  } else { p_dest = p }

  # sanity check for attempted append with stale times
  time_start = r_add |> terra::time() |> min()
  is_appended = all( p_time < time_start )
  msg_invalid = paste('cannot append times earlier than the latest existing:', time_start)
  if( is_update & !is_appended ) stop(msg_invalid)

  # write result
  terra::writeCDF(r_out, p_dest, varname=nm)
  if( is_update ) {

    # remove old file and rename new one to replace it
    unlink(p)
    file.rename(from=p_dest, to=p)
    p_dest = p
  }

  # update attributes on disk
  if( is_appended ) r_add = r_out
  p_dest |> write_time_json(r=r_add, append=is_appended)
  cat('\n')

  # remove unused SpatRaster objects from memory
  rm(r_add, r_existing, r_out)
  gc()

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
#' fields "time_na" and "time_obs" (see `?time_wx`).
#'
#' Call the function with default arguments to initialize the JSON for an existing
#' NetCDF file. This involves loading all data into memory, so it can be very slow when
#' there are numerous time points.
#'
#' If `r` is supplied along with `append=TRUE`, the function creates/updates the JSON
#' to make it consistent with the result of concatenating the layers in `nc_path` with
#' the layers in `r`. Only times occurring after the last time in `nc_path` are copied
#' from `r`
#'
#' Set `append=FALSE` to delete an existing JSON file and replace it with all of the
#' layer info for SpatRaster `r`,
#'
#' If `nc_path` doesn't exist, the function returns `FALSE` and writes nothing
#'
#' @param nc_path character path the NetCDF file
#' @param r SpatRaster or NULL, containing times to add
#' @param append logical indicating to append times rather than overwrite
#'
#' @return logical indicating if the file was modified
#' @export
write_time_json = function(nc_path, r=NULL, append=TRUE) {

  # nc_path must point to a single, existing file
  if( length(nc_path) > 1 ) stop('nc_path had length > 1')
  if( !file.exists(nc_path) ) return(FALSE)

  # fetch any existing JSON data in a list, and the file path
  json_data = time_json(nc_path)
  json_path = names(json_data)
  cat('\nwriting attributes to', json_path)

  # collapse the length-1 list
  json_data = json_data[[1]]

  # make the sub-directory if necessary
  json_dir = dirname(json_path)
  if( !dir.exists(json_dir) ) dir.create(json_dir)

  # check for times in input SpatRaster
  r_index = time_nc(r)
  r_exists = !anyNA(r_index)

  # overwrite request erases existing JSON data
  if( !append ) {

    if( !r_exists ) stop('append=FALSE but no data supplied')
    json_data = r_index

  } else {

    # if there is no existing JSON, compute NA index from scratch (slow)
    if( anyNA(json_data) ) json_data = time_nc(nc_path)[[1]]

    # append request merges new times with existing JSON data
    if( r_exists ) {

      # omit too-early times from r time (and NA) index
      final_time = json_data[['time']] |> max()
      is_early = r_index[['time']] <= final_time
      if( any(is_early) ) {

        r_index[['time']] = r_index[['time']][!is_early]
        r_index[['na']] = ( seq_along(is_early) %in% r_index[['na']] )[!is_early] |> which()
      }

      # append to existing times
      n_existing = length(json_data[['time']])
      json_data[['time']] = json_data[['time']] |> c(r_index[['time']])
      json_data[['na']] = json_data[['na']] |> c(n_existing + r_index[['na']])
    }
  }

  json_data[c('na', 'time')] |> jsonlite::toJSON(pretty=TRUE) |> writeLines(json_path)
  return(TRUE)
}
