#' Find index of NA layers and times available in a (set of) NetCDF files(s)
#'
#' The function returns a list of vectors indicating the times (layers) available
#' in one or more NetCDF files (at `nc_path`), along with information about NAs.
#'
#' If `nc_path` is a single path, or if `join=TRUE`, the function returns in a list:
#'
#' * "na": index of layers containing any number of NAs
#' * "time" : all times as POSIXct strings, in same order as layers
#' * "time_obs" : elements of "time" whose layers are not listed "na"
#' * "time_na" : elements of "time" whose layers are listed "na"
#'
#' If `nc_path` is a vector of multiple paths and `join=FALSE`, the function returns
#' a list of results of with the above structure (one per path).
#'
#' `join=TRUE` allows a time series to be split into different files, each with a
#' (possibly) different chunk of times. The function checks all paths and joins the
#' times that it finds, returning all results in a single list. If a time is observed
#' in any of the paths in `nc_path`, it is marked as observed in the output.
#'
#' If a NetCDF file is listed in `nc_path` but not not found on disk, it is dropped
#' from the results. When looking up times, the function calls `time_json` first, then
#' the much slower `time_nc` only if the JSON is not found. In that case the function
#' calls `write_time_json` to create it.
#'
#' If none of the paths in `nc_path` point to an existing file, the function returns
#' an empty list
#'
#' @param nc_path character vector, path(s) to the NetCDF file(s)
#' @param join logical, when TRUE the results from all `nc_path` are combined
#' @param collapse logical, if TRUE and output list has length 1, it is unlisted
#'
#' @return a list with vectors 'na' (integer), 'time', 'time_na', 'time_obs' (or a list of them)
#' @export
time_wx = function(nc_path, join=TRUE, collapse=TRUE) {

  # if requested, first attempt to get the info from the JSON (fast)
  time_result = time_json(nc_path)
  is_pending = is.na(time_result)

  # create any missing JSON files
  if( any(is_pending) ) {

    # write the JSONs
    is_json_new = nc_path[is_pending] |> sapply(write_time_json)
    if( any(is_json_new) ) {

      # load them
      idx_new = which(is_pending)[is_json_new]
      time_result[idx_new] = nc_path[idx_new] |> time_json()
    }
  }

  # at least one of the paths in `nc_path` should point to a valid file
  is_failed = is.na(time_result)
  if( all(is_failed) ) return( list() )

  # join lists from different chunks of the same time series
  if( join & ( sum(!is_failed) > 1 ) ) {

    # merge all times, handling degenerate case list(NULL, x) by dropping empty entries
    list_all = time_result[!is_failed] |> lapply(\(x) x[['time']])
    time_all = do.call(c, list_all[sapply(list_all, length) > 0]) |>
      unique() |>
      sort()

    # the same for observed times
    list_all_obs = time_result[!is_failed] |> lapply(\(x) x[['time_obs']])
    time_obs = do.call(c, list_all_obs[sapply(list_all_obs, length) > 0]) |>
      unique() |>
      sort()

    # update NA index for consistency
    na_all = which( !(time_all %in% time_obs) )
    return(list(na = na_all,
                time = time_all,
                time_obs = time_obs,
                time_na = time_all[na_all]))
  }

  # collapse length-1 lists
  list_out = time_result[!is_failed]
  if( collapse & ( length(list_out) == 1 ) ) list_out = list_out[[1]]
  return( list_out )
}


#' Read a NetCDF time series file to find times and layers with NAs
#'
#' For SpatRaster `r`, the function returns a list of four vectors:
#'
#' * "na": index of layers containing any number of NAs
#' * "time" : all times as POSIXct strings, in same order as layers
#' * "time_obs" : elements of "time" whose layers are not listed "na"
#' * "time_na" : elements of "time" whose layers are listed "na"
#'
#' or `NA`, if the file is not found. All times are assumed to be in the UTC time zone.
#'
#' If `r` is a vector of paths, the function opens each one in a loop
#' and returns results as a list (same length as `r`) of lists with
#' the above structure. Top-level names correspond to the path
#' that was read.
#'
#' This function can be slow with large files, so we recommend calling
#' once only, then caching results in a JSON (see `?write_time_json`)
#' where they can be retrieved much more quickly using `time_json`.
#'
#' @param r SpatRaster or character vector of paths to NetCDF file(s)
#'
#' @return a list of list(s) with vectors 'na' (integer), 'time', 'time_na', 'time_obs'
#' @export
time_nc = function(r) {

  # return layer info for the raster
  if( is.null(r) ) return( NA )
  if( is(r, 'SpatRaster') ) {

    # load the times then check data values for NAs
    cat('\nchecking', terra::nlyr(r), 'layers for NAs')
    n_obs = r |> terra::global('notNA') |> as.matrix() |> as.integer()
    is_na = n_obs < terra::ncell(r)
    r_time = terra::time(r)
    return(list(na = which(is_na),
                time = r_time,
                time_na = r_time[is_na],
                time_obs = r_time[!is_na]))
  }

  # character input interpreted as path to nc file
  if( is(r, 'character') ) {

    # loop for vectorized case
    r_result = r |> lapply(\(p) {

      if( file.exists(p) ) time_nc(terra::rast(p)) else NA

    }) |> stats::setNames(r)

    return(r_result)
  }

  stop('r had unrecognized class (expected character or SpatRaster)')
}


#' Open JSON file(s) associated with NetCDF time series
#'
#' To improve access and write times for large NetCDF files, we use a JSON file
#' to index times and layers containing NAs. This JSON is located in the "time"
#' sub-directory (relative to the NetCDF file) and has the same name up to the
#' file extension.
#'
#' This function returns the data from the JSON in a list, along with a pair of
#' derived time vectors (see `?time_wx`). If the JSON file is not found, the
#' function instead returns NA (see `?write_time_json`). If a vector of paths is
#' supplied, the function loops over them and returns all results in a list.
#'
#' For consistency with the vectorized case, the function returns a length-1
#' list when a single path is supplied to `nc_path`.
#'
#' @param nc_path character vector, path(s) to the NetCDF file(s)
#'
#' @return a list of list(s) with vectors 'na' (integer), 'time', 'time_na', 'time_obs'
#' @export
time_json = function(nc_path) {

  # expect like-named JSONs in subdirectory "time"
  json_nm = nc_path |> tools::file_path_sans_ext() |> basename() |> paste0('.json')
  json_path = file.path(dirname(nc_path), 'time', json_nm)

  # loop for vectorized case
  seq_along(json_path) |> lapply(\(i) {

    if( file.exists( json_path[i] ) ) {

      # load JSON data as list
      json_data = json_path[i] |> readLines() |> jsonlite::fromJSON()
      times = json_data[['time']]

      # split observed and unobserved
      is_na = seq_along(times) %in% json_data[['na']]
      json_data[['time_na']] = times[is_na]
      json_data[['time_obs']] = times[!is_na]

      # convert character to POSIXct
      nm_posix = c('time', 'time_na', 'time_obs')
      my_char2p = \(x) if( length(x) > 0 ) as.POSIXct(x, tz='UTC') else list()
      json_data[nm_posix] = json_data[nm_posix] |> sapply(my_char2p)
      json_data

    } else { NA }

  }) |> stats::setNames(json_path)
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
#' the layers in `r`. In this case all times in `r` must lie after the existing times
#' in `nc_path`.
#'
#' To delete an existing JSON file and replace it with layer info from SpatRaster `r`,
#' set `append=FALSE`.
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
  cat('\nindexing NAs and writing to', json_path)

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

      # validity check
      final_time = json_data[['time']] |> max()
      msg_final = paste('all times in r must lie after', final_time)
      if( !all( r_index[['time']] > final_time ) ) stop(msg_final)

      # update important fields
      n_existing = length(json_data[['time']])
      json_data[['time']] = json_data[['time']] |> c(r_index[['time']])
      json_data[['na']] = json_data[['na']] |> c(n_existing + r_index[['na']])
    }
  }

  json_data[c('na', 'time')] |> jsonlite::toJSON(pretty=TRUE) |> writeLines(json_path)
  cat(' \U2713')
  return(TRUE)
}
