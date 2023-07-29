#' Find index of NA layers and times available in a (set of) NetCDF or GeoTIFF files(s)
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
#' Some or all of the elements in `nc_path` may point to directories containing nc
#' files, in which case all of the contents are scanned (see `?nc_chunk`).
#'
#' If a NetCDF file is listed in `nc_path` but not not found on disk, it is dropped
#' from the results. When looking up times, the function calls `time_json` first, then
#' the much slower `time_nc` only if the JSON is not found. In that case the function
#' calls `write_time_json` to create it.
#'
#' If none of the paths in `nc_path` point to an existing file, the function returns
#' an empty list
#'
#' If `file_ext='tif'`, all of the above applies, but the function looks for an output
#' file of type GeoTIFF (rather than NetCDF), having extension '.tif'.
#'
#' @param nc_path character vector, path(s) to the NetCDF file(s)
#' @param join logical, when TRUE the results from all `nc_path` are combined
#' @param file_ext character, either 'tif' or 'nc'
#'
#' @return a list with vectors 'na' (integer), 'time', 'time_na', 'time_obs' (or a list of them)
#' @export
time_wx = function(nc_path, join=TRUE, file_ext='nc') {

  # expand paths to any files chunked by year. This checks existence
  nc_path = do.call(c, lapply(nc_path, \(p) nc_chunk(p, file_ext=file_ext)))
  if( anyNA(nc_path) ) nc_path = nc_path[!is.na(nc_path)]
  if( length(nc_path) == 0 ) return( list() )

  # first attempt to get the info from the JSON (fast)
  time_result = time_json(nc_path)
  is_pending = is.na(time_result)

  # create any missing JSON files
  if( any(is_pending) ) {

    # skip nc files not found on disk
    is_created = is_pending & file.exists(nc_path)
    if( any(is_created) ) {

      # write the JSON for existing nc file(s)
      is_json_new = nc_path[is_created] |> sapply(write_time_json)
      if( any(is_json_new) ) {

        # load them
        idx_new = which(is_pending)[is_json_new]
        time_result[idx_new] = nc_path[idx_new] |> time_json()
      }
    }
  }

  # at least one of the paths in `nc_path` should point to a valid file
  is_failed = is.na(time_result)
  if( all(is_failed) ) return( list() )

  # join lists from different chunks of the same time series
  if( join & ( sum(!is_failed) > 0 ) ) {

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

  return( time_result[!is_failed] )
}


#' Read a NetCDF or GeoTIFF time series file to find times and layers with NAs
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
#' @param r SpatRaster or character vector of paths to NetCDF or GeoTIFF file(s)
#'
#' @return a list with vectors 'na' (integer), 'time', 'time_na', 'time_obs'
#' @export
time_nc = function(r) {

  # return layer info for the raster
  if( is.null(r) ) return( NA )
  if( is(r, 'SpatRaster') ) {

    # load the times then check data values for NAs
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

    #
    if( length(r_result) == 1 ) r_result = r_result[[1]]

    return(r_result)
  }

  # Date vector
  if( is(r, 'Date') ) {

    r_time = seq.Date(min(r), max(r), by='day')
    r_obs = sort(r)
    is_obs = r %in% r_time
    return(list(na = which(!is_obs),
                time = r_time,
                time_na = r_time[!is_obs],
                time_obs = r_time[is_obs]))
  }

  stop('r had unrecognized class (expected Date, character, or SpatRaster)')
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

  # return NA for invalid input
  if( ( length(nc_path) == 0 ) | !is.character(nc_path) ) return(NA)

  # helper function to handle input that could be POSIX date or empty
  my_char2p = \(x) if( length(x) > 0 ) as.POSIXct(x, tz='UTC') else list()

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
      json_data[nm_posix] = json_data[nm_posix] |> sapply(my_char2p)
      json_data

    } else { NA }

  }) |> stats::setNames(json_path)
}
