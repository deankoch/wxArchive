#' Open a subset of times from one or more NetCDF files
#'
#' Returns a SpatRaster with time-indexed layers in chronological order,
#' loaded from one or several NetCDF files. This wraps `terra::rast`, setting
#' the `lyrs` argument by matching `times` to times listed in the JSON file(s)
#' in subdirectory "/time".  If `times` is `NULL`, the function returns all
#' times found in the first file `p[1]`.
#'
#'`p` can be a vector of paths, indicating to match `times` to the times in each
#' file. If a time is found in multiple files, the function returns the matching
#' layer from the first file (WRT to the order in `p`), and ignores the rest.
#' The returned SpatRaster will include one layer for all `times` that could be
#' matched, with different layers possibly coming from different files.
#'
#' This functionality allows you to split time series into temporal chunks,
#' and store them in multiple files. Note that the function assumes all files
#' in `p` contain data for the same variable (variable names are ignored).
#'
#' The function creates the JSON file(s) as needed, and assumes that
#' all of `times` can be found in the .nc file(s) at `p`. Available times
#' for the file at `p` can be listed with `terra::rast(p) |> terra::time()`,
#' and this should match the "time" entries of `time_wx(p)`, assuming the JSONs
#' are up to date.
#'
#' When `preload=TRUE` the function forces `terra` to copy all data values
#' into memory before returning the SpatRaster handle. Normally `terra` loads
#' values only when they are needed, but I have found that this can create
#' problems with file operations when there are a large number of times
#' (layers) in the file.
#'
#' When `na_rm=TRUE`, only layers with all grid points observed (non-NA) are
#' searched. In this case times that match a layers labelled as NA are reported
#' as "unmatched" (and no NA layers are returned).
#'
#' @param p character vector path to the nc file(s)
#' @param times vector of unique POSIXct times to match in the file(s)
#' @param preload logical indicating to load values into RAM
#' @param na_rm logical indicating to return only non-NA layers
#'
#' @return a SpatRaster with `length(t)` layers
#' @export
nc_layers = function(p, times=NULL, preload=TRUE, na_rm=FALSE) {

  # filter nonexistent files
  p_valid = file.exists(p)
  if( !any(p_valid) ) stop('file(s) not found: ', paste(p, collapse=', '))
  p = p[p_valid]

  # loop over vectorized input
  load_all_times = is.null(times)
  times = unique(times)
  r_out_list = seq_along(p) |> lapply(\(x) NULL)
  for(i in seq_along(p)) {

    # load/create attributes JSON and find matching times in this file
    attr_i = time_wx(p[i])
    time_available_i = if( na_rm ) attr_i[['time_obs']] else attr_i[['time']]
    if( load_all_times ) times = time_available_i
    is_i = times %in% time_available_i
    if( any(is_i) ) {

      # return the matching layers (+0 forces copying into memory)
      r_out_list[[i]] = terra::rast(p[i], lyrs=match(times[is_i], attr_i[['time']]))
      if( preload ) r_out_list[[i]] = r_out_list[[i]] + 0

      # remove from to-do stack
      times = times[!is_i]
    }
  }

  # report any unmatched times
  msg_unmatched = as.character(times, tz='UTC') |> paste(collapse=', ')
  if( length(times) > 0 ) cat('\nunmatched times:', msg_unmatched)

  # remove unmatched files outputs and sort by time
  r_out = r_out_list[!sapply(r_out_list, is.null)] |> terra::rast()
  r_out = r_out[[ order(terra::time(r_out)) ]]
  return(r_out)
}
