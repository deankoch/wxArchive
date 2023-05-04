#' DEPRACATED DO NOT USE
#' Get or set fields in a JSON file on disk representing attributes of a time series
#'
#' This uses `jsonlite` to write information about the observed times in the netCDF
#' files (`overwrite=TRUE`), or else reads it back to produce the following in a list:
#'
#' * "na": index of layers containing any number of NAs
#' * "time" : all times as POSIXct strings (in same order as layers)
#' * "time_obs" : elements of "time" whose layers are not listed "na"
#' * "time_na" : elements of "time" whose layers are listed "na"
#'
#' All times are assumed to be in the UTC time zone.
#'
#' In read mode (`overwrite=FALSE`), the function returns all fields in the JSON as
#' a list, along with "time_obs" and "time_na", which are computed from "na" and "time".
#' This happens without checking the contents of the nc file. If the file
#' doesn't exist in read mode, the function returns NULL.
#'
#' In write mode (`overwrite=TRUE`), only the fields "na" and "time" are written to
#' the JSON. This should leave any other existing fields alone unless there is a name
#' collision, or the structure is too complicated for a return trip through `fromJSON`
#' and `toJSON`. If `lazy=TRUE` the function will only write the JSON if it doesn't
#' already exist, and otherwise does not modify any existing ones.
#'
#' If a SpatRaster is passed to in `r`, the function will append the NA index and times
#' in `r` to the existing lists for the raster at `nc_path`. This is useful when
#' appending a small number of new times to a large existing time series, where it
#' would be slow and redundant to check for NAs in existing layers.
#'
#' Pass a vector to `nc_path` in read mode to get results for each file back in a list.
#' If `ch=TRUE`, the function treats these files as chunks of a single time series, and
#' the list output is collapsed to into a single list for all of the files - 'time' lists
#' times appearing in any of the files; 'time_obs' lists the times at which a non-NA layer
#' is found in at least one of the files; 'time_na' lists the times at which every
#' appearance in every file is an NA layer; and 'na' gives the index of `time_na` in `time`.
#'
#' By default the JSON is located in sub-directory 'time' and has the same name
#' as the nc (except for file extension). Change this with `json_path`.
#'
#' @param nc_path character path to the (.nc) time series data file
#' @param json_path character path to the JSON (.json) attributes file
#' @param overwrite logical indicating to compute fields and write to disk
#' @param lazy logical indicating to not modify any existing JSONs
#' @param ch logical indicates to collapse results for chunked files into a single list
#' @param r SpatRaster to check for appended times (avoids checking `rast(nc_path)`)
#'
#' @return a list, the contents of the JSON (after updating, if `overwrite=TRUE`)
#' @export
my_nc_attributes = function(nc_path, json_path=NULL, overwrite=FALSE, lazy=FALSE, ch=FALSE, r=NULL) {

  # by default we should have like-named JSONs in subdirectory "time"
  json_expected = nc_path |> tools::file_path_sans_ext() |> basename() |> paste0('.json')
  if( is.null(json_path) ) json_path = file.path(dirname(nc_path), 'time', json_expected)

  # length check
  if( length(json_path) != length(nc_path) ) stop('argument lengths did not match')
  is_added = !is.null(r) & file.exists(json_path)

  # loop for vectorized case
  list_out = seq_along(json_path) |> lapply(\(i) {

    # read mode
    if( !overwrite ) {

      if( !file.exists(json_path[i]) ) { NULL } else {

        # load associated JSON as list
        out_list = json_path[i] |> readLines() |> jsonlite::fromJSON()

        # create `POSIXct` times
        nm_posix = c('time', 'time_na', 'time_obs')
        is_na = seq_along(out_list[['time']]) %in% out_list[['na']]
        out_list[['time_na']] = out_list[['time']][is_na]
        out_list[['time_obs']] = out_list[['time']][!is_na]
        out_list[nm_posix] = out_list[nm_posix] |>
          sapply(\(t) { if( length(t) > 0 ) as.POSIXct(t, tz='UTC') else list() } )

        out_list
      }

    } else {

      # write mode: make the directory if necessary
      p = json_path[i]
      if( !file.exists(dirname(p)) ) dir.create(dirname(p))
      p_exists = file.exists(p)

      # create the JSON fields unless we're lazy and the file exists already
      if( !p_exists | !lazy | is_added[i] ) {

        # load existing JSON attributes as list (if any)
        out_list = list()
        if( p_exists ) out_list = my_nc_attributes(nc_path[i], overwrite=FALSE, lazy=FALSE)

        # if r not supplied as argument, open nc file and compute its NA index
        cat('\nindexing NA layers and writing to', p)
        if( !is_added[i] ) r = nc_path[i] |> terra::rast()
        time_r = terra::time(r)
        n_existing = length(out_list[['time']])
        n_na = r |> terra::global('isNA') |> as.matrix() |> as.numeric()
        idx_r_na = which( n_na > 0 )

        # build output NA index and times
        if( is_added[i] ) {

          # sanity check for append mode
          is_dupe = out_list[['time']] %in% time_r
          msg_dupe = out_list[['time']][is_dupe] |> paste(collapse=', ')
          if( any(is_dupe) ) warning( paste('time(s) found in both r and the nc file:', msg_dupe) )

          # times for every layer and index of NAs
          out_list[['na']] = c(out_list[['na']], idx_r_na + n_existing)
          out_list[['time']] = c(out_list[['time']], time_r)

        } else {

          # if r was not supplied or JSON not found these are already complete
          out_list[['na']] = idx_r_na
          out_list[['time']] = time_r
        }

        # omit names added by this function
        out_list = out_list[ !( names(out_list) %in% c('time_na', 'time_obs') ) ]

        # export to JSON and write changes to disk
        out_list |> jsonlite::toJSON(pretty=TRUE) |> writeLines(p)
        cat(' \U2713')
      }

      # return read mode output for the existing JSON
      my_nc_attributes(nc_path[i], overwrite=FALSE, lazy=FALSE)
    }

  }) |> stats::setNames( names(nc_path) )

  # remove NULL results from missing files
  is_missing = sapply(list_out, is.null)
  if( all(is_missing) ) { return(NULL) } else { list_out = list_out[!is_missing] }

  # collapse lists representing a time series in chunks
  if( ch & length(list_out) > 1 ) {

    # merge and sort all indices
    time_all_list = lapply(list_out, \(x) x[['time']])
    time_obs_list = lapply(list_out, \(x) x[['time_obs']])
    time_all = do.call(c, time_all_list[sapply(time_all_list, length) > 0]) |> unique() |> sort()
    time_obs = do.call(c, time_obs_list[sapply(time_obs_list, length) > 0]) |> unique() |> sort()
    na_all = which(!(time_all %in% time_obs))

    # simplified list
    list_out = list(na = na_all,
                    time = time_all,
                    time_obs = time_obs,
                    time_na = time_all[na_all])
  }

  # collapse length-1 lists
  if( length(list_out) == 1 ) list_out = list_out[[1]]
  return(list_out)
}



