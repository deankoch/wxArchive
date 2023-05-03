


#' Expands/subsets the tibble returned by `my_archive_lister` to get gap-less time series
#'
#' This creates a time series with interval `step_hours`, starting from the earliest time
#' (`times` column) found in `grib_df`. Rows inconsistent with the series are removed, and
#' new rows (with NA fields apart from `times`) are added to for missing times.
#'
#' @param grib_df data frame containing column `t`
#' @param times character column name containing `POSIXct` times.
#' @param quiet logical indicating to suppress console output
#'
#' @return a tibble with the same columns as `grib_df`
#' @export
my_archive_padder = function(grib_df, times='posix_pred', until=NULL, quiet=FALSE) {

  # times should be a column name of the first argument
  if( !(times %in% names(grib_df)) ) stop(paste(times, 'not found in grib_df'))
  grib_df = grib_df |> dplyr::arrange(get(times))

  # detect step size and identify gap sizes
  include_df = grib_df |> my_detect_step(times=times, more=TRUE)
  step_hours = grib_df |> my_detect_step(times=times, more=FALSE)
  n_omit = nrow(grib_df) - nrow(include_df)

  # identify times appearing more than once (eg forecast files from different hours)
  time_dupe = include_df |>
    dplyr::filter(interval == 0) |>
    dplyr::pull(get(times)) |>
    as.character() |>
    paste(collapse=', ')

  # this omits the second of any duplicates
  include_df[['interval']][1L] = step_hours
  include_df = include_df |> dplyr::filter(interval > 0)
  n_omit = nrow(grib_df) - nrow(include_df) - n_omit
  if(!quiet & (n_omit > 0) ) cat('omitted', n_omit, 'duplicate(s):', time_dupe, '\n')

  # filter date/times out of alignment with first, given interval step_hours
  include_df = include_df |> dplyr::filter( (ts_hours %% step_hours ) == 0 )
  n_omit = nrow(grib_df) - nrow(include_df) - n_omit
  if(!quiet & (n_omit > 0) ) cat('omitted', n_omit, 'file(s) with misaligned date/time\n')

  # set default end time
  to_grib = max(include_df[[times]])
  if( !is.null(until) ) {

    # extend time series as needed
    hours_added = ( as.integer(until) - as.integer(to_grib) ) / (60 * 60)
    if(hours_added > step_hours) {

      if(!quiet) cat('extending by', floor(hours_added/step_hours), 'time steps')
      to_grib = until
    }
  }

  # set up a data frame with regular sequence of times covering the input
  from_grib = min(include_df[[times]])
  missing_df = data.frame(seq(from_grib, to_grib, 60 * 60 * step_hours)) |> setNames(times)

  # join with existing data
  out_df = dplyr::right_join(include_df, missing_df, by=times) |> dplyr::arrange(get(times))

  # new column indicating gap length (or 0 for continuous)
  out_df = out_df |> dplyr::mutate(gap = my_look_ahead(interval-2L))

  # print info before returning the tibble
  n_miss = out_df[['ts_hours']] |> is.na() |> sum()

  # date range and time step
  msg_time = paste(format(from_grib, tz='UTC'), 'to', format(to_grib, tz='UTC'),
                   paste0('UTC (', step_hours, ' hour interval)\n'))

  # size and missingness
  msg_files = paste(nrow(out_df), 'time points', paste0('(', n_miss, ' missing)'))

  if( !quiet ) cat('\ntime series: ', paste0(msg_time, msg_files))
  return( dplyr::tibble(out_df) )
}


# find the index of first match to each element in pattern to the choices in available.
# pattern should be character vector of regular expressions
# if pattern is NULL the function instead returns all indices 1...length(available)
which_lyr = function(available, pattern=NULL, na.rm=FALSE, quiet=FALSE) {

  if(is.null(pattern)) return(stats::setNames(seq(length(available)), available))

  # literal matching with regex=FALSE
  pattern_map = match(pattern, available)

  # grep to get list of matches (expect 1 per nm)
  pattern_map_list = pattern |> lapply(\(x) grep(x, available))

  # deal with multiple matches
  is_replicate = sapply(pattern_map_list, length) > 1
  if( any(is_replicate) ) {

    # warn before discarding matches
    msg_pattern = paste(pattern[is_replicate], collapse=', ')
    msg_available = pattern_map_list[is_replicate] |> sapply(\(x) paste(x[-1], collapse=', '))
    msg_problem = paste('variable(s)', msg_pattern, 'also matched:', msg_available)
    if(!quiet) warning(msg_problem)

    # extract the first match only
    pattern_map = pattern_map_list |> sapply(\(x) x[1])

  } else { pattern_map = do.call(c, pattern_map_list) }

  if(na.rm) pattern_map = pattern_map[!is.na(pattern_map)]
  return(pattern_map)
}



# get or set attributes na, time, time_na, time_obs within large multilayer
# geoTIFF with associated JSON this is helpful for de-cluttering, and for
# speeding up loading files with lots of NA layers (as we don't request them
# in the rast(...)[] call)

# fills missing values (NAs) with the next non-NA value (based on zoo)
my_look_ahead = function(x) {

  is_obs = !is.na(x)
  rev(c(NA, rev(x[is_obs]))[ 1L + cumsum(rev(is_obs)) ])
}

#
#' Return a tibble of start/end indices for all NA sequences in a vector
#'
#' Identifies the start and end points of all gaps (contiguous stretches of `NA` values)
#' in a vector, along with their lengths. If `invert=TRUE` the function returns
#' information about the non-`NA` segments.
#'
#' `x` can be any vector containing `NA` entries. If `times` is supplied, The function
#' returns the value of `times` at the start/end index, instead of the index itself.
#'
#' If a data frame is passed to `x` then `times` should be a (character) column
#' name identifying the time values to use - the first (non-time) column in `x`
#' is then checked for `NA`s.
#'
#' If there are no `NA`s (or, with `invert=TRUE`, no non-NA's), the function returns NULL
#'
#' @param x
#' @param start_only
#' @param times
#' @param invert
#'
#' @return a tibble containing a row for each distinct gap
#' @export
#'
#' @examples
my_gap_finder = function(x, start_only=FALSE, times=NULL, invert=FALSE) {

  if( is.data.frame(x) ) {

    # if x is a data frame then times must give a column name
    is_time = names(x) == times[1]
    if( is.null(times) | !any(is_time) ) stop('expected "times" to be a column name of x')
    if( ncol(x) < 2 ) stop('expected a data frame with 2 or more columns')

    # proceed with first column of data frame (after removing time column)
    times = x[[ which(is_time)[1] ]]
    x = x[[ which(!is_time)[1] ]]
    return( my_gap_finder(x, start_only=start_only, times=times, invert=invert) )
  }

  # find non-NA entries
  is_obs = !is.na(x)
  if( invert ) is_obs = !is_obs
  if( all(is_obs) ) return(data.frame())

  # compute first index of all gaps
  i = which(!is_obs)
  i_start = head(i, 1) |> c( tail(i, -1)[ which( diff(i) != 1 ) ] )
  if( start_only ) return(i_start)

  # compute last index of all gaps
  i_end = if( all(!is_obs) ) length(is_obs) else {

    # create dummy vector with NAs in the right place (invert prevents using x)
    x_dummy = seq_along(is_obs) |> match(which(is_obs))

    # reverse on either end gets us sequence endpoints, but in reversed index
    i_end_inv = x_dummy |> rev() |> my_gap_finder(start_only=TRUE) |> rev()
    i_end = 1L + length(is_obs) - i_end_inv
  }

  # return in a tibble
  i_len = 1L+i_end-i_start
  out_df = data.frame(start=i_start, end=i_end, length=i_len) |> dplyr::tibble()

  # return start/end times too
  if( !is.null(times) ) {

    out_df[['start_time']] = times[i_start]
    out_df[['end_time']] = times[i_end]
  }

  return(out_df)
}


# estimate memory usage for a subset of the forecast archive
# gdim should be the grid dimensions (integer vector length 2)
estimate_memory = function(gdim, n_var=1L, n_file=1L, quiet=FALSE) {

  # memory requirement estimate (based on NCmisc)
  est_gb = 1.05 * prod(c(gdim, n_file, n_var) ) / 2^27
  if(!quiet) cat('\nexpected memory usage:', ifelse(est_gb < 1,
                                                    round(1e3*est_gb, 2) |> paste('MB'),
                                                    round(est_gb, 2) |> paste('GB')))

  return(est_gb)
}



#' Find the step size in a time series as minimal time difference between observations
#'
#' Returns the step size in hours, or if `more=TRUE`, appends the columns "ts_hour" (the
#' number of hours since the first observation) and "interval" (the difference in
#' "ts_hour" between a given observation and the previous). The step size is the least
#' nonzero "interval".
#'
#' @param grib_df
#' @param times
#' @param quiet
#' @param more
#'
#' @return
#' @export
#'
#' @examples
my_detect_step = function(grib_df, times='posix_pred', quiet=FALSE, more=FALSE) {

  # filter NA date/times
  include_df = grib_df |>
    dplyr::filter( !is.na(get(times)) )
  n_omit = nrow(grib_df) - nrow(include_df)
  if( !quiet & (n_omit > 0) ) cat('removed', n_omit, 'rows with NA', times, '\n')
  if( nrow(include_df) < 2 ) stop('grib_df must have at least two non-NA values of ', times)

  # append column for number of hours between a prediction time and the previous
  include_df = include_df |>
    dplyr::mutate(ts_hours = as.numeric(get(times) - min(get(times)), units='hours')) |>
    dplyr::mutate(interval = ts_hours - dplyr::lag(ts_hours))

  # extract time step as the least nonzero time difference
  step_hours = include_df |>
    dplyr::filter(interval > 0) |>
    dplyr::pull(interval) |> min(na.rm=TRUE)

  if( !more ) return(step_hours)
  return( include_df )
}

# sample contiguous blocks of times at random, excluding sets with too many NAs
my_sample_na = function(t_obs, n, p_max=0.5, step_hours=NULL, na_rm=TRUE, i_max=1e5) {

  # extract completed time series and na index
  df_pad = data.frame(time=t_obs) |> my_archive_padder(t='time', step_hours, quiet=TRUE)
  all_times = df_pad[['time']]
  is_na = df_pad[['interval']] |> is.na()
  n_times = length(all_times)

  # loop until a valid sample is found
  i_max=1e3
  i = 0
  while(i < i_max) {

    i = i + 1
    idx_start = seq(n_times - n + 1) |> sample(1)
    idx_sample = idx_start - 1 + seq(n)
    n_na = sum( is_na[idx_sample] )
    if( n_na <= ceiling( n * p_max ) ) i_max = 0
  }

  if(i == i_max) stop('i_max reached. Try decreasing n or increasing p_max')
  t_out = all_times[seq(n) + idx_start - 1L]
  if(na_rm) t_out = t_out[ !is_na[seq(n) + idx_start - 1L] ]
  return(t_out)
}


#' Make an empty NetCDF file covering time period prior to a cutoff date
#'
#' Note that this overwrites any existing data in `output_nc`. The result is a file
#' containing a layer for each time in `grib_df[['posix_pred']]` prior to `date_cutoff`,
#' with all data values set to `NA`.
#'
#' At least one of the files listed in `grib_df` must exist so that this function can
#' read the grid dimensions etc. It is assumed that all files in `grib_df` are of the
#' same spatial grid.
#'
#' `pcp_total` is not available prior to 2016-10-01 (the default for `date_cutoff`).
#' This function allows us to trick `my_update_nc` into not checking for this variable
#' on prior dates which speeds loading.
#'
#' @param grib_df data frame returned by `my_archive_lister(..., dupe=FALSE)`
#' @param aoi geometry object passed to `my_grib_idx` (area of interest)
#' @param output_nc the output file name
#'
#' @return nothing, but overwrites the file `output_nc`
#' @export
my_dummy_nc = function(output_nc, grib_df, aoi, date_cutoff=as.Date('2016-10-01')) {

  # all rows in grib_df prior to date_cutoff
  t_NA = grib_df |> dplyr::filter(date_rel < date_cutoff) |> pull(posix_pred)
  if( length(t_NA) > 0 ) {

    # open first available file in grib_df
    cat('\ncreating empty file', basename(output_nc))
    r_dummy = my_grib_idx(grib_df, aoi=aoi, quiet=TRUE, try_again=TRUE)[['r_aoi']] |>
      terra::rast(nlyrs=length(t_NA)) |>
      stats::setNames(t_NA)

    # write dummy file to disk (and attributes JSON)
    r_dummy[] = NA
    terra::time(r_dummy) = t_NA
    my_nc_write(r_dummy, output_nc, overwrite=TRUE, append=FALSE)
  }
}

# Return a data frame of information about time coverage of different time series
my_archive_dates = function(base_dir, var_nm=NULL, sub_dir=NULL) {

  # by default check all subdirectories except the GRIB source files
  if( is.null(sub_dir) ) {

    base_files = base_dir |> list.files()
    base_info = base_dir |> file.path(base_files) |> file.info()
    sub_dir = base_files[ ( base_files != 'grib' ) & ( base_info[['isdir']] ) ]
  }

  # loop over sub directories
  nc_info_list = sub_dir |> stats::setNames(sub_dir) |> lapply(\(s) {

    # list all files, filter to nc files
    s_files = base_dir |> file.path(s) |> list.files()
    nc_path = base_dir |> file.path(s, s_files[ endsWith(s_files, '.nc') ])
    s_nm = nc_path |> tools::file_path_sans_ext() |> basename()
    names(nc_path) = s_nm

    # filter to relevant nm
    if( !is.null(var_nm) ) s_nm = s_nm[ s_nm %in% var_nm ]
    if( length(s_nm) == 0 ) return( NULL )

    # read all attributes JSON files to find dates info
    s_index_path = wx_file('index', base_dir, s, s_nm)
    var_df_list = s_nm |> stats::setNames(s_nm) |> lapply(\(nm) {

      var_attr = my_nc_attributes(nc_path[nm])
      is_indexed = !is.null(var_attr)
      df_out = data.frame(name=nm, indexed=is_indexed)
      t_all = if(is_indexed) var_attr[['time']] else as.POSIXct(NA)
      t_obs = if(is_indexed) var_attr[['time_obs']] else as.POSIXct(NA)
      data.frame(name = nm,
                 indexed = is_indexed,
                 n = sum(!is.na(t_all)),
                 n_obs = sum(!is.na(t_obs)),
                 n_na = sum(!is.na(t_all)) - sum(!is.na(t_obs)),
                 start = min(t_all),
                 end = max(t_all),
                 start_obs = min(t_obs),
                 end_obs = max(t_obs))


    })

    # collapse all variables into a single data frame
    data.frame(sub=s) |> cbind( do.call(rbind, var_df_list) )

  })

  # collapse all sub-directories into a single data frame
  do.call(rbind, nc_info_list) |> dplyr::tibble()
}
