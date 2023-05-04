
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
  df_pad = data.frame(time=t_obs) |> archive_pad(t='time', step_hours, quiet=TRUE)
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
#' This function allows us to trick `nc_update` into not checking for this variable
#' on prior dates which speeds loading.
#'
#' @param grib_df data frame returned by `grib_list(..., dupe=FALSE)`
#' @param aoi geometry object passed to `grib_idx` (area of interest)
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
    r_dummy = grib_idx(grib_df, aoi=aoi, quiet=TRUE, try_again=TRUE)[['r_aoi']] |>
      terra::rast(nlyrs=length(t_NA)) |>
      stats::setNames(t_NA)

    # write dummy file to disk (and attributes JSON)
    r_dummy[] = NA
    terra::time(r_dummy) = t_NA
    nc_write(r_dummy, output_nc, overwrite=TRUE, append=FALSE)
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
    s_index_path = file_wx('index', base_dir, s, s_nm)
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
