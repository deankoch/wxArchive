#' Expands/subsets the tibble returned by `grib_list` to get gap-less time series
#'
#' This creates a time series with interval `step_hours`, starting from the earliest time
#' (`times` column) found in `grib_df`. Rows inconsistent with the series are removed, and
#' new rows (with NA fields apart from `times`) are added to for missing times.
#'
#' @param grib_df data frame containing column `t`
#' @param times character column name containing `POSIXct` times.
#' @param quiet logical indicating to suppress console output
#' @param until POSIXct time, the series is truncated/padded to this time point
#'
#' @return a tibble with the same columns as `grib_df`
#' @export
archive_pad = function(grib_df, times='posix_pred', until=NULL, quiet=FALSE) {

  # times should be a column name of the first argument
  if( !(times %in% names(grib_df)) ) stop(paste(times, 'not found in grib_df'))
  grib_df = grib_df |> dplyr::arrange(get(times))

  # detect step size and identify gap sizes
  include_df = grib_df |> time_step(times=times, more=TRUE)
  step_hours = grib_df |> time_step(times=times, more=FALSE)
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

  # Fill missing values (NAs) with the next non-NA value
  look_ahead = \(x) {

    # based on code in the `zoo` package
    is_obs = !is.na(x)
    rev(c(NA, rev(x[is_obs]))[ 1L + cumsum(rev(is_obs)) ])
  }

  # new column indicating gap length (or 0 for continuous)
  out_df = out_df |> dplyr::mutate(gap = look_ahead(interval-2L))

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
