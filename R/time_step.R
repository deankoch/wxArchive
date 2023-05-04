#' Find the minimal time difference between subsequent observations in a time series
#'
#' Returns the step size in hours or, if `more=TRUE`, appends the columns "ts_hour" (the
#' number of hours since the first observation) and "interval" (the difference in
#' "ts_hour" between a given observation and the previous). The step size is the least
#' nonzero "interval".
#'
#' NA times are removed from the results. The default `quiet=FALSE` warns if this happens.
#
#' @param grib_df data frame with POSIXct times in column 'times'
#' @param times character column name of times
#' @param quiet logical, if TRUE the function warns about NAs
#' @param more logical, if TRUE the function returns a data frame
#'
#' @return numeric (hours) or copy of `grib_df` with NA times removed and two columns appended
#' @export
time_step = function(grib_df, times='posix_pred', quiet=FALSE, more=FALSE) {

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

