#' Return a data matrix of seasonal linear predictors for a time series
#'
#' This creates a data matrix for daily and yearly cycles in a regular time series
#' observed at times `t_obs`. With `na_rm=TRUE`, the output has a row for each
#' `t_obs`. Otherwise it has a row for each time in the sequence of frequency
#' `freq_year` spanning `t_obs`.
#'
#' @param t_obs POSIXct vector, the observed times in the time series
#' @param freq_year integer, the number of observations per year (365.25 days)
#' @param intercept logical indicating to include an intercept column (of `1`s)
#' @param daily_n integer, the number of fourier term pairs to use for daily cycles
#' @param yearly_n integer, the number of fourier term pairs to use for yearly cycles
#' @param default_n integer, copied to `daily_n` and/or `yearly_n` when they are NA
#' @param na_rm logical, omits rows from the output where the corresponding time is unobserved
#'
#' @return a matrix of covariates
#' @export
time_X = function(t_obs,
                  intercept = TRUE,
                  na_rm = FALSE,
                  knots_t = NULL,
                  knots_n = 5L,
                  daily_n = 5L,
                  yearly_n = 5L) {

  # unsorted times are invalid
  if( !identical(t_obs, sort(t_obs)) ) stop('expected t_obs to be sorted')

  # pad with NAs to make data frame with regular time steps covering t_obs
  t_df = data.frame(posix_pred=t_obs) |> archive_pad(quiet=TRUE)
  t_pad = t_df[['posix_pred']]
  t_len = nrow(t_df)
  is_obs = !is.na(t_df[['ts_hours']])

  # convert times to integer
  t_int = as.integer(t_pad)
  day_in_sec = 24*60*60
  t_len_sec = t_pad |> range() |> diff()

  # loop over periods to make fourier columns
  p_nm = stats::setNames(nm=c('day', 'year'))
  p = c(day = 60 * 60 * 24, year = 60 * 60 * 24 * 365.25)
  n = c(day = daily_n, year = yearly_n)
  X_fourier_list = p_nm |> lapply(\(nm) {

    # times modulo period (in seconds) and normalize to [0, 1]
    t_norm = as.integer(t_pad) %% p[nm] / p[nm]
    if(n[nm] == 0) { NULL } else {

      # loop to get covariates for each order
      do.call(cbind, lapply(seq(n[nm]), \(ord) {

        X_n = c(sin( 2 * pi * t_norm * ord),
                cos( 2 * pi * t_norm * ord)) |> matrix(ncol=2)

        # S = sin, C = cos, integer order is number of cycles per period
        colnames(X_n) = paste0(nm, '_', c('s', 'c'), ord)
        X_n

      }) )
    }
  })

  # make a spline basis if requested
  X_spline = NULL
  if( !any(is.na(knots_t)) ) {

    # set knots_t (if not supplied) based on input range
    if( is.null(knots_t) ) {

      n_inner = knots_n * ceiling( t_len_sec / (day_in_sec * 365.25) )
      knots_t = t_pad[ seq(1, t_len, length = 2 + n_inner) ] |> unique()
    }

    # assume UTC time zone for character input
    if( is.character(knots_t) ) knots_t = as.POSIXct(knots_t, tz='UTC')

    # splines::ns has boundary knots argument
    knots_inner = knots_t |> as.integer() |> head(-1) |> tail(-1)
    knots_outer = knots_t |> as.integer() |> range()
    X_spline = t_int |> splines::ns(knots=knots_inner, Boundary.knots=knots_outer)
    colnames(X_spline) = paste0('spline_', seq(ncol(X_spline)))
  }

  # make an intercept column if requested
  X_intercept = NULL
  if( intercept ) X_intercept = data.frame(intercept = rep(1, t_len)) |> as.matrix()

  # combine everything and optionally remove unobserved rows
  X_out = cbind(X_intercept, X_spline, do.call(cbind, X_fourier_list))
  if( na_rm ) X_out = X_out[is_obs, ]

  # NOTE na_rm doesn't change the attributes of the output, only the rows
  # returned. You still get the padded time series info in the attributes

  # add attributes mapping rows to input times before returning
  attr(X_out, 'knots') = knots_t
  attr(X_out, 'time') = t_pad
  attr(X_out, 'idx_obs') = which(is_obs)
  attr(X_out, 'na_map') = seq(t_len) |> match(attr(X_out, 'idx_obs'))
  return(X_out)
}
