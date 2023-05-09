#' Fill temporal gaps using AR(2) model and seasonal predictor
#'
#' @param var_nm character vector or list of them, the input variables to process
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm character vector or list, subdirectories containing the .nc files to process
#' @param output_nm character, sub-directory to write the output nc files
#' @param model_nm character, sub-directory to look for fitted model files
#' @param n_max integer, maximum number of time points to use (NULL for all)
#' @param until POSIXct time, the series is truncated/padded to this time point (passed to `archive_pad`)
#' @param quiet logical indicating to suppress console output
#'
#' @return nothing, but writes output to nc files in `output_nm`
#' @export
time_impute = function(var_nm,
                       base_dir,
                       input_nm = 'fine',
                       output_nm = 'completed',
                       model_nm = input_nm[1],
                       n_max = NULL,
                       until = NULL,
                       quiet = FALSE) {

  # data input/output paths (var_nm is list to ensure output paths are in list)
  input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_nc = file_wx('nc', base_dir, output_nm, var_nm_list, make_dir=TRUE)

  # paths to fitted parameter nc files are in this file
  pars_json = file_wx('temporal_index', base_dir, model_nm[1], var_nm_list)

  # load time coverage of each variable
  cat('reading times and grid information for', paste(names(var_nm), collapse=', '))
  var_info = input_nc |> lapply(\(p) time_wx(p))

  # loop over variable names
  for(nm in names(var_nm)) {

    cat('\n\nprocessing', nm, '...')
    t1 = proc.time()

    # get daily_n and yearly_n from JSON
    pars_info = pars_json[[nm]] |> readLines() |> jsonlite::fromJSON()
    daily_n = pars_info |> sapply(\(x) x[['daily_n']]) |> unique()
    yearly_n = pars_info |> sapply(\(x) x[['yearly_n']]) |> unique()
    if( ( length(daily_n) > 1 ) ) stop('inconsistent "daily_n"')
    if( ( length(yearly_n) > 1 ) ) stop('inconsistent "yearly_n"')

    # load latest fitted parameter matrices
    #pars_r = Reduce('+', lapply(pars_nc, terra::rast)) / length(pars_info)
    pars_nc = pars_info |> sapply(\(x) x[['file']])
    pars_r = tail(pars_info, 1)[[1]][['file']] |> terra::rast()
    t_knots = tail(pars_info, 1)[[1]][['knots']] |> as.POSIXct(tz='UTC')
    gdim = dim( pars_r[[1]] )[1:2]

    # first two parameters are the AR2 coefficients
    pars_X = pars_r[][, -seq(2)]
    pars_AR2 = pars_r[][, seq(2)]

    # check for existing imputed times
    t_out = time_wx(output_nc[[nm]])[['time_obs']]

    # extract observed times, compute step size, compute default n_max (2 years)
    t_obs = var_info[[nm]][['time_obs']]
    step_hours = data.frame(posix_pred=t_obs) |> time_step()
    freq_year = as.integer(24 * 365.25 / step_hours)
    if( is.null(n_max) ) n_max = 2 * freq_year

    # make data frame of times observed padded with NAs for gaps
    ts_df = data.frame(posix_pred=t_obs) |>
      archive_pad(quiet=quiet, until=until) |>
      dplyr::mutate( missing = is.na(ts_hours) ) |>
      dplyr::mutate( done = posix_pred %in% t_out ) |>
      dplyr::mutate( changed = done & !missing )

    # no gaps found in input
    if( !any(ts_df[['missing']]) ) {

      cat('\ntime series is complete')
      next
    }

    # warn of times that were previously missing (and imputed) but are now observed
    if( any( ts_df[['changed']] ) ) {

      output_dir = dirname(output_nc[[nm]])
      t_changed = ts_df[['posix_pred']][ ts_df[['changed']] ]
      msg_info = paste('\n(fix by deleting', output_dir, 'and running time_impute again)')
      cat('\nNOTE:', output_nc[[nm]], 'contains observed time(s):', msg_info)
    }

    # collect missing times marked for processing
    t_update = ts_df |> dplyr::filter(missing) |> dplyr::filter(!done) |> dplyr::pull(posix_pred)
    if( length(t_update) == 0 ) {

      cat('\nup to date')
      next
    }

    # data frame of gaps to impute, filtered to time range requiring update
    t_start = min(t_update)
    gap_table = ts_df |> gap_finder(times='posix_pred') |> dplyr::filter(end_time >= t_start)
    n_gap = nrow(gap_table)

    # list with index of all gap times, and another that includes preceding two observations
    gap_seq = Map(\(j, k) j:k, j=gap_table[['start']], k=gap_table[['end']])
    gap_seq_pre = Map(\(j, k) j:k, j=gap_table[['start']]-2L, k=gap_table[['start']]-1L)

    # times of gaps,times of observations preceding gaps, and both combined
    t_all = ts_df[['posix_pred']]
    t_gap = t_all[ unlist(gap_seq) ]
    t_pre = t_all[ unlist(gap_seq_pre) ]
    t_pre = t_pre[ !(t_pre %in% t_gap) ]
    t_both = c(t_pre, t_gap) |> sort()

    # load all observed data into a matrix with NAs for unobserved to be filled below
    cat('\nloading observed data into memory')
    r_obs = input_nc[[nm]] |> nc_layers(t_pre)
    mat_both = r_obs[][, match(t_both, t_pre)]

    # loop over gaps
    if(n_gap > 1) {

      cat('\nlooping over', n_gap, 'gaps to impute missing layers...\n')
      pb = utils::txtProgressBar(max=n_gap, style=3)
    }
    for(g in seq(n_gap)) {

      # copy relevant times
      t_pred = t_all[ gap_seq[[g]] ]
      t_given = t_all[ gap_seq_pre[[g]] ]
      t_g = c(t_given, t_pred)

      # copy observed data to matrix with NAs for missing times
      z_obs = mat_both[, match(c(t_given, t_pred), t_both)]

      # build matrix of seasonal predictors for observed times
      X_time = t_g |> time_X(daily_n=daily_n, yearly_n=yearly_n, knots_t=t_knots)

      # compute linear predictor for each observed time then residuals
      X_pred = nrow(pars_X) |> seq() |> sapply(\(j) X_time %*% c(pars_X[j,]) ) |> t()

      #X_pred = nrow(X_time) |> seq() |> sapply(\(j) pars_X %*% X_time[j,] )
      z_res = z_obs - X_pred

      # recurse forward to compute expected values, filling residuals matrix
      for( j in seq(1, length(t_pred)) ) {

        # sum of lag 1 and 2 terms
        z_res[, j+2] = ( pars_AR2[,1] * z_res[, j+1] ) + ( pars_AR2[,2] * z_res[, j] )
      }

      # add back linear predictor to get predictions, copy to storage
      mat_both[, match(t_pred, t_both)] = (z_res + X_pred)[, -seq(2)]
      if(n_gap > 1) utils::setTxtProgressBar(pb, g)
    }
    if(n_gap > 1) close(pb)

    # export stack to SpatRaster in memory
    cat('\nexporting', length(t_gap), 'imputed layer(s) to SpatRaster')
    r_pred = r_obs |> terra::rast(nlyrs=length(t_both))
    r_cells = terra::ncell(r_pred) |> prod() |> seq()
    terra::set.values(r_pred, cells=r_cells, values=mat_both)
    terra::time(r_pred) = t_both

    # omit observed times
    r_pred = r_pred[[ (t_both %in% t_gap) ]]
    t_overwrite = terra::time(r_pred)

    # append to nc file
    r_pred |> nc_write(output_nc[[nm]])
    t2 = proc.time()
    cat('\n\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.')
  }
}


#' Return a tibble of start/end indices for all NA sequences in a vector
#'
#' A helper function for time_impute that identifies the start and end points of all gaps
#' (contiguous stretches of `NA` values) in a vector, along with their lengths. If `invert=TRUE`
#' the function returns information about the non-`NA` segments.
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
#' @param x a vector with at least one `NA` value
#' @param start_only logical for internal use
#' @param times POSIXct vector or character giving a column name in `x`
#' @param invert logical indicating to return the inverted result
#'
#' @return a tibble with a row of information about each distinct gap
#' @export
gap_finder = function(x, start_only=FALSE, times=NULL, invert=FALSE) {

  if( is.data.frame(x) ) {

    # if x is a data frame then times must give a column name
    is_time = names(x) == times[1]
    if( is.null(times) | !any(is_time) ) stop('expected "times" to be a column name of x')
    if( ncol(x) < 2 ) stop('expected a data frame with 2 or more columns')

    # proceed with first column of data frame (after removing time column)
    times = x[[ which(is_time)[1] ]]
    x = x[[ which(!is_time)[1] ]]
    return( gap_finder(x, start_only=start_only, times=times, invert=invert) )
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
    i_end_inv = x_dummy |> rev() |> gap_finder(start_only=TRUE) |> rev()
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

