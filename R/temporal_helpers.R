#' Fill temporal gaps using AR(2) model and seasonal predictor
#'
#' @param var_nm character vector or list of them, the input variables to process
#' @param base_dir 
#' @param input_nm 
#' @param output_nm 
#' @param model_nm 
#' @param n_max 
#' @param until 
#' @param quiet 
#'
#' @return
#' @export
#'
#' @examples
my_impute_temporal = function(var_nm,
                              base_dir,
                              input_nm = 'fine',
                              output_nm = 'completed',
                              model_nm = input_nm[1],
                              n_max = NULL,
                              until = NULL,
                              quiet = FALSE) {
  
  # data input/output paths (var_nm is list to ensure output paths are in list)
  input_nc = my_file_path('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_nc = my_file_path('nc', base_dir, output_nm, var_nm_list, make_dir=TRUE)
  
  # paths to fitted parameter nc files are in this file
  pars_json = my_file_path('temporal_index', base_dir, model_nm[1], var_nm_list)
  
  # load time coverage of each variable 
  cat('\nreading times and grid information for', paste(names(var_nm), collapse=', '))
  var_info = input_nc |> lapply(\(p) my_nc_attributes(p, ch=TRUE))

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
    
    # load latest fitted parameter matrices and take their element-wise mean
    pars_nc = pars_info |> sapply(\(x) x[['file']])
    
    #pars_r = Reduce('+', lapply(pars_nc, terra::rast)) / length(pars_info)
    pars_r = tail(pars_info, 1)[[1]][['file']] |> terra::rast()
    t_knots = tail(pars_info, 1)[[1]][['knots']] |> as.POSIXct(tz='UTC')
    gdim = dim( pars_r[[1]] )[1:2]
    
    # first two parameters are the AR2 coefficients
    pars_X = pars_r[][, -seq(2)]
    pars_AR2 = pars_r[][, seq(2)]
    
    # check for existing imputed times
    t_out = my_nc_attributes(output_nc[[nm]])[['time_obs']]
    
    # extract observed times, compute step size, compute default n_max (2 years)
    t_obs = var_info[[nm]][['time_obs']]
    step_hours = data.frame(posix_pred=t_obs) |> my_detect_step()
    freq_year = as.integer(24 * 365.25 / step_hours)
    if( is.null(n_max) ) n_max = 2 * freq_year

    # make data frame of times observed padded with NAs for gaps
    ts_df = data.frame(posix_pred=t_obs) |> 
      my_archive_padder(quiet=quiet, until=until) |>
      dplyr::mutate( missing = is.na(ts_hours) ) |>
      dplyr::mutate( done = posix_pred %in% t_out ) |>
      dplyr::mutate( changed = done & !missing )
    
    # no gaps found in input
    if( !any(ts_df[['missing']]) ) {
      
      cat('\ntime series is complete \U2713')
      next
    }
    
    # warn of times that were previously missing (and imputed) but are now observed
    if( any( ts_df[['changed']] ) ) {
      
      t_changed = ts_df[['posix_pred']][ ts_df[['updated']] ]
      msg_changed = t_changed |> paste(collapse=', ')
      output_dir = dirname(output_nc[[nm]])
      msg_info_1 = paste('\n(fix by deleting', output_dir, 'and running my_impute_temporal again)')
      cat(output_nc[[nm]], 'contains observed time(s):', msg_changed, msg_info_1)
    }
    
    # collect missing times marked for processing
    t_update = ts_df |> dplyr::filter(missing) |> dplyr::filter(!done) |> dplyr::pull(posix_pred)
    if( length(t_update) == 0 ) {
      
      cat('\nall missing times have been imputed already \U2713')
      next
    }
    
    # data frame of gaps to impute, filtered to time range requiring update
    t_start = min(t_update)
    gap_table = ts_df |> my_gap_finder(times='posix_pred') |> filter(end_time >= t_start)
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
    r_obs = input_nc[[nm]] |> my_nc_layers(t_pre)
    mat_both = r_obs[][, match(t_both, t_pre)]
    cat(' \U2713')
    
    # loop over gaps
    if(n_gap > 1) {
      
      cat('\nlooping over', n_gap, 'gaps to impute missing layers...\n')
      pb = txtProgressBar(max=n_gap, style=3)
    }
    for(g in seq(n_gap)) {
      
      # copy relevant times
      t_pred = t_all[ gap_seq[[g]] ]
      t_given = t_all[ gap_seq_pre[[g]] ]
      t_g = c(t_given, t_pred)
      
      # copy observed data to matrix with NAs for missing times
      z_obs = mat_both[, match(c(t_given, t_pred), t_both)]
      
      # build matrix of seasonal predictors for observed times
      X_time = t_g |> my_time_X(daily_n=daily_n, yearly_n=yearly_n, knots_t=t_knots)
      
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
      if(n_gap > 1) setTxtProgressBar(pb, g)
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
    r_pred |> my_nc_write(output_nc[[nm]], overwrite=TRUE, append=TRUE)
    t2 = proc.time()
    cat('\n\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.')
  }
}

#' Fit a temporal model to a NetCDF time series
#'
#' @param var_nm list of character vectors, the variable names to fit
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm character vector or list, subdirectories containing the .nc files to fit
#' @param p_max numeric the maximum proportion of missing points in the series 
#' @param n_max the maximum number of points to sample (NA to select all)
#'
#' @return
#' @export
my_fit_temporal = function(var_nm,
                           base_dir,
                           input_nm = 'fine',
                           daily_n = 5L,
                           yearly_n = 5L,
                           knots_n = 5L,
                           p_max = 0.1,
                           n_max = NA) {
  
  # input/output paths (var_nm is list to ensure output paths are in list)
  input_nc = my_file_path('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_json = my_file_path('temporal_index', base_dir, input_nm[1], var_nm_list, make_dir=TRUE)
  output_nc = my_file_path('temporal_nc', base_dir, input_nm[1], var_nm_list, make_dir=TRUE)
  
  # load time coverage of each variable 
  cat('\nreading times and grid information for', paste(names(var_nm), collapse=', '))
  var_info = input_nc |> lapply(\(p) my_nc_attributes(p, ch=TRUE))
  cat(' \U2713')
  
  # loop over variables
  for(v in seq_along(var_nm)) {
    
    nm = names(var_nm)[v]
    cat('\n\nprocessing', nm, '...\n')
    t1 = proc.time()
    
    # add POSIXct time of function call to result
    call_time = Sys.time() |> as.character(tz='UTC')
    
    # extract observed times, compute step size, set default n_max (all times)
    t_all = var_info[[nm]][['time_obs']]
    step_hours = data.frame(posix_pred=t_all) |> my_detect_step()
    cat('\nstep size:', step_hours, 'hours')
    freq_year = as.integer(24 * 365.25 / step_hours)
    if( is.na(n_max) ) n_max = length(t_all)

    # select a subsample of times at random
    n_fit = length(t_all) |> pmin(n_max)
    t_fit = my_sample_na(t_all, n_fit, p_max=p_max, step_hours=step_hours)
    
    # copy data to memory as sk object
    cat('\nselected subset', as.character(min(t_fit)), 'to', as.character(max(t_fit)))
    cat('\ncopying', n_fit, 'layers to memory')
    g_obs = input_nc[[nm]] |> my_nc_layers(t_fit) |> snapKrig::sk()

    # build matrix of seasonal predictors and compute QR decomposition
    cat('\ncomputing temporal covariates')
    X_time = my_time_X(t_fit, daily_n=daily_n, yearly_n=yearly_n, knots_n=knots_n, na_rm=FALSE)
    X_time_obs_qr = qr(X_time[attr(X_time, 'idx_obs'),])
    cat(' \U2713')
    
    # make storage matrix for results, then loop over pixels
    n_g = prod(dim(g_obs))
    fit_mat = matrix(NA_real_, n_g, 2L + ncol(X_time))
    cat('\nfitting seasonal trend and AR2 model separately to', n_g, 'time series...\n')
    pb = txtProgressBar(max=n_g, style=3)
    for(k in seq(n_g)) {
      
      # time series vector with NAs for missing times
      z_pad = c(g_obs[k, attr(X_time, 'na_map')]) |> stats::ts(frequency=freq_year)
      
      # compute coefficients from lm(z~X_time) directly using QR, then linear predictor
      temporal_betas = X_time_obs_qr |> solve.qr(z_pad[attr(X_time, 'idx_obs')])
      z_lm = c(X_time %*% temporal_betas) |> stats::ts(frequency=freq_year) 
      
      # attempt to fit AR(2) model to residuals (assign zeros on failure)
      temporal_alphas = tryCatch({ 
        
        AR2_result = (z_pad - z_lm) |> stats::arima(order=c(2L, 0L, 0L), include.mean=FALSE)
        AR2_result[['coef']] |> c() 
  
      }, error = function(err) c(0,0))
      
      # copy to storage
      fit_mat[k,] = c(temporal_alphas, temporal_betas)
      setTxtProgressBar(pb, k)
    }
    
    close(pb)
    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
    
    # write all parameter data to nc on disk
    cat('\nwriting results to', output_nc[[nm]])
    g_obs |> sk(gval=fit_mat) |> sk_export() |> terra::writeCDF(output_nc[[nm]])
    cat(' \U2713')
    
    # load any existing model fitting info and compile the new addition
    json_exists = file.exists(output_json[[nm]])
    json_old = if(!json_exists) NULL else output_json[[nm]] |> readLines() |> jsonlite::fromJSON()
    json_add = list(var_nm = var_nm[[v]],
                    sub_dir = input_nm,
                    at = call_time,
                    file = output_nc[[nm]], 
                    daily_n = daily_n,
                    yearly_n = yearly_n,
                    knots = attr(X_time, 'knots'),
                    train = list(start = min(t_fit),
                                 end = max(t_fit),
                                 length = n_fit,
                                 n = length(t_fit)))
    
    # append the training set and file info to JSON
    cat('\nupdating', output_json[[nm]])
    json_out = json_old |> c(list(json_add))
    names(json_out) = paste0('fit_', seq_along(json_out))
    json_out |> jsonlite::toJSON(pretty=TRUE) |> writeLines( output_json[[nm]] )
    cat(' \U2713')
    
  }
}

#' Return a data matrix of seasonal linear predictors for a time series
#' 
#' # loosely based on `forecast:::...fourier`
#' 
#' This creates a data matrix for daily and yearly cycles in a regular time series
#' observed at times `t_obs`. With `na_rm=TRUE`, the output has a row for each
#' `t_obs`. Otherwise it has a row for each time in the sequence of frequency
#' `freq_year` spanning `t_obs`.
#' 
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
my_time_X = function(t_obs,
                     intercept = TRUE,
                     na_rm = FALSE,
                     knots_t = NULL,
                     knots_n = 5L,
                     daily_n = 5L,
                     yearly_n = 5L) {
  
  # unsorted times are invalid
  if( !identical(t_obs, sort(t_obs)) ) stop('expected t_obs to be sorted')
  
  # pad with NAs to make data frame with regular time steps covering t_obs
  t_df = data.frame(posix_pred=t_obs) |> my_archive_padder(quiet=TRUE)
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


