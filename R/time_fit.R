#' Fit a temporal model to a NetCDF time series
#'
#' @param var_nm list of character vectors, the variable names to fit
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm character vector or list, subdirectories containing the .nc files to fit
#' @param p_max numeric the maximum proportion of missing points in the series
#' @param n_max integer maximum number of points to sample (NA to select all)
#' @param daily_n integer number of Fourier term pairs for daily cycles
#' @param yearly_n integer number of Fourier term pairs for yearly cycles
#' @param knots_n integer number of knots to use for spline
#'
#' @return returns nothing but modifies the files in "model" subdirectory of "input_nm"
#' @export
time_fit = function(var_nm,
                    base_dir,
                    input_nm = 'fine',
                    daily_n = 5L,
                    yearly_n = 5L,
                    knots_n = 5L,
                    p_max = 0.1,
                    n_max = NA) {

  # input/output paths (var_nm is list to ensure output paths are in list)
  input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_json = file_wx('temporal_index', base_dir, input_nm[1], var_nm_list, make_dir=TRUE)
  output_nc = file_wx('temporal_nc', base_dir, input_nm[1], var_nm_list, make_dir=TRUE)

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
    g_obs = input_nc[[nm]] |> nc_layers(t_fit) |> snapKrig::sk()

    # build matrix of seasonal predictors and compute QR decomposition
    cat('\ncomputing temporal covariates')
    X_time = time_X(t_fit, daily_n=daily_n, yearly_n=yearly_n, knots_n=knots_n, na_rm=FALSE)
    X_time_obs_qr = qr(X_time[attr(X_time, 'idx_obs'),])
    cat(' \U2713')

    # make storage matrix for results, then loop over pixels
    n_g = prod(dim(g_obs))
    fit_mat = matrix(NA_real_, n_g, 2L + ncol(X_time))
    cat('\nfitting seasonal trend and AR2 model separately to', n_g, 'time series...\n')
    pb = utils::txtProgressBar(max=n_g, style=3)
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
      utils::setTxtProgressBar(pb, k)
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
