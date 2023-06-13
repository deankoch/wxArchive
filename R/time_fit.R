#' Fit a temporal model to a NetCDF time series
#'
#' This loops over `var_nm`, loading the time series data from disk in sub-folder(s)
#' `input_nm` of `base_dir`, and a temporal model to it. The model first regresses on
#' a matrix of seasonal predictors (see `?time_X`) by OLS before fitting an AR(2) model
#' to the residuals.
#'
#' The model is fitted independently to each spatial grid point using all time points
#' and a common covariates matrix. Fitted parameters are saved as raster layers in a NetCDF
#' file named using `file_wx('temporal_nc', base_dir, model_nm, ...)`. This file name
#'  - which changes every time `file_wx` is called - is then appended to the JSON file
#' `file_wx('temporal_index', base_dir, model_nm, ...)`. This JSON keeps track of all
#' previous model fits.
#'
#' Arguments `daily_n`, `yearly_n`, and `knots_n` are passed to `time_X` to make the
#' seasonal predictors. These are also copied to the JSON, along with information about
#' the knot locations. This information should all be loaded and passed to `time_X` when
#' predicting on new times.
#'
#' @param var_nm list of character vectors, the variable names to fit
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm character vector or list, subdirectories containing the .nc files to fit
#' @param model_dir character, name of parent directory of `model_nm`
#' @param model_nm character vector or list, sub-directories to write output files
#' @param daily_n integer number of Fourier term pairs for daily cycles (passed to `time_X`)
#' @param yearly_n integer number of Fourier term pairs for yearly cycles (passed to `time_X`)
#' @param knots_n integer number of knots to use for spline (passed to `time_X`)
#'
#' @return returns nothing but modifies the files in "model" subdirectory of "input_nm"
#' @export
time_fit = function(var_nm,
                    base_dir,
                    input_nm = 'fine',
                    model_dir = base_dir,
                    model_nm = .nm_temporal_model,
                    daily_n = 5L,
                    yearly_n = 5L,
                    knots_n = 5L) {

  # input/output paths (var_nm is list to ensure output paths are in list)
  input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_json = file_wx('temporal_index', model_dir, model_nm, var_nm_list, make_dir=TRUE)
  output_nc = file_wx('temporal_nc', model_dir, model_nm, var_nm_list, make_dir=TRUE)

  # load time coverage of each variable
  cat('\nreading times and grid information for', paste(names(var_nm), collapse=', '))
  var_info = input_nc |> lapply(\(p) time_wx(p))

  # loop over variables
  for(v in seq_along(var_nm)) {

    nm = names(var_nm)[v]
    cat('\n\nprocessing', nm, '...\n')
    t1 = proc.time()

    # add POSIXct time of function call to result
    call_time = Sys.time() |> as.character(tz='UTC')

    # extract observed times, compute step size
    t_all = var_info[[nm]][['time_obs']]
    step_hours = data.frame(posix_pred=t_all) |> time_step()
    cat('\nstep size:', step_hours, 'hours')
    freq_year = as.integer(24 * 365.25 / step_hours)

    # sample all observed times
    t_fit = t_all
    n_fit = length(t_fit)

    # copy data to memory as sk object
    cat('\nselected subset', as.character(min(t_fit)), 'to', as.character(max(t_fit)))
    cat('\ncopying', n_fit, 'layers to memory')
    g_obs = input_nc[[nm]] |> nc_layers(t_fit, na_rm=TRUE) |> snapKrig::sk()

    # build matrix of seasonal predictors and compute QR decomposition
    cat('\ncomputing temporal covariates')
    X_time = time_X(t_fit, daily_n=daily_n, yearly_n=yearly_n, knots_n=knots_n, na_rm=FALSE)
    X_time_obs_qr = qr(X_time[attr(X_time, 'idx_obs'),])

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
    g_obs |> snapKrig::sk(gval=fit_mat) |> snapKrig::sk_export() |> terra::writeCDF(output_nc[[nm]])

    # remove temporary raster data from memory
    rm(g_obs)
    gc()

    # load any existing model fitting info and compile the new addition
    json_exists = file.exists(output_json[[nm]])
    json_old = if(!json_exists) NULL else output_json[[nm]] |> readLines() |> jsonlite::fromJSON()
    json_add = list(var_nm = var_nm[[v]],
                    sub_dir = input_nm,
                    at = call_time,
                    file = basename(output_nc[[nm]]),
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
  }
}
