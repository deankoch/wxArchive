#' Create/update a netCDF file with missing pcp_total times filled using pcp_large and pcp_small
#'
#' This creates a copy of the "pcp_total" time series, named `pcp_nm`, containing only
#' times imputed using the sum of "pcp_large" (large scale precipitation) and "pcp_small"
#' (convective precipitation) whenever both are available, but "pcp_total" is not.
#'
#' `input_nm` should be a named list character vector elements 'coarse' and/or 'fine', each
#' representing a batch of files to read a common variable from each resolution. Typical usage
#' has at least two source files at each resolution - a long term storage file, and a smaller file
#' with more recent times.
#'
#' The function writes one output file for each of the vectors in `input_nm`, to the
#' location(s) in `output_nm`. The output file name is `paste0(pcp_nm, 'nc')`, so the
#' function will not allow you to set `pcp_nm` to any of the (in-use) names to avoid
#' confusion downstream (list them with `file_wx('nc', base_dir, output_nm)`)
#'
#' Layers with negatives are set to `NA` in the output file.
#'
#' A JSON for the output is written to the "time" subdirectory. The function will update
#' existing files by only copying layers for times not already listed in the file.
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param pcp_nm a name for the output variable and file
#' @param input_nm list of character vectors, naming input subdirectories in `base_dir`
#' @param output_nm character vectors, naming an output subdirectory for each `input_nm` vector
#'
#' @return nothing, but possibly modifies the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
my_pcp_total = function(base_dir,
                        pcp_nm = 'pcp',
                        input_nm = list(coarse='coarse', fine='fine'),
                        output_nm = c('coarse', 'fine')) {

  # paths to expected inputs
  pcp_var_nm = list('pcp_total', 'pcp_large', 'pcp_small')
  input_path = file_wx('nc', base_dir, input_nm, var_nm=pcp_var_nm)

  # sanity check for pcp_nm
  nm_var_all = stats::setNames(nm=names(input_path[[1]]))
  if( pcp_nm %in% nm_var_all ) stop('pcp_nm cannot be an existing variable name')

  # new files to write
  output_nc_path = file_wx('nc', base_dir, as.list(output_nm), var_nm=pcp_nm)

  # find time coverage of each variable at both resolutions
  cat('\nscanning available times for', paste(nm_var_all, collapse=', '))
  var_info = lapply(input_path, \(r) lapply(r, \(p) my_nc_attributes(p, ch=TRUE)) )

  # loop over resolutions
  for( nm_res in output_nm ) {

    cat('\n\nchecking', nm_res, 'grids')
    t1 = proc.time()

    # check for available times by variable
    t_all = do.call(c, lapply(var_info[[nm_res]], \(v) v[['time']])) |> unique() |> sort()
    is_obs = var_info[[nm_res]] |> lapply(\(v) t_all %in% v[['time_obs']])
    is_relevant = is_obs[['pcp_large']] & is_obs[['pcp_small']] & !is_obs[['pcp_total']]
    time_relevant = t_all[is_relevant]

    # check for existing data (creating JSON as needed)
    p_out = output_nc_path[[nm_res]]
    time_done = character(0) |> as.POSIXct()
    if( file.exists(p_out) ) time_done = my_nc_attributes(p_out, overwrite=TRUE, lazy=TRUE)[['time']]

    # select required layers from each source
    is_needed = !( time_relevant %in% time_done )
    time_add = time_relevant[is_needed]

    # if there's nothing to add then we are finished
    if( length(time_add) == 0 ) {

      cat('\nup to date \U2713')
      next
    }

    # compute totals from large and small source
    r_add = NULL
    if( length(time_add) > 0 ) {

      cat('\ncopying sum of pcp_large and pcp_small')
      r_large = input_path[[nm_res]][['pcp_large']] |> nc_layers(time_add)
      r_small = input_path[[nm_res]][['pcp_small']] |> nc_layers(time_add)
      r_add = r_large + r_small
    }

    # handle NAs encoded as negatives
    cat('\nchecking for negatives... ')
    is_na = c(terra::global(r_add, 'min')[[1]]) < 0
    n_na = sum(is_na)
    if( n_na > 0 ) {

      r_add[seq(nrow(r_add)), seq(ncol(r_add)), which(is_na)] = NA
      cat('omitting', n_na, 'NA layers')

    } else { cat('none found') }

    # append to existing data file (or create the file and write to it)
    r_add |> nc_write(p_out, overwrite=TRUE, append=TRUE)
    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  }
}


#' Create/update a NetCDF archive of resampled coarse resolution grids
#'
#' This resamples the coarse resolution time series to get data at fine resolution,
#' writing results to nc files in subdirectory `output_nm`. Only those times found in
#' the coarse series but not the fine series are included.
#'
#' Specify the coarse and fine sub-directories as named elements in list `input_nm`. Each
#' should be a character vector of sub-directories to check for data on `var_nm`. When
#' multiple sub-directories contain data for a given time, the function extracts the data
#' from the first (according to the order in `input_nm`)
#'
#' `var_nm` can be a list of character vectors, where each vector gives a set of variable
#' names to be treated as equivalent (eg 'pcp' and 'pcp_total', which are the same variable,
#' compiled from different sources). The output variable is assigned the corresponding name
#' from `var_nm`. If `var_nm` is unnamed, the function names it using the first element of
#' each vector.
#'
#' The sub-directory `output_nm` will be created if it doesn't exist. If the file in
#' `output_nm` exists already, only those times available in `input_nm`, but not already
#' present in `output_nm`, are appended.
#'
#' Coarse resolution grids are by default resampled by bilinear averaging (with GDAL).
#' Change this by setting named arguments in `...` to pass to `terra::project`.
#'
#' @param var_nm character vector or list of them, names of the variables to process
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm list of character vectors, namess of the source file subdirectories in `base_dir`
#' @param output_nm character name for the sub-directory to create.
#' @param SpatRaster optional target grid for resampling
#' @param ... additional arguments passed to `terra::project`
#'
#' @return nothing, but possibly modifies the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
my_resample = function(var_nm,
                       base_dir,
                       input_nm = list(coarse='coarse', fine='fine'),
                       output_nm = 'coarse_resampled',
                       r_fine = NULL,
                       from = NULL,
                       append = TRUE,
                       ...) {

  # paths to expected inputs and outputs
  input_nc = file_wx('nc', base_dir, as.list(input_nm), as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc[[1]]))
  output_nc = file_wx('nc', base_dir, output_nm, as.list(names(var_nm)), make_dir=TRUE)

  # find time coverage of each variable at both resolutions
  cat('\nscanning available times for', paste(names(var_nm), collapse=', '))
  var_info = lapply(input_nc, \(r) lapply(r, \(p) my_nc_attributes(p, ch=TRUE)) )

  # set default start times and check consistency of two arguments
  if( is.null(from) ) from = var_info[['coarse']] |> sapply(\(v) min(v[['time_obs']]) )
  if( length(from) == 1 ) from = rep(from, length(var_nm))
  if( length(from) != length(var_nm) ) stop('"from" must have the same length as "var_nm" (or 1)')
  names(from) = names(var_nm)

  # loop over variable names to check for existing resampled data
  time_add = stats::setNames(nm=names(var_nm)) |> lapply(\(nm) {

    # copy existing output times (or nothing if overwriting)
    time_done = if( !append | !file.exists(output_nc[[nm]]) ) as.POSIXct(integer(0)) else {

      #open JSON or create it as needed
      my_nc_attributes(output_nc[[nm]], overwrite=TRUE, lazy=TRUE)[['time_obs']]
    }

    # find times in coarse series that haven't been processed yet
    time_coarse = var_info[['coarse']][[nm]][['time_obs']]
    time_fine = var_info[['fine']][[nm]][['time_obs']]
    time_pending = time_coarse[ !( time_coarse %in% c(time_fine, time_done) ) ]

    # filter to supplied start time (default includes all)
    time_pending = time_pending[ time_pending >= from[nm] ]
  })

  # finished if no updates to write
  needs_update = sapply(time_add, length) > 0
  if( !any(needs_update) ) {

    cat('\nall variables are up to date \U2713')
    return(invisible())
  }

  # a template SpatRaster at fine resolution from first file
  if( is.null(r_fine) ) r_fine = input_nc[['fine']][[1]][1] |> terra::rast(lyrs=1) |> terra::rast()

  # loop over variables
  for(nm in names(var_nm)) {

    # skip variables that are up to date
    if( needs_update[nm] ) {

      cat('\n\nprocessing', nm, '...')
      t1 = proc.time()

      # `project` calls GDAL
      cat('\nresampling', length(time_add[[nm]]), 'coarse resolution layers')
      r_resample = input_nc[['coarse']][[nm]] |>
        nc_layers(time_add[[nm]]) |>
        terra::project(r_fine)
        #terra::project(r_fine, ...)

      # check for problem layers
      n_na = r_resample |> terra::global('isNA') |> as.matrix() |> as.numeric()
      if( any(n_na == terra::ncell(r_resample)) ) stop('one or more layers had all NA cells')

      # clean up NAs around borders with moving window mean
      is_na = n_na > 0
      sanity_i = 1
      while( any(is_na) & (sanity_i < 1e3) ) {

        r_resample[[which(is_na)]] = r_resample[[which(is_na)]] |>
          terra::focal(fun=mean, na.rm=TRUE, na.policy='only')

        # update index of NA layers
        n_na[is_na] = r_resample[[is_na]] |> terra::global('isNA') |> as.matrix() |> as.numeric()
        is_na = n_na > 0
        sanity_i = sanity_i + 1
      }

      # append to existing data file (or create the file and write to it)
      r_resample |> nc_write(output_nc[[nm]], overwrite=TRUE, append=append)

      t2 = proc.time()
      cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')

    } else { cat(paste0('\n', nm), 'is up to date \U2713') }
  }
}



