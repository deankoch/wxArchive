#' Create/update a NetCDF archive of resampled coarse resolution grids
#'
#' This resamples the coarse resolution time series to get data at fine resolution,
#' writing results to NetCDF files in subdirectory `output_nm`. Only times found in
#' the coarse series but not the fine series are processed
#'
#' Specify the coarse and fine sub-directories as named elements in list `input_nm`. Each
#' should be a character vector of sub-directories. These are checked for time series files
#' named `var_nm` (see `?nc_layers`)
#'
#' `var_nm` can be a list of character vectors, where each vector gives a set of variable
#' names to be treated as equivalent (eg 'pcp' and 'pcp_total'). The output variables are
#' assigned names from `names(var_nm)`.
#'
#' The function creates sub-directory `output_nm` if it doesn't exist already. Existing
#' files in `output_nm` are updated by appending all times in `input_nm` occurring after
#' the last time in `output_nm`.
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
nc_resample = function(var_nm,
                       base_dir,
                       input_nm = list(coarse='coarse', fine='fine'),
                       output_nm = 'coarse_resampled',
                       r_fine = NULL,
                       from = NULL,
                       ...) {

  # paths to expected inputs and outputs
  input_nc = file_wx('nc', base_dir, as.list(input_nm), as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc[[1]]))
  output_nc = file_wx('nc', base_dir, output_nm, as.list(names(var_nm)), make_dir=TRUE)

  # get a template SpatRaster at fine resolution from first available file
  if( is.null(r_fine) ) {

    r_fine_path = input_nc[['fine']][[1]][file.exists(input_nc[['fine']][[1]])]
    if( length(r_fine_path) == 0 ) stop('no files found in "fine" sub-directory. Try supplying r_fine')
    r_fine = r_fine_path |> terra::rast(lyrs=1) |> terra::rast()
  }

  # find time coverage of each input variable at both resolutions
  cat('checking available times for', paste(names(var_nm), collapse=', '))
  var_info = lapply(input_nc, \(r) lapply(r, \(p) time_wx(p)) )

  # loop over variable names to check for existing resampled data
  time_add = stats::setNames(nm=names(var_nm)) |> lapply(\(nm) {

    # copy existing output times (or nothing if overwriting)
    time_done = if( file.exists(output_nc[[nm]]) ) {

      # open JSON or create it as needed
      time_wx(output_nc[[nm]])[['time']]

    } else { NULL }

    # find times in coarse series that haven't been processed yet
    time_coarse = var_info[['coarse']][[nm]][['time_obs']]
    time_fine = var_info[['fine']][[nm]][['time_obs']]
    time_pending = time_coarse[ !( time_coarse %in% c(time_fine, time_done) ) ]

    # filter to supplied start time (default includes all)
    from_i = if( is.null(from[nm]) ) min(time_coarse) else as.POSIXct(from[nm], tz='UTC')
    time_pending = time_pending[ time_pending >= from_i ]
  })

  # finished if no updates to write
  needs_update = sapply(time_add, length) > 0
  if( !any(needs_update) ) {

    cat('\nall variables are up to date\n')
    return(invisible())
  }

  # loop over variables
  for(nm in names(var_nm)) {

    # skip variables that are up to date
    if( needs_update[nm] ) {

      cat('\n\nprocessing', nm, '...')
      t1 = proc.time()

      # append results to existing data file (or create the file and write to it)
      input_nc[['coarse']][[nm]] |>
        nc_project(r_fine, times=time_add[[nm]], ...) |>
        nc_write(output_nc[[nm]])

      t2 = proc.time()
      cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')

    } else { cat(paste0('\n', nm), 'is up to date\n') }
  }
}


#' Wrapper for terra::project with gap-filling
#'
#' This loads the NetCDF data in `nc_path` and resamples to match the grid
#' in `target` using `terra::project`. NAs are then imputed using with a moving
#' window average using `terra::focal`.
#'
#' This addresses a problem where two spatial datasets are cropped to the same
#' bounding box in different projections. Resampling one to match the other can
#' result in NA pixels around edges of your grid. This can be avoided in the first
#' place by expanding the bounding box of the source grid, but this isn't always
#' easy or practical.
#'
#' After projecting, the function calls `terra::focal`, with default settings
#' and `fun=mean`, repeatedly until all NAs are filled (to a maximum of 1e3
#' iterations).
#'
#' Only the times listed in `times` are loaded and processed. Set `times=NULL`
#' to process all. Both `target` and the file(s) at `nc_path` must have a well
#' defined CRS (eg check `terra::crs(target)`).
#'
#' @param nc_path character vector, the path(s) to the input NetCDF file(s)
#' @param target a SpatRaster
#' @param times vector of POSIXct times to load from `nc_path`
#' @param ... named arguments to pass to `terra::project`
#'
#' @return SpatRaster with grid matching `target` and data layers resampled from `nc_path`
#' @export
nc_project = function(nc_path, target, times=NULL, ...) {

  # load all data into memory and project with GDAL via terra (fast)
  r_resampled = nc_path |> nc_layers(times) |> terra::project(target, ...)

  # check for problem layers
  not_na = r_resampled |> terra::global('notNA') |> as.matrix() |> as.numeric()
  if( any(not_na == 0) ) stop('one or more layers had all NA cells')

  # clean up NAs around borders with moving window mean
  cat('\ncleaning edges...')
  is_na = not_na < terra::ncell(r_resampled)
  insanity = 1
  while( any(is_na) & (insanity < 1e3) ) {

    # resample then update index of NA layers
    r_resampled = r_resampled |> terra::focal(fun=mean, na.rm=TRUE, na.policy='only')
    not_na = r_resampled |> terra::global('notNA') |> as.matrix() |> as.numeric()
    is_na = not_na < terra::ncell(r_resampled)
    insanity = insanity + 1
  }

  return(r_resampled)
}


