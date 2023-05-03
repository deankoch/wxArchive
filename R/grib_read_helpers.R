

#' Find layer and sub-grid indices for calls to `terra::rast` and `terra::values`
#'
#' This uses `terra` to read the layer names and grid dimensions from the GRIB at
#' `grib_path` (but not its data values) in order to learn (1) the index of layer
#' names matching the regular expressions in `regex`, and (2) the dimensions and
#' position of the sub-grid inscribing the area of interest `aoi`. These vectors can
#' then be passed to `terra::rast` and `terra::values` to speed up loading the data
#' (see `my_grib_read`).
#'
#' The function returns a list with the following
#'
#'  * 'k' = index of requested layers
#'  * 'ij' = row/column of top-left corner point in AOI sub-grid
#'  * 'dim' = dimensions of AOI sub-grid
#'  * 'r_aoi' = empty SpatRaster defining the AOI
#'
#' If `regex=NULL`, all layers are returned. If any of the `regex` match more than
#' one layer, only the first match is returned. If `aoi=NULL`, the entire grid is
#' requested. Otherwise `aoi` should be coercible to 'SpatVector' and must overlap
#' with some part of the grid in `grib_path` (after transforming to a common CRS).
#'
#' `grib_path` can be a data-frame containing multiple paths, in which case the first
#' is loaded. If this load attempt fails, and `try_again=TRUE`, the function will then
#' attempt to load the second file, and so on.
#'
#' @param grib_path character (or dataframe with 'path' column). The path to the GRIB file
#' @param regex character vector, regular expressions to match to layer names (using `grepl`)
#' @param aoi geometry object coercible to `SpatVector`, the area of interest
#' @param quiet logical indicating to suppress console output
#' @param try_again logical indicating to attempt to load another file if the first fails
#'
#' @return list of indices (see details)
#' @export
my_grib_idx = function(grib_path, regex=NULL, aoi=NULL, quiet=FALSE, try_again=FALSE) {

  # by default pick first file
  if(is.data.frame(grib_path)) grib_path = grib_path[['path']]
  grib_path = grib_path[1]

  if(!quiet) cat('\nchecking names in ', basename(grib_path))

  # attempt to load the file
  r = tryCatch({ grib_path |> terra::rast() |> terra::rast() }, error = function(err) err)
  if(is(r, 'error')) {

    if( !try_again | (length(grib_path) < 2) ) stop('failed to load ', grib_path)

    # on failures, attempt to load next file in the list
    return( my_grib_idx(tail(grib_path, -1), regex, aoi, quiet, TRUE) )
  }

  # identify layer(s) of interest
  k = names(r) |> which_lyr(regex, quiet=quiet, na.rm=TRUE)

  # default aoi is the whole grid
  if( is.null(aoi) ) { r_aoi = r } else {

    aoi_projected = as(aoi, 'SpatVector') |> terra::project(terra::crs(r))
    r_aoi = terra::crop(terra::rast(r), aoi_projected)
  }

  # find the index of first column and row of aoi in outer grid
  i = terra::rowFromY(r, terra::yFromRow(r_aoi, 1))
  j = terra::colFromX(r, terra::xFromCol(r_aoi, 1))

  # strip the inner grid of values and unneeded layers
  r_aoi = terra::rast(r_aoi[[1]])

  # count number of cells along each dimension then return
  return( list(k=k, ij=c(i, j), dim=dim(r_aoi)[1:2], r_aoi=r_aoi) )
}


#' GRIB loader to catch errors with missing/damaged files
#'
#' This is a wrapper for `terra::rast` that loads only the sub-grid and
#' layers specified in `regex` and `aoi` (using `my_grib_idx`) for the GRIB
#' file at `grib_path`.
#'
#' This catches loading errors - instead of halting, the function
#' returns the error object (and prints it as a warning if `quiet=FALSE`).
#'
#' @param grib_path character path to the GRIB file to (attempt to) load
#' @param regex character vector passed to `my_grib_idx` (layer names)
#' @param aoi geometry object passed to `my_grib_idx` (area of interest)
#' @param quiet logical indicating to suppress console output
#'
#' @return SpatRaster containing the requested layers and sub-grid
#' @export
my_grib_read = function(grib_path, regex=NULL, aoi=NULL, quiet=FALSE) {

  # attempt to load the file
  t_mat = tryCatch({

    # find indices to load
    r_info = my_grib_idx(grib_path, regex=regex, aoi=aoi, quiet=quiet)

    # copies data as matrix
    terra::rast(grib_path, lyrs = r_info[['k']]) |>
      terra::values(row = r_info[['ij']][1],
                    nrows = r_info[['dim']][1],
                    col = r_info[['ij']][2],
                    ncols = r_info[['dim']][2])

  }, error = function(err) err)

  if(is(t_mat, 'error') & !quiet) {

    message('\nfailed to load grib ', t_path, '\n')
    warning(as.character(t_mat))
  }

  return(t_mat)
}


#' Extract the specified area of interest (sub-grid) and variables (layers) from GRIB archive
#'
#' This opens (the subset `file_idx` of) GRIB files listed in `grib_df`, loading the
#' sub-grid and layers specified by `aoi` and `regex` (respectively). Results are returned
#' as a list of SpatRasters, one for each variable in `regex`, each with a layer for every
#' file (ie time).
#'
#' With a large number of GRIBs, memory usage can get high and everything slows down
#' overall. I recommend loading at most 100s to 1000s of files at once, and merging results
#' as you go. `memory_limit` sets up a limit (in GB) as a sanity check.
#'
#' The function calls `my_grib_idx` on `grib_df` (and not the subset in `file_idx`)
#' to initially get basic grid information like dimensions, so make sure all the GRIBs
#' in `grib_df` use the same grid.
#'
#' @param grib_df data frame, a subset of the rows by `grib_list`
#' @param file_idx integer vector, the rows of `grib_df` to load
#' @param regex character vector passed to `my_grib_idx` (layer names)
#' @param aoi geometry object coercible to `SpatVector`, the area of interest
#' @param memory_limit integer (in GB) the maximum memory usage allowed of the output
#'
#' @return a list of multi-layer SpatRasters (one per variable in `regex`)
#' @export
my_grib_extract = function(grib_df, file_idx=NULL, regex=.rr_regex, aoi=NULL, memory_limit=8L) {

  # by default open all files in grib_df
  if(is.null(file_idx)) file_idx = nrow(grib_df) |> seq()

  # get grid dimensions by looking at first file
  r_info = my_grib_idx(grib_df, regex, aoi=aoi, quiet=TRUE, try_again=TRUE)

  # sanity check for memory usage
  n_space = prod(r_info[['dim']])
  n_var = length(regex)
  n_file = length(file_idx)
  est_gb = estimate_memory(r_info[['dim']], n_var, n_file)
  if(est_gb > memory_limit) stop('memory limit exceeded.')

  # initialize storage matrices for numeric data from all times
  z_out = lapply(regex, \(x) matrix(NA_real_, n_space, n_file))

  # loop over files, copying only the requested variables to storage
  cat('\nloading', n_file, 'files...\n')
  if(n_file > 1) pb = txtProgressBar(max=n_file, style=3)
  for(t in seq(n_file)) {

    # progress bar in console updated with current filename
    if(n_file > 1) setTxtProgressBar(pb, t)
    t_path = grib_df[['path']][file_idx[t]]
    cat('\r|=', basename(t_path), '')

    # read in the requested layers, handling load failures
    t_mat = my_grib_read(t_path, regex=regex, aoi=aoi, quiet=TRUE)

    # copy all data in loop over variables
    if( !is(t_mat, 'error') ) {

      # loop over raster layer names, filling columns of output matrix
      j_out = colnames(t_mat) |> which_lyr(regex, na.rm=FALSE, quiet=TRUE)
      for(nm in names(j_out)) z_out[[nm]][, t] = t_mat[, j_out[nm]]
    }
  }
  if(n_file > 1) close(pb)

  # make a template raster for output
  r_template = my_grib_idx(grib_df, aoi=aoi, try_again=TRUE, quiet=TRUE)[['r_aoi']] |>
    terra::rast(nlyrs=length(regex))

  # initialize SpatRaster list in memory, split by variable name
  r_out = r_template |>
    lapply(\(r) terra::rast(r, nlyrs=n_file)) |>
    stats::setNames(nm=names(regex))

  # loop over variables, copying matrix data to SpatRaster
  idx_write = seq(n_space)
  cat('\ncopying', n_file, 'layers to', paste(names(regex), collapse=', '), '... ')
  for(v in seq_along(r_out)) {

    # assign timestamps then data values in a loop over layers (times)
    terra::time(r_out[[v]]) = grib_df[['posix_pred']][file_idx]
    names(r_out[[v]]) = terra::time(r_out[[v]])

    # assign values to raster cells
    terra::set.values(r_out[[v]], cells = idx_write, values = z_out[[v]][])
  }

  cat('done.\n')
  return(r_out)
}
