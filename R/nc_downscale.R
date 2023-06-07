#' Down-scale a spatio-temporal time series in a set of NetCDF files
#'
#' This is a wrapper for `.nc_downscale` that down-scales - ie predicts at finer
#' resolution - the gridded time series data in `input_nm`. By default it writes
#' a set of yearly NetCDF files in the sub-directory `output_nm` of `base_dir`,
#' with one set for each of the variables named in `var_nm`.
#'
#' Down-scaling happens by universal kriging using the fitted spatial covariance
#' parameters saved in the folder named in `model_nm` (generate this using
#' `space_fit`), and using covariates generated from the digital elevation model
#' in `dem` with values in metres.
#'
#' Argument `down` is the factor by which to reduce grid line spacing. This should
#' be an integer greater than 1. eg. `down=2` "doubles the resolution" by inserting
#' a new grid line in between each existing pair of adjacent ones. This produces an
#' output with approximately 4X (ie `down^2`) as many grid points as the input.
#'
#' The function expects `down` to be used consistently within a given output
#' subdirectory. If you want to change the output resolution, you will need to delete
#' any existing output or else specify a new `output_nm`.
#'
#' Use `from` and `to` to specify a date range to update (inclusive), or leave them `NULL`
#' to use a default range. The default range is meant to cover layers originating from GFS
#' (allowing them to be can be replaced by newly added RAP layers, or more recently released
#' GFS forecasts).
#'
#' The default for `to` is always the latest available date in the input. The default for
#' `from` is 10 days before the latest date found in the existing output files. If there
#' are no existing outputs, the default is set to the earliest available date in the input.
#'
#' Define an area of interest with `poly_out`. The function conditions its
#' predictions on a subset of the input grid covering `poly_out`, buffered by
#' the distance `edge_buffer` (in m). By default this is set to the diagonal length
#' of a single grid cell. Set this to a larger number, like the side length of your
#' domain, to ensure that all observations are considered in predictions, reducing
#' edge effects.
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param dem SpatRaster of elevation data at finer resolution than the input data
#' @param down integer, the down-scaling factor
#' @param input_nm character, the sub-directory name of input to process
#' @param model_nm character, the sub-directory name where the model files can be found
#' @param output_nm character, the sub-directory name for output
#' @param var_nm character vector, the name(s) of the variable to process
#' @param poly_out sfc object, polygon of interest
#' @param edge_buffer numeric, length in metres to buffer input data (see details)
#' @param from Date, the first date of the sequence to process
#' @param to Date, the last date of the sequence to process
#'
#' @return returns nothing but possibly modifies the NetCDF data in `output_nm`
#' @export
nc_downscale = function(base_dir,
                        dem,
                        down = 100,
                        input_nm = .nm_daily,
                        model_nm = .nm_model,
                        output_nm = .nm_down,
                        var_nm = .var_daily,
                        poly_out = NULL,
                        edge_buffer = NULL,
                        from = NULL,
                        to = NULL,
                        write_nc = TRUE) {

  # var_nm is list to ensure input/output paths are lists too
  input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_nc = file_wx('nc', base_dir, output_nm, var_nm_list, make_dir=TRUE)

  # check for available input and existing output dates
  input_var_info = lapply(input_nc, \(p) time_wx(p))
  output_var_info = lapply(output_nc, \(p) time_wx(p))

  # read fitted parameter values of existing models for the variable(s)
  pars_json = file_wx('spatial', base_dir, model_nm, var_nm_list, make_dir=TRUE)
  pars_json = pars_json[ sapply(pars_json, file.exists) ]
  if( length(pars_json) == 0 ) stop('model file(s) not found in ', model_nm)
  pars_all =  pars_json |> lapply(jsonlite::fromJSON)

  # get input grid info from first nc file
  r_grid_in = nc_chunk(input_nc[[1]][1])[1] |> terra::rast(lyrs=1)
  crs_in = terra::crs(r_grid_in)
  aoi_in = r_grid_in |> terra::ext() |> sf::st_bbox() |> sf::st_as_sfc()

  # apparently setting "unnamed" crs works this way but not with sf::st_as_sfc(crs=crs_in)?
  sf::st_crs(aoi_in) = crs_in

  # set default polygon (whole extent)
  if( is.null(poly_out) ) poly_out = aoi_in
  poly_out = poly_out |> sf::st_geometry() |> sf::st_transform(crs_in)

  # set default buffer size in m (length of grid cell diagonal in source)
  if( is.null(edge_buffer) ) edge_buffer = sum( terra::res(r_grid_in)^2 ) |> sqrt()
  cat('\nusing', round(edge_buffer/1e3, 2), 'km buffer to set bounding box for analysis')

  # transform output polygons to input projection and find their bounding box
  bbox_out = sf::st_bbox(poly_out) |> sf::st_as_sfc()
  bbox_out_big = bbox_out |> sf::st_buffer(dist=edge_buffer)

  # set cell values to pixel key
  cat('\nsetting up output grid')
  r_grid_in[] = terra::ncell(r_grid_in) |> seq()

  # define output grid: down-scale (no data) version of input via snapKrig then crop to bbox
  r_grid_out = r_grid_in |>
    snapKrig::sk() |>
    snapKrig::sk_rescale(down=down) |>
    snapKrig::sk_export() |>
    terra::crop(bbox_out_big, snap='out')

  # reshape as sk object and find mapping from un-cropped SpatRaster
  g_input = r_grid_out |> snapKrig::sk()
  is_obs = !is.na(g_input)
  idx_obs = g_input[is_obs]

  # garbage collection
  rm(r_grid_in)
  gc()

  # loop over variables
  for( v in seq_along(input_nc) ) {

    paste0('\n\nprocessing ', var_nm_list[[v]], '...') |> cat()
    from_v = from
    to_v = to

    # use the latest parameter fit
    pars_v = pars_all[[v]][[1]]

    # bilinear averaging to get DEM points on same grid as r_grid_out
    X_out = r_grid_out |> space_X(dem,
                                  dem_knots = pars_v[['knots']],
                                  X_center = pars_v[['center']],
                                  X_scale = pars_v[['scale']],
                                  intercept = FALSE)

    # garbage collection
    rm(r_grid_out)
    gc()

    # copy existing dates
    time_out = output_var_info[[v]][['time_obs']] |> as.Date()
    time_in = input_var_info[[v]][['time_obs']] |> as.Date()
    if( is.null(time_in) )  {

      cat('\nup to date')
      next
    }

    # set default starting/ending times
    is_initial = length(time_out) == 0
    if( is.null(to_v) ) to_v = max(time_in)
    if( is.null(from_v) ) {

      # on first call this writes everything
      from_v = min(time_in)

      # subsequently the default start time is 10 days before latest time
      if( !is_initial ) from_v = as.POSIXct(max(time_out)) - ( 10 * (60*60*24) )
    }

    # silently fix invalid start/end times
    from_v = as.Date(from_v)
    to_v = as.Date(to_v)
    if( from_v < min(time_in) ) from_v = min(time_in)
    if( to_v > max(time_in) ) to_v = max(time_in)

    # filter to requested range
    time_in = time_in[ ( time_in >= from_v ) & ( time_in <= to_v ) ]
    if( length(time_in) == 0 ) {

      cat('nothing to write\n')
      next
    }

    # split by year
    yr_v = time_in |> format('%Y', tz=tz)
    yr_unique = yr_v |> unique()

    # loop over years
    for( y in seq_along(yr_unique) ) {

      yr_y = yr_unique[y]
      cat('\nyear', yr_y)
      time_y = time_in[yr_v == yr_y]

      # load all required grid values to list of vectors (each one corresponds to a date)
      z_obs = input_nc[[v]] |> nc_layers(times=time_y) |> lapply(\(x) x[][idx_obs])

      # downscale in a loop over dates (returns in a list one SpatRaster per polygon)
      pars = pars_v[['pars']]
      r_output = wxArchive::.nc_downscale(g = g_input,
                                          z_obs = z_obs,
                                          pars = pars,
                                          X = X_out,
                                          poly_out = poly_out)

      # assign dates to layers then collect any garbage from .nc_downscale env
      terra::time(r_output) = time_y
      gc()

      # write outputs via `nc_write` directly to save memory
      cat('\nupdating .nc file(s)')
      p = nc_write_chunk(r=r_output, p=output_nc[[v]], path_only=TRUE)
      nc_write(r=r_output, p=p, insert=TRUE)
      cat('done\n')

      # remove the large source raster from memory
      rm(r_output)
      gc()
    }
  }
}


#' Down-scale a spatio-temporal time series
#'
#' This uses `snapKrig` and a previously fitted model (parameters list `pars`) to interpolate
#' from a grid of observed data onto a grid of finer resolution.
#'
#' The target grid and the location of the observations within it are specified by `g`, in
#' the same way as they are specified in `snapKrig::sk_cmean`; NAs for unobserved grid points
#' and non-NA values otherwise. The values themselves are ignored, and the data from `z_obs`
#' are substituted in their place.
#'
#' `z_obs` should be a list of vectors, one for each time layer to process. These vectors
#' should contain only the observed (non-NA) data, in the vectorized ordering used by `snapKrig`.
#' For example to downscale the data layer in `g`, you could set `z_obs = list(g[!is.na(g)])`.
#'
#' The data matrix `X` should be ordered the same way, but it must contain a row for every point
#' in the target grid, including unobserved ones. Typical usage is to create `X` with a call to
#' `space_X` (with `intercept=FALSE`).
#'
#' Output layers are returned in a SpatRaster. If `poly_out` is supplied, the function crops
#' all layers to the bounding box of the polygons as a whole.
#'
#' @param g sk grid containing the locations of observed grid points (see `?snapKrig::sk_cmean`)
#' @param z_obs list of vectors, the observed data at different times
#' @param pars list of spatial model parameters (eg as returned by `snapKrig::sk_fit`)
#' @param X numeric matrix of covariates passed (see `?snapKrig::sk_cmean`)
#' @param poly_out sf POLYGON objects to use for cropping output (or NULL for no cropping)
#'
#' @return either a SpatVector, a list of them
#' @export
.nc_downscale = function(g, z_obs, pars, X, poly_out=NULL) {

  is_obs = !is.na(g)
  n_layer = length(z_obs)
  is_cropped = !is.null(poly_out)

  # loop over layers
  cat('\ndownscaling', n_layer, 'layer(s)...\n')
  if(n_layer > 1) pb = utils::txtProgressBar(max=n_layer, style=3)
  list_out = vector(mode='list', length=n_layer)
  for( j in seq(n_layer) ) {

    # copy observations to snapKrig grid
    if(n_layer > 1) utils::setTxtProgressBar(pb, j)
    g[['gval']][is_obs] = z_obs[[j]]

    # universal kriging to fill the grid with predictions
    g_out = g |> snapKrig::sk_cmean(pars, X, quiet=TRUE)

    # convert to SpatRaster
    r_out = g_out |> snapKrig::sk_export()

    # optionally crop output to polygon
    if( is_cropped ) r_out = terra::crop(r_out, poly_out, snap='out')

    # copy result to storage
    list_out[[j]] = r_out
  }

  # merge layers into a single SpatRaster
  if(n_layer > 1) close(pb)
  r_out = list_out |> terra::rast()
  return(r_out)
}


