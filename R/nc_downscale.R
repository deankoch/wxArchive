# TODO: wrapper for down-scaling
#
# Since downscaling is slow, this will by default only update existing files,
# and write its outputs as nc split by year, similar to nc_update
#
# optionally also write plaintext output of weather aggregate values (compatible
# with rswat? or with SWAT+ directly?)
#
#

nc_downscale = function(base_dir,
                        dem,
                        input_nm = .nm_daily,
                        model_nm = .nm_model,
                        poly_out = NULL,
                        output_nm = .nm_export,
                        var_nm = .var_daily,
                        down = 100,
                        edge_buffer = NULL,
                        fun = NULL,
                        dates = NULL) {

  # var_nm is list to ensure input/output paths are lists too
  input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc))
  var_nm_list = names(var_nm) |> as.list()
  output_nc = file_wx('nc', base_dir, output_nm, var_nm_list, make_dir=TRUE)

  # read fitted parameter values of existing models for the variable(s)
  pars_json = file_wx('spatial', base_dir, model_nm, var_nm_list, make_dir=TRUE)
  pars_json = pars_json[ sapply(pars_json, file.exists) ]
  if( length(pars_json) == 0 ) stop('model file(s) not found in ', model_nm)
  pars_all =  pars_json |> lapply(jsonlite::fromJSON)

  # check for available input and existing output dates
  input_var_info = lapply(input_nc, \(p) time_wx(p))
  output_var_info = lapply(output_nc, \(p) time_wx(p))

  # get input grid info from first nc file
  r_grid_in = nc_chunk(input_nc[[1]][1])[1] |> terra::rast(lyrs=1)
  crs_in = terra::crs(r_grid_in)
  aoi_in = r_grid_in |> terra::ext() |> sf::st_bbox() |> sf::st_as_sfc()
  sf::st_crs(aoi_in) = crs_in

  # set default polygon (whole extent)
  if( is.null(poly_out) ) poly_out = aoi_in
  poly_out = poly_out |> sf::st_geometry() |> sf::st_transform(crs_in)
  poly_list = poly_out |> split( seq_along(poly_out) )

  # set default buffer size based on source resolution
  if( is.null(edge_buffer) ) edge_buffer = sum(terra::res(r_grid_in)^2) |> sqrt()

  # transform output polygons to input projection and find their bounding box
  bbox_out = sf::st_bbox(poly_out) |> sf::st_as_sfc()
  bbox_out_big = bbox_out |> sf::st_buffer(dist=edge_buffer)

  # DEBUGGING
  plot(aoi_in)
  plot(bbox_out_big, add=T)
  plot(bbox_out, add=T)
  plot(poly_out, add=T)

  # set cell values to pixel key
  r_grid_in[] = terra::ncell(r_grid_in) |> seq()

  # get down-scaled (no data) version of input via snapKrig then crop to bbox
  r_grid_out = r_grid_in |>
    snapKrig::sk() |>
    snapKrig::sk_rescale(down=down) |>
    snapKrig::sk_export() |>
    terra::crop(bbox_out_big, snap='out')

  # snapKrig version and mapping from uncropped SpatRaster
  g_input = r_grid_out |> snapKrig::sk()
  is_obs = !is.na(g_input)
  idx_obs = g_input[is_obs]

  # select a variable
  v = 4

  # copy existing dates
  time_out = output_var_info[[v]][['time_obs']]
  time_in = input_var_info[[v]][['time_obs']]
  if( is.null(time_in) ) next

  # set default dates to modify
  time_mod = time_in
  if( !is.null(time_out) ) {

    # existing dates are left alone by default
    if( is.null(dates) ) dates = time_mod[ !( time_mod %in% time_in ) ]
  }

  # filter to the user-specified range
  if( !is.null(dates) ) time_mod = time_mod[ time_mod %in% dates ]
  if( length(time_mod) ) next

  # use the latest parameter fit
  pars_v = pars_all[[v]][[1]]

  # bilinear averaging to get DEM points on same grid as r_grid_out
  X_out = r_grid_out |> space_X(dem,
                                dem_knots = pars_v[['knots']],
                                X_center = pars_v[['center']],
                                X_scale = pars_v[['scale']],
                                intercept = FALSE)

  # loop over years
  time_split = time_mod |> split(format(time_mod, '%Y'))
  for( i in seq_along(time_split) ) {

    # load all requested layers and copy the required grid values to list of vectors
    z_obs = input_nc[[v]] |> nc_layers(times=time_split[[i]]) |> lapply(\(x) x[][idx_obs])

    # downscale layers in a loop
    g_output = .nc_downscale(g_input, z_obs, pars_v[['pars']], X_out, poly_list)
    if( !is.list(g_output) ) g_output = list(g_output)
    for(j in length(g_output) ) terra::time( g_output[[j]] ) = time_split[[i]]
  }






  g_output[[12]][[1]] |> plot()
  plot(poly_out, add=T)



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
#' Output layers are returned in a SpatRaster. If `poly_list` is supplied, the function returns
#' a list of SpatRasters cropped to each of the polygons' extents.
#'
#' @param g sk grid containing the locations of observed grid points (see `?snapKrig::sk_cmean`)
#' @param z_obs list of vectors, the observed data at different times
#' @param pars list of spatial model parameters (eg as returned by `snapKrig::sk_fit`)
#' @param X numeric matrix of covariates passed (see `?snapKrig::sk_cmean`)
#' @param poly_list list of sf POLYGON objects to use for cropping output (or NULL for no cropping)
#'
#' @return either a SpatVector, a list of them
#' @export
.nc_downscale = function(g, z_obs, pars, X, poly_list=NULL) {

  is_obs = !is.na(g)
  n_layer = length(z_obs)
  is_by_poly = !is.null(poly_list)
  is_aggregate = !is.null(fun)

  # loop over layers
  cat('\nprocessing', n_layer, 'layer(s)...\n')
  if(n_layer > 1) pb = utils::txtProgressBar(max=n_layer, style=3)
  list_out = vector(mode='list', length=n_layer)
  for( j in seq(n_layer) ) {

    # copy observations to snapKrig grid
    if(n_layer > 1) utils::setTxtProgressBar(pb, j)
    g[['gval']][is_obs] = z_obs[[j]]

    # universal kriging
    g_out = g |> snapKrig::sk_cmean(pars, X, quiet=TRUE)
    r_out = g_out |> snapKrig::sk_export()

    # optionally split output by polygon
    if( is_by_poly ) {

      # set default polygon list
      if( !is_by_poly ) poly_list = terra::ext(r_out) |> sf::st_bbox() |> sf::st_as_sfc() |> list()

      # crop to individual polygons
      r_out = poly_list |> lapply(\(p) terra::crop(r_out, p, snap='out'))
    }

    # copy result to storage
    list_out[[j]] = r_out
  }

  # aggregate values returned in a wide matrix
  if(n_layer > 1) close(pb)

  # reshape output
  if(is_by_poly) {

    # else make a time series for each polygon
    n_poly = length(poly_list)
    r_out = seq(n_poly) |> lapply(\(i) {

      seq(n_layer) |> lapply(\(j) list_out[[j]][[i]] ) |> terra::rast()

    })

  } else {

    # bind all layers into a single SpatRaster time series
    r_out = list_out |> terra::rast()
  }

  return(r_out)
}




#' Aggregate raster data over a set of polygons
#'
#' This aggregates the values in SpatRaster `r` by masking to the polygon(s) in `poly_list`
#' and applying the function `fun` to each subset.
#'
#' `r` can have multiple layers. A single numeric is returned for each polygon, forming
#' a vector of aggregate values at different times. When `poly_list` has length 1 this is
#' returned as a vector, otherwise it is returned as a matrix with one row for each polygon.
#'
#' `fun` specifies the name of the aggregation function to use - typical usage will set one
#' of 'mean' 'max' and 'min' (see `?terra::extract`). Alternatively, set `fun=NULL` to return
#' a list of SpatRasters cropped to `poly_list` (one per polygon). When `fun` is not `NULL`,
#' this list is passed (along with `fun`) to `terra::extract` and the result is reshaped as
#' described above.
#'
#' @param r SpatRaster with extent covering `poly_list`
#' @param poly_list list of sf POLYGON objects to use for cropping output (or NULL for no cropping)
#' @param fun character name of the function to use for aggregation (or NULL for no aggregation)
#'
#' @return either a SpatVector, a list of them
#' @export
.nc_spatial_aggregate = function(r, fun='mean', poly_list=NULL) {

  # set default polygon list
  if( is.null(poly_list) ) poly_list = terra::ext(r) |> sf::st_bbox() |> sf::st_as_sfc() |> list()

  # crop to individual polygons
  r = poly_list |> lapply(\(p) terra::crop(r, p, snap='out'))
  if( is.null(fun) ) return(r)

  # loop over polygons, results in a list of data frames
  r_stat_list = Map(\(x, p) terra::extract(x = x,
                                           y = as(p, 'SpatVector'),
                                           fun = fun,
                                           raw = TRUE,
                                           touches = TRUE,
                                           ID = FALSE), x=r, p=poly_list)

  # reshape as vector
  v_out = do.call(c, lapply(r_stat_list, unlist)) |> unname()
  return(v_out)
}

