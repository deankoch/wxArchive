nc_downscale = function() {



}


#' Down-scale a spatio-temporal time series and optionally aggregate output over polygons
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
#' in the target grid (including unobserved).
#'
#' Output layers are returned in a SpatRaster. If `poly_list` is supplied, the function returns
#' a list of SpatRasters, one per output polygon, cropped to the polygon extent. If `fun` is
#' supplied, the function aggregates the data from (each) SpatRaster over the corresponding
#' polygon(s) to return a single number for each polygon/time.
#'
#' @param g sk grid containing the locations of observed grid points (see `?snapKrig::sk_cmean`)
#' @param z_obs list of vectors, the observed data at different times
#' @param pars list of spatial model parameters (eg as returned by `snapKrig::sk_fit`)
#' @param X numeric matrix of covariates passed (see `?snapKrig::sk_cmean`)
#' @param poly_list list of sf POLYGON objects to use for cropping output (or NULL for no cropping)
#' @param fun character name of the function to use for aggregation (or NULL for no aggregation)
#'
#' @return either a SpatVector, a list of them, or a numeric matrix (when `fun` is specified)
#' @export
.nc_downscale = function(g, z_obs, pars, X, poly_list=NULL, fun=NULL) {

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
    if( is_by_poly | is_aggregate ) {

      # set default polygon list
      if( !is_by_poly ) poly_list = terra::ext(r_out) |> sf::st_bbox() |> sf::st_as_sfc() |> list()

      # crop to individual polygons
      r_out = poly_list |> lapply(\(p) terra::crop(r_out, p, snap='out'))

      # optionally aggregate within polygons
      if( is_aggregate ) {

        # results in a list
        r_stat_list = Map(\(r, p) terra::extract(x = r,
                                                 y = as(p, 'SpatVector'),
                                                 fun = fun,
                                                 raw = TRUE,
                                                 touches = TRUE,
                                                 ID = FALSE), r=r_out, p=poly_list)

        # replace with vector
        r_out = do.call(c, lapply(r_stat_list, unlist)) |> unname()
      }
    }

    # copy result to storage
    list_out[[j]] = r_out
  }

  # aggregate values returned in a wide matrix
  if(n_layer > 1) close(pb)
  if(is_aggregate) return( do.call(cbind, list_out) )

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
