#' Return a data matrix of spatial linear predictors using a DEM
#'
#' This creates predictors for Northing and Easting (or latitude and longitude),
#' as well as a spline basis for the lapse rate effect of elevation. All columns
#' are centered and scaled by standard deviation, in that order.
#'
#' The output matrix has a row for each (spatial) point in `r` and a column for
#' every predictor (`1 + length(dem_knots)`, with default `intercept=FALSE`).
#'
#' Rows are in the column major ordering used by `snapKrig` (but not `terra`).
#'
#' The centering/scaling constants and the spline knots are provided in attributes
#' of the returned matrix. These can be copied and passed back to this function in
#' subsequent calls with a different `dem` or `r`. This will ensure that your fitted
#' `betas` are scaled properly when predicting at unseen locations.
#'
#' It is assumed that `dem` covers the extent of `r`. `dem` can have a different
#' projection than `r`, but for best results it should have a much finer resolution.
#'
#' The spline basis is created using `splines::ns`. If `dem_knots` is supplied,
#' it should include boundary points. If it is not supplied, the function puts
#' knots at the default `quantile` values of `dem`, after cropping `dem` to the
#' extent of `r` (but before resampling).
#'
#' @param r SpatRaster providing the target grid (data layers are ignored)
#' @param dem SpatRaster of elevation data at finer resolution than `r`
#' @param dem_knots numeric vector of elevation knot locations (in metres)
#' @param X_center numeric vector of centering constants
#' @param X_scale numeric vector of scaling constants
#' @param intercept logical indicating to include an intercept column (of `1`s)
#'
#' @return a matrix of covariates
#' @export
space_X = function(r, dem, dem_knots=NULL, X_center=NULL, X_scale=NULL, intercept=FALSE) {

  # use bilinear averaging to project onto forecast grid
  dem_warp = dem |> terra::project(r) |> snapKrig::sk()

  # get dem_knots (if not supplied) from quantiles of dem
  if( is.null(dem_knots) ) {

    bbox_crop = sf::st_bbox(r) |> sf::st_as_sfc() |> sf::st_transform(terra::crs(dem))
    dem_crop = dem |>  terra::crop(as(bbox_crop, 'Spatial')) |> snapKrig::sk()
    dem_knots = dem_crop[] |> stats::quantile()
  }

  # make a spline basis for the DEM
  dem_knots_inner = dem_knots |> head(-1) |> tail(-1)
  dem_basis = dem_warp[] |> splines::ns(knots=dem_knots_inner, Boundary.knots=range(dem_knots))
  colnames(dem_basis) = paste0('dem_spline_', colnames(dem_basis))

  # other simple covariates
  northing = snapKrig::sk(dem_warp, gval=snapKrig::sk_coords(dem_warp, out='list', quiet=T)[['y']])
  easting = snapKrig::sk(dem_warp, gval=snapKrig::sk_coords(dem_warp, out='list', quiet=T)[['x']])

  # combine into spatial covariates matrix, scale columns
  X_space_unscaled = cbind(y = northing[], x = easting[], dem_basis)
  if( is.null(X_center) ) X_center = X_space_unscaled |> apply(2, mean, na.rm=TRUE)
  X_space = X_space_unscaled |> sweep(2, X_center)
  if( is.null(X_scale) ) X_scale = X_space |> apply(2, sd, na.rm=TRUE)
  X_space = X_space |> sweep(2, X_scale, '/')
  colnames(X_space) = colnames(X_space_unscaled)

  # add intercept column if requested
  if( intercept ) X_space = cbind(rep(1, nrow(X_space)), X_space)

  # copy attributes from scaling and splines to output
  attr(X_space, 'Boundary.knots') = attr(dem_basis, 'Boundary.knots')
  attr(X_space, 'knots') = attr(dem_basis, 'knots')
  attr(X_space, 'center') = X_center
  attr(X_space, 'scale') = X_scale
  return(X_space)
}
