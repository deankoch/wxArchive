#' fun character, naming a function to use for spatial aggregation (NA to disable)
#' #'
#' writes a separate NetCDF file for each one, using an integer-valued suffix
#' in the filename corresponding to the order of the polygons in `poly_out`.
#'
#' If `fun` is supplied, the function uses the named function to aggregate the
#' data values over the supplied polygons (within each date), and writes its
#' results in a plain text table on disk. Set `write_nc=FALSE` to skip writing
#' the NetCDF files (the function will still write the aggregate data). Set
#' `fun=NA` to skip aggregation.
#'
#' With default `fun=NULL` the function expects variable names of the form
#' "var_fun"; Names are split at the "_", and the second string is copied to `fun`.
#' eg if `fun=NULL` and `var_nm = 'tmp_max'` then the function will compute the
#' `max` on each polygon, on each date.




# # set default polygon (whole extent)
# if( is.null(poly_out) ) poly_out = aoi_in
# poly_out = poly_out |> sf::st_geometry() |> sf::st_transform(crs_in)
# poly_list = poly_out |> split( seq_along(poly_out) )

# # DEBUGGING
# plot(aoi_in)
# plot(bbox_out_big, add=T)
# plot(bbox_out, add=T)
# plot(poly_out, add=T)
#fun = 'mean'




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
