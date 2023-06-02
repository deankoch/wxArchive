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
