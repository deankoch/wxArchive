#' Expand a directory path to get a set of paths to the NetCDF files within
#'
#' It can be convenient to split a long time series into year-long chunks
#' and write each chunk to its own NetCDF (.nc) file. For variable "var", these
#' files are written to a directory named "var.nc". Pass the full directory path to
#' this function to get a vector of paths to the yearly chunks.
#'
#' If `year` is supplied, the function returns corresponding file paths for each
#' of `unique(year)`. With default `year=NULL` the function scans for files in `p`
#' and returns what it finds (or `character(0)` if empty, or NA if not a directory).
#'
#' The expected file name for a chunk is "var_year.nc", with "year" the 4-digit year
#' of all times in the file (eg "tmp.nc/tmp_2005.nc" holds observations of temperature
#' from 2005).
#'
#' @param p character path to a file or directory with extension '.nc'
#' @param year (optional) character or integer vector of (4-digit) years
#'
#' @return
#' @export
nc_chunk = function(p, year=NULL) {


}
