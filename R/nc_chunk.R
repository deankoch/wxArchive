#' Expand a directory path to get a set of paths to the NetCDF files within
#'
#' It can be convenient to split a long time series into year-long chunks
#' and write each chunk to its own NetCDF (.nc) file. For variable "var", these
#' files are written to a directory named "var.nc". This function looks for files
#' in this directory and returns their paths.
#'
#' With the default `year=NULL` the function scans for files at `p` and returns what
#' it finds there: either `p` itself if the path points to a NetCDF file, or a vector
#' of paths to all NetCDF files inside the directory `p`. If the directory exists but is
#' empty, the function returns `character(0)`.
#'
#' If `year` is supplied, the function returns file paths for each of `unique(year)`,
#' without checking if the parent directory exists.
#'
#' The expected file name for a chunk is "var_year.nc", with "year" the 4-digit year
#' of all times in the file. For example "tmp.nc/tmp_2005.nc" holds observations of
#' temperature from 2005.
#'
#' @param p character path to a file or directory with extension '.nc'
#' @param year (optional) character or integer vector of (4-digit) years
#'
#' @return character vector of file path(s) to NetCDF file(s) associated with `p`
#' @export
nc_chunk = function(p, year=NULL) {

  # validity check for arguments
  if( !is.character(p) ) stop('is.character(p) was FALSE. p must be a string')
  if( endsWith(p, '.nc') ) stop('path ', p, ' should have extension ".nc"')
  if( is.null(year) ) { if( !file.exists(p) ) stop('path ', p, ' not found on disk') } else {

    # handle simple calls with known years
    year = as.integer(year)
    f = gsub('.nc', paste0('_', as.character(year)), basename(p)) |> paste0('.nc')
    return( file.path(p, f) )
  }

  # if the file is a directory scan its contents
  if( dir.exists(p) ) {

    p_contents = list.files(p)
    f = p_contents[ endsWith(p_contents, '.nc') ]
    return( file.path(p, f) )
  }

  # p points to a single NetCDF file
  return(p)
}
