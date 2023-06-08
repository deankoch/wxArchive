#' Expand a directory path to get a set of paths to the NetCDF files within
#'
#' It can be convenient to split a long time series into year-long chunks
#' and write each chunk to its own NetCDF (.nc) file. For variable "var", these
#' files are written to a directory named "var.nc". This function looks for files
#' in this directory and returns their paths.
#'
#' The function scans for files at `p` and returns what it finds there: either `p`
#' itself if the path points to a NetCDF file, or a vector of paths to all NetCDF
#' files inside the directory `p`. If the directory exists but is empty, the function
#' returns `character(0)`. If the directory doesn't exist the function returns `NA`.
#'
#' The expected file name for a chunk is "var_year.nc", with "year" the 4-digit year
#' of all times in the file. For example "tmp.nc/tmp_2005.nc" holds observations of
#' temperature from 2005.
#'
#' If `file_ext='tif'`, all of the above applies, but the function looks for output
#' files of type GeoTIFF (rather than NetCDF), having extension '.tif' (rather than
#' '.nc'). Note that the directory name may still end with '.nc' in this case, but
#' only filenames ending in 'tif' are returned.
#'
#' @param p character path to a file or directory with extension ".nc"
#' @param file_ext character, either 'tif' or 'nc'
#'
#' @return character vector of file path(s) to NetCDF file(s) associated with `p`
#' @export
nc_chunk = function(p, file_ext='.nc') {

  # validity check for arguments
  if( !is.character(p) ) stop('is.character(p) was FALSE. p must be a string')
  if( !file.exists(p) ) return(NA)

  # if the file is a directory scan its contents
  if( dir.exists(p) ) {

    p_contents = list.files(p)
    f = p_contents[ endsWith(p_contents, file_ext) ]
    return( file.path(p, f) )
  }

  # else p points to a single existing file
  return(p)
}


#' Write time series data to NetCDF files that are chunked by year
#'
#' A wrapper for `nc_write` with multiple output files. This splits the data
#' from `r` by year and writes a separate file for each year in directory `p`
#' according to the naming scheme in `?nc_chunk`.
#'
#' `p` must be a directory and its name must end with ".nc"
#'
#' Set `name_only` to return the file names that would be written, but not
#' actually write anything.
#'
#' Note that the year for a given time can depend on the time zone. Years are
#' delineated by this function in UTC.
#'
#' @param r SpatRaster to write
#' @param p path to the output directory
#' @param path_only logical, if TRUE the function creates the directory but writes nothing to it
#' @param insert logical, enables replacement of existing times (passed to `nc_write`)
#'
#' @return a list of times (the result of `nc_write` for each year in `r`)
#' @export
nc_write_chunk = function(r, p, path_only=FALSE, insert=FALSE) {

  # create/scan output directory for chunked files
  p_chunk = nc_chunk(p)
  if( anyNA(p_chunk) ) dir.create(p, recursive=TRUE)
  p_chunk = nc_chunk(p)
  if( anyNA(p_chunk) ) stop('there was a problem creating the directory ', p)

  # get variable name from path
  nm = gsub('.nc', '', basename(p))

  # check input times and get output filenames
  all_years = terra::time(r) |> format('%Y', tz='UTC')
  if( is.null(all_years) ) stop('terra::time(r) was NULL (expected POSIXct times or Dates)')
  years = unique(all_years) |> sort()
  f = paste0(nm, '_', years, '.nc') |> stats::setNames(years)
  if( path_only ) return(file.path(p, f))

  # create/overwrite chunks in a loop over years
  cat('\n')
  for(yr in years) {

    # explicitly remove the large raster from RAM after the write completes
    nc_write(r=r[[all_years == yr]], p=file.path(p, f[yr]), insert=insert)
    gc()
  }

  return( invisible() )
}

