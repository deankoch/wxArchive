#' List the expected NetCDF files in the archive
#'
#' This returns a list of paths to the NetCDF data files created by this package
#' with default settings. If either of `rap` or `gfs` is NULL, files from that
#' model are omitted.
#'
#' Variable names are passed in a list, and can split into groups of files (not
#' necessarily having the same name) that are considered equivalent when loading
#' data. eg we consider 'pcp' and 'pcp_total' to be equivalent by including
#' `c('pcp', 'pcp_total')` in the first element of `var_rap` (see `?file_wx`).
#'
#' @param base_dir path to parent directory of GRIB storage subfolders
#' @param rap character vector of sub-directory names for RAC/RUC files
#' @param gfs character vector of sub-directory names for GFS files
#' @param var_rap list of character vectors, naming the RAP/RUC variables to include
#' @param var_gfs list of character vectors, naming the GFS variables to include
#'
#' @return character vector or list of paths (depends if `rap` and/or `gfs` are NULL)
#' @export
nc_list = function(base_dir,
                   rap = .nm_rap_export,
                   gfs = .nm_gfs_export,
                   var_rap = .var_rap_export,
                   var_gfs = .var_gfs_export) {

  # NULL specifies to omit from results
  is_rap = !is.null(rap)
  if( is_rap ) {

    # base directories for all RAP/RUC files
    base_dir_rap = base_dir |> file.path('rap')
    rap_nc_path = file_wx('nc', base_dir_rap, rap, var_rap)
  }

  # repeat for GFS paths
  is_gfs = !is.null(gfs)
  if( is_gfs ) {

    # base directories for all GFS files
    base_dir_gfs = base_dir |> file.path('gfs')
    gfs_nc_path = file_wx('nc', base_dir_gfs, gfs, var_gfs)
  }

  # return paths from GFS or RAP/RUC (not both)
  if( !(is_gfs|is_rap) ) return( character(0) )
  if( is_rap & !is_gfs ) return( rap_nc_path )
  if( is_gfs & !is_rap ) return( gfs_nc_path )

  # return both (order establishes preference for RAP/RUC archive over GFS)
  p_all = Map(\(rap, gfs) c(rap, gfs), rap = rap_nc_path, gfs = gfs_nc_path)
  return(p_all)
}
