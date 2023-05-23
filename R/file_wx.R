#' Build file paths for writing outputs
#'
#' This organizes the input/output files for various steps in the RAP workflow. When
#' `make_dir=TRUE`, the function creates the parent directory containing the requested
#' file name (if it doesn't already exist).
#'
#' The argument `sub_dir` defines the directory structure. eg by default we have:
#'
#' * .../grib/ stores the many source files (.grb and .grb2)
#' * .../coarse/ stores NetCDF versions (.nc) for each variable at coarse resolution
#' * .../fine/ stores NetCDF versions (.nc) for each variable at fine resolution
#'
#' Note that a "file" may refer to a directory with the suffix ".nc" in its name. Such
#' directories should contain a set of yearly NetCDF files named as defined by `nc_chunk`.
#'
#' Argument `what` specifies the path(s) to return:
#'
#' 1. 'grib' : character, the source GRIB file directory
#' 2. 'nc' : a list of character vectors, paths to the .nc file
#' 3. 'index' : a list of character vectors, paths to JSON files storing time index
#' 4. 'spatial' : a list of character vectors, paths to JSON files storing fitted spatial parameters
#' 5. 'temporal_index' : character, path to a JSON file storing information about temporal fit
#' 6. 'temporal_nc' : character, path to a new .nc file for storing fitted temporal parameters
#'
#' With options (1) and (5-6) the function returns a character path, whereas with (2-4) it
#' returns a list of character path vectors - one list entry per resolution, and one vector
#' entry per variable name. If `sub_dir` has length 1 (one resolution only) then the result is
#' automatically unlisted and the function returns a character vector.
#'
#' With option (6), the file name will contain the current system time - use this only when
#' creating a new file, then append that file name to the JSON from (5) for later use.
#'
#' When `collapse=TRUE`, the top level of the list result is collapsed by concatenating
#' its contents using the names from the first list entry. This has no effect when
#' the output would otherwise have length 1, or when `sub_dir` has length 1.
#'
#' Both `sub_dir` and `var_nm` can be vectors or a lists of vectors. Lists are used to group
#' like subdirectories (eg different pieces of the same time series) or variable names that
#' should be viewed as equivalent. The function loops over lists, passing each element to itself
#' in a recursive call and collapsing the results. This makes a (possibly nested) list of
#' vectors, the result of looping over `sub_dir` first (outer), then `var_nm` (inner).
#'
#' @param base_dir character path the parent directory for RAP file storage
#' @param sub_dir character vector or list of them, names of sub-directories in `base_dir`
#' @param what one of 'grib', 'csv', 'nc', 'index', 'spatial', 'temporal' (see details)
#' @param var_nm character vector of variable names, or a list of them
#' @param make_dir logical indicating to create directories where needed
#'
#' @return a named character vector, or a named list of them (depending on `what` and `sub_dir`)
#' @export
file_wx = function(what,
                   base_dir = '',
                   sub_dir = list('coarse', 'fine'),
                   var_nm = as.list(names(.rap_regex)),
                   make_dir = FALSE) {

  # set default variable names (first element of each vector in list)
  var_nm = setNames(var_nm, nm = sapply(var_nm, \(x) x[1]))
  sub_dir = setNames(sub_dir, nm = sapply(sub_dir, \(x) x[1]))

  # build lists of paths using filename(s) (p) and one or more directories (d)
  my_dir = function(d, f=NULL) {

    # create the directory, return it if no file names supplied
    for(x in d) if( !file.exists(x) & make_dir ) dir.create(x, recursive=TRUE)
    if( is.null(f) ) return(d)
    out_list = lapply(d, \(x) file.path(x, f)) |> unlist()
    return(out_list)
  }

  # grib and metadata storage locations
  if( what == 'grib' ) return( my_dir(file.path(base_dir, 'grib')) )
  if( what == 'csv' ) return( my_dir(file.path(base_dir, 'grib'), 'grib_df.csv') )

  # recursion for list input to sub_dir
  if( is.list(sub_dir) ) {

    # loop over sub-directory lists and collapse result
    out_list = sub_dir |> lapply(\(s) file_wx(what = what,
                                              base_dir = base_dir,
                                              sub_dir = s,
                                              var_nm = var_nm,
                                              make_dir = make_dir))


  } else {

    # loop over variable name lists and collapse result
    if( is.list(var_nm) ) {

      out_list = var_nm |> lapply(\(v) file_wx(what = what,
                                               base_dir = base_dir,
                                               sub_dir = sub_dir,
                                               var_nm = v,
                                               make_dir = make_dir))


    } else {

      # all recursive function calls should end up here

      # get a string representing current date and time
      call_time = gsub('[\\: ]', '.', as.POSIXct(Sys.time(), tz='UTC'), perl=TRUE)

      # sub-directories
      sub_path = file.path(base_dir, sub_dir)
      dir_nm = list(grib = file.path(base_dir, 'grib'),
                    nc = stats::setNames(sub_path, sub_dir),
                    index = stats::setNames(file.path(sub_path, 'time'), sub_dir),
                    spatial = stats::setNames(sub_path, sub_dir),
                    temporal_index = stats::setNames(sub_path, sub_dir),
                    temporal_nc = stats::setNames(file.path(sub_path, 'temporal'), sub_dir))

      # file names
      file_nm = list(nc = stats::setNames(paste0(var_nm, '.nc'), var_nm),
                     index = stats::setNames(paste0(var_nm, '.json'), var_nm),
                     spatial = stats::setNames(paste0(var_nm, '_spatial.json'), var_nm),
                     temporal_index = paste0(var_nm, '_temporal.json'),
                     temporal_nc = paste0(var_nm, '_', call_time, '.nc'),
                     csv = 'grib_df.csv')

      # nc file locations
      out_list = my_dir(dir_nm[[what]], file_nm[[what]])
      len_match = length(var_nm) == length(out_list)
      if( is.null(names(out_list)) & len_match ) names(out_list) = var_nm
    }
  }

  return(out_list)
}
