#' Create/update a netCDF file with wind speed, computed from its components
#'
#' This creates a new variable `wnd_nm` reporting the norm (speed) of the
#' wind vector, computed from the two components (u and v). These components should
#' be found in the length-2 names vector `uv_nm`. Specify the directory names in
#' which to look for these variables with `input_nm` (see also `?pcp_update`).
#'
#' `output_nm` may contain more than one directory, so that the output can be split
#' into a long-term storage file and a recent updates file. The first directory named
#' in `output_nm` is always used for writing changes.
#'
#' The function creates/modifies a NetCDF file in the directory named `output_nm[1]`.
#' This file is named `paste0(wnd_nm, '.nc')` similar the other variables, so the
#' function will not allow you to set `wnd_nm` to any of the (in-use) variable
#' names (list them with `file_wx('nc', base_dir, output_nm)`).
#'
#' An output time (layer) is generated for each time in which both wind component
#' variables are observed in at least one the files in `input_nm`. The function only
#' appends those times not already found in the output file when overwriting it.
#' A JSON file in the "time" subdirectory is also created/updated (see `?nc_write`).
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param wnd_nm a name for the output variable and file
#' @param input_nm list of character vectors, naming input subdirectories in `base_dir`
#' @param output_nm character vector, naming the output subdirectories
#'
#' @return nothing, but possibly modifies the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
wnd_update = function(base_dir,
                      wnd_nm = .var_wnd,
                      uv_nm = .var_wnd_uv,
                      input_nm = .nm_complete_rap,
                      output_nm = wnd_nm) {

  # paths to expected inputs
  input_path = file_wx('nc', base_dir, input_nm, var_nm=as.list(uv_nm))

  # sanity check for wnd_nm and uv_nm
  if( length(uv_nm) != 2 ) stop('uv_nm must name exactly two variables')
  nm_var_all = stats::setNames(nm=names(input_path))
  if( wnd_nm %in% nm_var_all ) stop('wnd_nm cannot be ',  uv_nm[1], ' or ', uv_nm[2])

  # new file to write
  output_nc_path = file_wx('nc', base_dir, output_nm, var_nm=wnd_nm, make_dir=TRUE)

  # find common time coverage (creating JSON as needed)
  cat('checking available times for', paste(nm_var_all, collapse=', '))
  var_info = lapply(input_path, \(p) time_wx(p) )
  if( any( sapply(var_info, length) == 0 ) ) stop('one or both variables had no observed data')
  t_relevant = Reduce(intersect, lapply(var_info, \(x) x[['time_obs']])) |> as.POSIXct(tz='UTC')

  # check for existing data (creating JSON as needed)
  cat('\nchecking for existing file', output_nc_path)
  t_done = character(0) |> as.POSIXct()
  if( file.exists(output_nc_path) ) t_done = time_wx(output_nc_path)[['time']]

  # new times to add to output file
  t_add = t_relevant[ !( t_relevant %in% t_done ) ]

  # if there's nothing to add then we are finished
  if( length(t_add) == 0 ) {

    cat('\nup to date')
    return( invisible() )
  }

  # load into RAM as matrix and compute wind speed (m/s) as norm of vector
  cat('\nloading', length(t_add), 'layers from', uv_nm[1], 'and', uv_nm[2])
  r_u = nc_layers(input_path[[ uv_nm[1] ]], t_add, na_rm=TRUE)
  r_v = nc_layers(input_path[[ uv_nm[2] ]], t_add, na_rm=TRUE)

  cat('\ncomputing wind speed')
  r = sqrt(r_u^2 + r_v^2)

  # append to existing data file (or create the file and write to it)
  r |> nc_write(output_nc_path)
  cat('\n')
  return( invisible() )
}
