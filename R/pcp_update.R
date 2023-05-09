#' Create/update a netCDF file with missing pcp_total times filled using pcp_large and pcp_small
#'
#' This creates a copy of the "pcp_total" time series, named `pcp_nm`, containing only
#' times imputed using the sum of "pcp_large" (large scale precipitation) and "pcp_small"
#' (convective precipitation) whenever both are available, but "pcp_total" is not.
#'
#' `input_nm` should be a named list character vector elements 'coarse' and/or 'fine', each
#' representing a batch of files to read a common variable from each resolution. Typical usage
#' has at least two source files at each resolution - a long term storage file, and a smaller file
#' with more recent times.
#'
#' The function writes one output file for each of the vectors in `input_nm`, to the
#' location(s) in `output_nm`. The output file name is `paste0(pcp_nm, 'nc')`, so the
#' function will not allow you to set `pcp_nm` to any of the (in-use) names to avoid
#' confusion downstream (list them with `file_wx('nc', base_dir, output_nm)`)
#'
#' Layers with negatives are set to `NA` in the output file.
#'
#' A JSON for the output is written to the "time" subdirectory. The function will update
#' existing files by only copying layers for times not already listed in the file.
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param pcp_nm a name for the output variable and file
#' @param input_nm list of character vectors, naming input subdirectories in `base_dir`
#' @param output_nm character vectors, naming an output subdirectory for each `input_nm` vector
#'
#' @return nothing, but possibly modifies the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
pcp_update = function(base_dir,
                      pcp_nm = 'pcp',
                      input_nm = list(coarse='coarse', fine='fine'),
                      output_nm = c('coarse', 'fine')) {

  # paths to expected inputs
  pcp_var_nm = list('pcp_total', 'pcp_large', 'pcp_small')
  input_path = file_wx('nc', base_dir, input_nm, var_nm=pcp_var_nm)

  # sanity check for pcp_nm
  nm_var_all = stats::setNames(nm=names(input_path[[1]]))
  if( pcp_nm %in% nm_var_all ) stop('pcp_nm cannot be an existing variable name')

  # new files to write
  output_nc_path = file_wx('nc', base_dir, as.list(output_nm), var_nm=pcp_nm)

  # find time coverage of each variable at both resolutions (creating JSON as needed)
  cat('checking available times for', paste(nm_var_all, collapse=', '))
  var_info = lapply(input_path, \(r) lapply(r, \(p) time_wx(p)) )

  # loop over resolutions
  for( nm_res in output_nm ) {

    paste0('\n\n', nm_res, ' grids...') |> cat()
    t1 = proc.time()

    # check for available times by variable
    t_all = do.call(c, lapply(var_info[[nm_res]], \(v) v[['time']])) |> unique() |> sort()
    is_obs = var_info[[nm_res]] |> lapply(\(v) t_all %in% v[['time_obs']])
    is_relevant = is_obs[['pcp_large']] & is_obs[['pcp_small']] & !is_obs[['pcp_total']]
    time_relevant = t_all[is_relevant]

    # check for existing data (creating JSON as needed)
    p_out = output_nc_path[[nm_res]]
    time_done = character(0) |> as.POSIXct()
    if( file.exists(p_out) ) time_done = time_wx(p_out)[['time']]

    # select required layers from each source
    is_needed = !( time_relevant %in% time_done )
    time_add = time_relevant[is_needed]

    # if there's nothing to add then we are finished
    if( length(time_add) == 0 ) {

      cat('\nup to date \U2713')
      next
    }

    # compute totals from large and small source
    r_add = NULL
    if( length(time_add) > 0 ) {

      cat('\ncopying sum of pcp_large and pcp_small')
      r_large = input_path[[nm_res]][['pcp_large']] |> nc_layers(time_add)
      r_small = input_path[[nm_res]][['pcp_small']] |> nc_layers(time_add)
      r_add = r_large + r_small
    }

    # handle NAs encoded as negatives
    cat('\nchecking for negatives... ')
    is_na = unlist(terra::global(r_add, 'min')) < 0
    n_na = sum(is_na)
    if( n_na > 0 ) {

      r_add[seq(nrow(r_add)), seq(ncol(r_add)), which(is_na)] = NA
      cat('omitting', n_na, 'NA layers')

    } else { cat('none found') }

    # append to existing data file (or create the file and write to it)
    r_add |> nc_write(p_out)
    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.')
  }
  cat('\n')
}
