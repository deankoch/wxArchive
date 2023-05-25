#' Create or update a NetCDF version of a local RUC/RAP/GFS archive
#'
#' Most GRIB files fetched by `archive_update` are around 20-30 MB compressed,
#' so loading tens of thousands at once is extremely slow. To make the time series
#' easier to work with, this function extracts the subset that we need and saves it in
#' a format that can be loaded more quickly.
#'
#' The data are merged into a single (multi-layer) NetCDF file for each of the variables
#' selected by `regex`. Only the sub-grid overlapping with `aoi` is copied.
#'
#' Output files are given base names `names(regex)`, and extension '.nc', and
#' have layers named after the POSIXct time (as character, in GMT) at which the
#' forecast is valid. To see the paths of the output nc files, do:
#'
#' `file_wx('nc', base_dir, output_nm, regex)`
#'
#' RAP/RUC GRIBS come in two resolutions, so these are processed separately and saved
#' to separate sub-directories, named in `output_nm`. This means there are two output
#' nc files per variable: the primary 13km grid series ("fine"), and the much smaller
#' 25km series ("coarse").
#'
#' `output_nm` should be a list with entries 'coarse' and/or 'fine', each containing
#' one or more subdirectory names. This is to allow users to store a small nc file
#' for relatively new data separately from one or more "archived" nc files that are
#' unlikely to change very often. The function always writes to the first subdirectory
#' named in each of the elements of `output_nm`, but all subdirectories are checked
#' when determining if a time has been processed yet.
#'
#' A date range for processing can be specified by passing either or both of `from`
#' and `to`. GRIB files outside of this date range are ignored. Different date ranges
#' can be specified for each variable, in which case `length(from)` and `length(to)`
#' should match `length(regex)`. By default, `from` and `to` are set to the the
#' earliest and latest GRIB available dates.

#' The script will be very slow on the initial run with many input GRIBs, but subsequent
#' calls to update an existing set of nc files (and JSONs) will be much faster, as
#' only the missing layers are read and copied, and the files to modify on disk are
#' relatively small.
#'
#' @param aoi geometry object passed to `grib_idx` (area of interest)
#' @param base_dir character path to parent directory of `names(output_nm)`
#' @param output_nm list of character vectors, sub-directories in `base_dir` for the nc files
#' @param regex character vector passed to `grib_idx` (layer names)
#' @param from Date or vector of them, GRIB files for all earlier times are ignored
#' @param to Date or vector of them, GRIB files for all later times are ignored
#' @param grib_dir character path to GRIB directory (default is "base_dir/grib")
#'
#' @return nothing, but possibly writes to the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
nc_update = function(aoi,
                     base_dir,
                     output_nm = list(coarse=c('coarse'), fine=c('fine')),
                     regex = .rap_regex,
                     from = NULL,
                     to = NULL,
                     grib_dir = file_wx('grib', base_dir)) {

  # new data written to the first of the directories listed in elements of `output_nm`
  var_nm = names(regex)
  output_path = file_wx('nc', base_dir, lapply(output_nm, \(x) x[1]), as.list(var_nm))

  # parse filenames of existing archive files to get times
  all_gribs = grib_list(grib_dir, dupe=FALSE)
  if( nrow(all_gribs) == 0 ) stop('no GRIB files found in ', grib_dir)

  # default start time is earliest GRIB file, default end time is latest GRIB file
  from = if( is.null(from) ) { as.Date( min(all_gribs[['posix_pred']]) ) } else { as.Date(from) }
  to = if( is.null(to) ) { as.Date( max(all_gribs[['posix_pred']]) ) } else { as.Date(to) }

  # check consistency of vector start/end arguments
  if( length(from) == 1 ) from = rep(from, length(regex))
  if( length(from) != length(regex) ) stop('"from" must have the same length as "regex" (or 1)')
  if( length(to) == 1 ) to = rep(to, length(regex))
  if( length(to) != length(regex) ) stop('"to" must have the same length as "regex" (or 1)')
  names(from) = names(regex)
  names(to) = names(regex)

  # filter to relevant GRIBs
  all_gribs = all_gribs |> dplyr::filter(date_rel >= min(from))
  if( nrow(all_gribs) == 0 ) stop('no GRIB files found at or after ', as.character(min(from)))
  all_gribs = all_gribs |> dplyr::filter(date_rel <= max(to))
  if( nrow(all_gribs) == 0 ) stop('no GRIB files found at or before ', as.character(max(to)))

  # files at different resolutions processed separately
  for( nm_res in names(output_nm) ) {

    paste0('\n\n', nm_res, ' grids...') |> cat()
    path_nc = output_path[[nm_res]]

    # identify all files at this resolution
    is_coarse = nm_res == 'coarse'
    is_included = all_gribs[['coarse']] == is_coarse
    grib_df = all_gribs[is_included,]
    if( nrow(grib_df) == 0 ) {

      cat('\nnone found')
      next
    }

    # matrix indicating for each variable (column) whether the file (row) needs to be loaded
    time_available = grib_df[['posix_pred']]
    is_new_mat = names(output_path[[nm_res]]) |> sapply(\(nm) {

      # select all times beyond cutoff date
      p = output_path[[nm_res]][[nm]]
      is_eligible = ( as.Date(time_available) >= from[nm] ) & ( as.Date(time_available) <= to[nm] )

      # of those, select all times not already processed
      time_done = time_wx(p)[['time']]
      is_eligible & !( time_available %in% time_done )

    }) |> matrix(nrow=length(time_available))

    # which files and times are missing from the nc files
    idx_new = do.call(c, apply(is_new_mat, 2, which, simplify=FALSE)) |> unique()
    time_add = time_available[idx_new]
    if( length(time_add) == 0 ) {

      cat('\nup to date')
      next
    }

    # updates done one year at a time
    yr_add = time_add |> format('%Y')
    yr_unique = yr_add |> unique()
    cat('\nprocessing', length(yr_unique), 'year(s)...\n')
    for(yr in yr_unique) {

      cat('\nyear', yr)
      file_idx = match(time_add[yr_add == yr], time_available)
      n_load = length(file_idx)

      # slow loading part starts
      if( n_load == 0 ) {

        if( all( !file.exists( unlist(input_path[[nm_res]]) ) ) ) stop('no data to write')
        r_from_gribs = input_path[[nm_res]] |> lapply(\(x) NULL)

      } else {

        # check if we only need a subset of variables
        is_pending = as.matrix(is_new_mat[file_idx,], ncol=length(var_nm)) |> apply(2, any)

        # load GRIB data
        r_from_gribs = grib_extract(grib_df,
                                    file_idx = file_idx,
                                    regex = regex[is_pending],
                                    aoi = aoi)
      }

      # slow writing part starts
      if( n_load == 0 ) { cat('\nup to date') } else {

        cat('\nupdating .nc files:')
        for( nm in names(r_from_gribs) ) {

          cat(paste0('\n\n', nm), '...')
          nc_write_chunk(r=r_from_gribs[[nm]], p=output_path[[nm_res]][[nm]])
        }
      }
    }
  }
  cat('\n')
}
