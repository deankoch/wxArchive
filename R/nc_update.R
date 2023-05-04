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
#' If `from` is supplied, the function only copies new times greater than or equal to
#' `from`. A different `from` can be specified for each variable, in which case
#' `length(from)` should match `length(regex)`. By default, `from` is set to the
#' earliest available time.
#'
#' Within each of these sub-directories, another sub-directory, 'time', is created
#' to store JSON files, one per variable (ie one per nc file). These hold the indices
#' of NA layers, and the observed times, for quicker loading later on. To see their
#' output paths. do:
#'
#' `file_wx('index', base_dir, output_nm, regex)`
#'
#' The script will be very slow on the initial run with many input GRIBs, but subsequent
#' calls to update an existing set of nc files (and JSONs) will be much faster, as
#' only the missing layers are read and copied, and the files to modify on disk are
#' relatively small.
#'
#' Note that total precipitation `.rap_regex['pcp_total']`is not available for the
#' first decade or so. By default the function creates a dummy nc file (empty of
#' data) for this period to avoid unnecessarily loading GRIBs to look for it. This
#' behaviour can be switched off with `make_dummy=FALSE`.
#'
#' @param aoi geometry object passed to `grib_idx` (area of interest)
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param output_nm list of character vectors, sub-directories in `base_dir` for the nc files
#' @param regex character vector passed to `grib_idx` (layer names)
#' @param n_chunk number of files to load before saving intermediate results to disk
#' @param memory_limit integer (GB) maximum memory usage passed to `grib_extract`
#' @param make_dummy logical, indicates to omit "pcp_total" for early years (see details)
#' @param from POSIXct time or vector of them, GRIB files for this and all earlier times are ignored
#'
#' @return nothing, but possibly writes to the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
nc_update = function(aoi,
                     base_dir,
                     output_nm = list(coarse=c('coarse'), fine=c('fine')),
                     regex = .rap_regex,
                     n_chunk = 5e3,
                     memory_limit = 8L,
                     make_dummy = FALSE,
                     from = NULL) {

  # pass lists to file_wx to loop over the variable/directory and set names properly
  var_nm = names(regex)
  grib_dir = file_wx('grib', base_dir)
  input_path = file_wx('nc', base_dir, as.list(output_nm), as.list(var_nm), make_dir=TRUE)

  # new data written to the first of the files listed in `output_nm`
  output_path = file_wx('nc', base_dir, lapply(output_nm, \(x) x[1]), as.list(names(regex)))
  output_json = file_wx('index', base_dir, lapply(output_nm, \(x) x[1]), as.list(names(regex)))

  # parse filenames of existing archive files to get times
  all_gribs = grib_list(grib_dir, dupe=FALSE)
  if( nrow(all_gribs) == 0 ) stop('no GRIB files found in ', grib_dir)

  # set default start times
  if( is.null(from) ) { from = min(all_gribs[['posix_pred']]) } else {

    all_gribs = all_gribs |> dplyr::filter(posix_pred >= min(from))
    if( nrow(all_gribs) == 0 ) stop('no GRIB files found at or after', as.character(min(from)))
  }

  # check consistency of two arguments
  if( length(from) == 1 ) from = rep(from, length(regex))
  if( length(from) != length(regex) ) stop('"from" must have the same length as "regex" (or 1)')
  names(from) = names(regex)

  # files at different resolutions processed separately
  for( nm_res in c('coarse', 'fine') ) {

    cat('\n\n', nm_res, 'grids')
    path_nc = output_path[[nm_res]]

    # identify all files at this resolution
    is_coarse = nm_res == 'coarse'
    is_included = all_gribs[['coarse']] == is_coarse
    grib_df = all_gribs[is_included,]
    if( nrow(grib_df) == 0 ) {

      cat('\nnone found!')
      next
    }

    # create empty pcp_total file for times where it was not offered
    is_missing = ifelse(is.null(path_nc[['pcp_total']]), TRUE, !file.exists(path_nc[['pcp_total']]))
    make_dummy = make_dummy & ('pcp_total' %in% var_nm) & is_missing
    if( make_dummy ) my_dummy_nc(path_nc[['pcp_total']], grib_df, aoi)

    # matrix indicating for each variable (column) whether the file (row) needs to be loaded
    is_new_mat = names(input_path[[nm_res]]) |> sapply(\(nm) {

      # select all times beyond cutoff date
      p = input_path[[nm_res]][[nm]]
      is_eligible = grib_df[['posix_pred']] >= from[nm]

      # of those, select all times not already processed
      time_done = time_wx(p)[['time']]
      is_eligible & !( grib_df[['posix_pred']] %in% time_done )
    })

    # load chunks in a loop until nothing left to load
    is_new = TRUE
    while( any(is_new) ) {

      # which files are missing from the nc
      is_new = is_new_mat |> apply(1, any)

      # select the first n_chunk for loading
      file_idx = which(is_new) |> head(n_chunk)
      n_new = sum(is_new)
      n_load = length(file_idx)

      t1 = proc.time()

      # slow loading part starts
      if( n_load == 0 ) {

        if( all( !file.exists( unlist(input_path[[nm_res]]) ) ) ) stop('no data to write')
        r_from_gribs = input_path[[nm_res]] |> lapply(\(x) NULL)
        t2 = proc.time()

      } else {

        # check if we only need a subset of variables
        is_pending = as.matrix(is_new_mat[file_idx,], ncol=length(var_nm)) |> apply(2, any)
        msg_subset = paste(names(is_pending)[is_pending], collapse=', ') |> paste('from')
        cat('\nloading', msg_subset,
            paste0(ifelse(n_load < n_new, 'first ', ''), n_load),
            'of', n_new , 'gribs...')

        # load GRIB data
        r_from_gribs = grib_extract(grib_df,
                                    file_idx = file_idx,
                                    regex = regex[is_pending],
                                    aoi = aoi,
                                    memory_limit = memory_limit)

        # finished the slow part
        t2 = proc.time()
        cat('finished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
      }

      # slow writing part starts
      if( n_new == 0 ) { cat('\nup to date \U2713') } else {

        cat('\nupdating .nc files:')
        for( nm in names(r_from_gribs) ) {

          cat('\n\n', nm, '...')
          nc_write(r = r_from_gribs[[nm]], p = output_path[[nm_res]][[nm]])
        }

        # finished the slow part
        t3 = proc.time()
        cat('\nfinished in', round((t3-t2)['elapsed'] / 60, 2), 'minutes.\n')
      }

      # update counters
      is_new_mat[file_idx,] = FALSE
    }
  }
}
