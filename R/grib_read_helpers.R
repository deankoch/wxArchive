#' Create/update a CSV with info on GRIB files retrieved with `my_archive_getter`.
#' 
#' The CSV is saved to `file.path(base_dir, csv)`, overwriting anything already there,
#' except when `csv=NA`, in which case the function returns its usual output but writes
#' nothing to disk.
#' 
#' The data frame reports information about the GRIB files in `grib_dir` parsed from their
#' file names. Parsing is slow with tens of thousands of files, so the results are saved in
#' a CSV on disk for quick access later on. Each time the function is called it checks the
#' directory for changes, updating the CSV on disk as needed (unless `csv=NA`).
#' 
#' Output fields describe the file, the model it came from, and the grid resolution.
#' Three fields describe the (two) temporal dimensions of the data:
#' 
#'  * posix_pred : the POSIXct time at which the forecast file is valid
#'  * date_rel : the Date on which the forecast was created
#'  * hour_rel : the (integer) hour of the day when the forecast was created
#'  
#' Output rows are first sorted by `posix_pred` (chronological order), then by `date_rel`,
#' then `hour_rel`. Duplicate `posix_pred` values are dealt with according to the value of
#' `dupe`. The default `NA` returns all duplicates; `FALSE` returns only the first of every
#' duplicate (ie the most recently released file); and `TRUE` returns the complement of
#' `FALSE` (ie all duplicates except the most recent).
#' 
#' Call `wx_file('csv')` to get the default `csv`; and `wx_file('grib', base_dir)`
#' to get the default `grib_dir` for project directory `base_dir`.
#'
#' @param grib_dir path to directory with GRIB files
#' @param csv character, file name for the CSV (set `NULL` to use default in `wx_file`)
#' @param dupe logical, indicating to keep only the first of every duplicate time
#'
#' @return a tibble, the contents of the CSV
#' @export
my_archive_lister = function(grib_dir, csv=NULL, dupe=NA, quiet=FALSE) {
  
  # default paths from another helper function
  if( is.null(csv) ) csv = wx_file('csv') |> basename()
  
  if( !quiet ) cat('scanning for GRIB files in', grib_dir)
  files_new = files_existing = list.files(grib_dir)
  
  # define output classes in data frame with 0 rows
  empty_df = data.frame(name = character(0),
                        grib2 = logical(0),
                        posix_pred = character(0) |> as.POSIXct(),
                        hour_pred = integer(0),
                        hour_rel = integer(0),
                        date_rel = character(0) |> as.Date(),
                        size = numeric(0) |> units::set_units('MB'),
                        file = character(0),
                        path = character(0),
                        coarse = logical(0))
  
  # compare to existing data frame loaded from CSV file on disk (if any)
  csv_path = file.path(grib_dir, csv)
  has_csv = file.exists(csv_path)
  if( has_csv ) {
    
    # read the CSV, setting column classes and removing nonexistent files
    if( !quiet ) cat('\nreading existing files list in', csv)
    grib_df = csv_path |> read.csv(header=TRUE) |> 
      dplyr::mutate(posix_pred = as.POSIXct(posix_pred, tz='UTC')) |>
      dplyr::mutate(date_rel = as.Date(date_rel)) |>
      dplyr::mutate(size = units::set_units(size, megabytes)) |>
      dplyr::tibble()
    
    # copy file names not found in existing list, ignore the CSV itself
    files_new = files_new[ !(files_new %in% grib_df[['file']]) ]
    files_new = files_new[ !endsWith(files_new, '.csv') ]
    
  } else { grib_df = empty_df }
  
  # remove nonexistent files
  is_still_there = grib_df[['file']] %in% files_existing
  any_dropped = any(!is_still_there)
  if( any_dropped ) {
    
    cat('\nomitting', sum(!is_still_there), 'missing files \U2713')
    grib_df = grib_df[is_still_there,]
  }
  
  # process any new GRIBs found in storage directory
  any_changes = length(files_new) > 0 
  if( any_changes ) {
    
    cat('\nparsing', length(files_new), 'new file name(s)')
    file_df = data.frame(file = files_new) |>
      dplyr::filter(grepl('\\.grb2*$', file)) |>
      dplyr::mutate(grib2 = grepl('\\.grb2$', file)) |>
      dplyr::mutate(path = file.path(grib_dir, file)) |>
      dplyr::mutate(size = file.size(file.path(grib_dir, file))) |>
      dplyr::mutate(size = units::set_units(units::set_units(size, bytes), megabytes))
    
    if(nrow(file_df) == 0) {
      
      cat(ifelse(any_dropped, '\nno further', '\nno'), 'changes detected')
      return(file_df)
    }
    
    # append columns for date/time parsed from file name
    time_df = do.call(rbind, file_df[['file']] |> strsplit('[_\\.]') |> lapply(\(x) {
      
      data.frame(name = paste(x[1:2], collapse='_'),
                 date_rel = as.Date(x[3], '%Y%m%d'),
                 hour_rel = as.integer(as.integer(x[4])/1e2),
                 hour_pred = as.integer(x[5]))
    }) )
    
    # add POSIXct field to establish time zone
    date_df = time_df |>
      dplyr::mutate(date_rel_string = paste(date_rel, '00:00:00')) |>
      dplyr::mutate(posix_rel = as.POSIXct(date_rel_string, tz='UTC') + (60*60*hour_rel)) |>
      dplyr::mutate(posix_pred = posix_rel + (60*60*hour_pred))
    
    # combine all, sort into chronological order add empty columns for grid dimensions
    new_grib_df = cbind(file_df, date_df) |>
      dplyr::arrange(posix_pred) |>
      dplyr::tibble() |>
      dplyr::select(name,
                    grib2,
                    posix_pred,
                    posix_rel,
                    hour_pred,
                    hour_rel,
                    date_rel,
                    size,
                    file,
                    path)
    
    # match file name prefix with expected grid dimensions
    new_grib_df[['coarse']] = new_grib_df[['name']] |> strsplit('_') |> 
      sapply(\(x) tail(x, 1L) != '130')

    # merge with any existing results
    grib_df = rbind(grib_df, new_grib_df) |> dplyr::arrange(posix_pred)
    
  } else { 
    
    # either the CSV is up to date or there were no GRIBs in the directory
    if( has_csv ) {
      
      if( !quiet ) cat(ifelse(any_dropped, '\nno further', '\nno'), 'changes detected')
      
    } else { if( !quiet ) { message('\nno gribs found in ', grib_dir) } }
  }
  
  # ordered by prediction time, with ties sorted to show later release times first
  grib_df = grib_df |> dplyr::arrange(posix_pred, dplyr::desc(posix_rel))
  
  # write changes to csv on disk
  if( !is.na(csv) & ( any_changes | any_dropped ) ) {
    
    if( !quiet ) cat('\nupdating', csv)
    write.csv(grib_df, csv_path, row.names=FALSE)
    if( !quiet ) cat(' \U2713')
  }
  
  # remove or return only the duplicates as requested
  if( !is.na(dupe) ) {
    
    is_dupe = grib_df[['posix_pred']] |> duplicated()
    n_dupe = sum(is_dupe)
    if( !dupe ) { msg_verb = '\nremoving' } else {
      
      is_dupe = !is_dupe
      msg_verb = '\nreturning'
    }
    
    if( !quiet & any(is_dupe) ) cat(msg_verb, n_dupe, 'duplicate prediction time(s)')
    grib_df = grib_df |> dplyr::slice( which(!is_dupe) )
  }
  return(grib_df)
}


#' Find layer and sub-grid indices for calls to `terra::rast` and `terra::values`
#' 
#' This uses `terra` to read the layer names and grid dimensions from the GRIB at
#' `grib_path` (but not its data values) in order to learn (1) the index of layer
#' names matching the regular expressions in `regex`, and (2) the dimensions and
#' position of the sub-grid inscribing the area of interest `aoi`. These vectors can
#' then be passed to `terra::rast` and `terra::values` to speed up loading the data
#' (see `my_grib_read`).
#' 
#' The function returns a list with the following
#' 
#'  * 'k' = index of requested layers
#'  * 'ij' = row/column of top-left corner point in AOI sub-grid
#'  * 'dim' = dimensions of AOI sub-grid
#'  * 'r_aoi' = empty SpatRaster defining the AOI
#'  
#' If `regex=NULL`, all layers are returned. If any of the `regex` match more than
#' one layer, only the first match is returned. If `aoi=NULL`, the entire grid is
#' requested. Otherwise `aoi` should be coercible to 'SpatVector' and must overlap
#' with some part of the grid in `grib_path` (after transforming to a common CRS). 
#' 
#' `grib_path` can be a data-frame containing multiple paths, in which case the first
#' is loaded. If this load attempt fails, and `try_again=TRUE`, the function will then
#' attempt to load the second file, and so on.
#'
#' @param grib_path character (or dataframe with 'path' column). The path to the GRIB file 
#' @param regex character vector, regular expressions to match to layer names (using `grepl`)
#' @param aoi geometry object coercible to `SpatVector`, the area of interest
#' @param quiet logical indicating to suppress console output
#' @param try_again logical indicating to attempt to load another file if the first fails
#'
#' @return list of indices (see details)
#' @export
my_grib_idx = function(grib_path, regex=NULL, aoi=NULL, quiet=FALSE, try_again=FALSE) {
  
  # by default pick first file
  if(is.data.frame(grib_path)) grib_path = grib_path[['path']]
  grib_path = grib_path[1]
  
  if(!quiet) cat('\nchecking names in ', basename(grib_path))
  
  # attempt to load the file
  r = tryCatch({ grib_path |> terra::rast() |> terra::rast() }, error = function(err) err)
  if(is(r, 'error')) {
    
    if( !try_again | (length(grib_path) < 2) ) stop('failed to load ', grib_path)
    
    # on failures, attempt to load next file in the list
    return( my_grib_idx(tail(grib_path, -1), regex, aoi, quiet, TRUE) )
  }
  
  # identify layer(s) of interest
  k = names(r) |> which_lyr(regex, quiet=quiet, na.rm=TRUE)
  
  # default aoi is the whole grid
  if( is.null(aoi) ) { r_aoi = r } else { 
    
    aoi_projected = as(aoi, 'SpatVector') |> terra::project(terra::crs(r))
    r_aoi = terra::crop(terra::rast(r), aoi_projected)
  }
  
  # find the index of first column and row of aoi in outer grid 
  i = terra::rowFromY(r, terra::yFromRow(r_aoi, 1)) 
  j = terra::colFromX(r, terra::xFromCol(r_aoi, 1))
  
  # strip the inner grid of values and unneeded layers
  r_aoi = terra::rast(r_aoi[[1]])
  
  # count number of cells along each dimension then return
  return( list(k=k, ij=c(i, j), dim=dim(r_aoi)[1:2], r_aoi=r_aoi) )
}


#' GRIB loader to catch errors with missing/damaged files
#' 
#' This is a wrapper for `terra::rast` that loads only the sub-grid and
#' layers specified in `regex` and `aoi` (using `my_grib_idx`) for the GRIB
#' file at `grib_path`.
#' 
#' This catches loading errors - instead of halting, the function
#' returns the error object (and prints it as a warning if `quiet=FALSE`).
#'
#' @param grib_path character path to the GRIB file to (attempt to) load
#' @param regex character vector passed to `my_grib_idx` (layer names)
#' @param aoi geometry object passed to `my_grib_idx` (area of interest)
#' @param quiet logical indicating to suppress console output
#'
#' @return SpatRaster containing the requested layers and sub-grid
#' @export
my_grib_read = function(grib_path, regex=NULL, aoi=NULL, quiet=FALSE) {
  
  # attempt to load the file
  t_mat = tryCatch({
    
    # find indices to load
    r_info = my_grib_idx(grib_path, regex=regex, aoi=aoi, quiet=quiet)
    
    # copies data as matrix
    terra::rast(grib_path, lyrs = r_info[['k']]) |> 
      terra::values(row = r_info[['ij']][1],
                    nrows = r_info[['dim']][1],
                    col = r_info[['ij']][2],
                    ncols = r_info[['dim']][2])
    
  }, error = function(err) err)
  
  if(is(t_mat, 'error') & !quiet) {
    
    message('\nfailed to load grib ', t_path, '\n')
    warning(as.character(t_mat))
  }
  
  return(t_mat)
}


#' Extract the specified area of interest (sub-grid) and variables (layers) from GRIB archive
#'
#' This opens (the subset `file_idx` of) GRIB files listed in `grib_df`, loading the
#' sub-grid and layers specified by `aoi` and `regex` (respectively). Results are returned
#' as a list of SpatRasters, one for each variable in `regex`, each with a layer for every
#' file (ie time).
#' 
#' With a large number of GRIBs, memory usage can get high and everything slows down
#' overall. I recommend loading at most 100s to 1000s of files at once, and merging results
#' as you go. `memory_limit` sets up a limit (in GB) as a sanity check.
#' 
#' The function calls `my_grib_idx` on `grib_df` (and not the subset in `file_idx`)
#' to initially get basic grid information like dimensions, so make sure all the GRIBs
#' in `grib_df` use the same grid.
#'
#' @param grib_df data frame, a subset of the rows by `my_archive_lister`
#' @param file_idx integer vector, the rows of `grib_df` to load
#' @param regex character vector passed to `my_grib_idx` (layer names)
#' @param aoi geometry object coercible to `SpatVector`, the area of interest
#' @param memory_limit integer (in GB) the maximum memory usage allowed of the output
#'
#' @return a list of multi-layer SpatRasters (one per variable in `regex`)
#' @export
my_grib_extract = function(grib_df, file_idx=NULL, regex=.rr_regex, aoi=NULL, memory_limit=8L) {
  
  # by default open all files in grib_df
  if(is.null(file_idx)) file_idx = nrow(grib_df) |> seq()
  
  # get grid dimensions by looking at first file
  r_info = my_grib_idx(grib_df, regex, aoi=aoi, quiet=TRUE, try_again=TRUE)
  
  # sanity check for memory usage
  n_space = prod(r_info[['dim']])
  n_var = length(regex)
  n_file = length(file_idx)
  est_gb = estimate_memory(r_info[['dim']], n_var, n_file)
  if(est_gb > memory_limit) stop('memory limit exceeded.')
  
  # initialize storage matrices for numeric data from all times
  z_out = lapply(regex, \(x) matrix(NA_real_, n_space, n_file))
  
  # loop over files, copying only the requested variables to storage
  cat('\nloading', n_file, 'files...\n')
  if(n_file > 1) pb = txtProgressBar(max=n_file, style=3)
  for(t in seq(n_file)) {
    
    # progress bar in console updated with current filename
    if(n_file > 1) setTxtProgressBar(pb, t)
    t_path = grib_df[['path']][file_idx[t]]
    cat('\r|=', basename(t_path), '')
    
    # read in the requested layers, handling load failures 
    t_mat = my_grib_read(t_path, regex=regex, aoi=aoi, quiet=TRUE)
    
    # copy all data in loop over variables
    if( !is(t_mat, 'error') ) {
      
      # loop over raster layer names, filling columns of output matrix
      j_out = colnames(t_mat) |> which_lyr(regex, na.rm=FALSE, quiet=TRUE)
      for(nm in names(j_out)) z_out[[nm]][, t] = t_mat[, j_out[nm]]
    }
  }
  if(n_file > 1) close(pb)
  
  # make a template raster for output
  r_template = my_grib_idx(grib_df, aoi=aoi, try_again=TRUE, quiet=TRUE)[['r_aoi']] |> 
    terra::rast(nlyrs=length(regex))
  
  # initialize SpatRaster list in memory, split by variable name
  r_out = r_template |> 
    lapply(\(r) rast(r, nlyrs=n_file)) |> 
    stats::setNames(nm=names(regex))
  
  # loop over variables, copying matrix data to SpatRaster
  idx_write = seq(n_space)
  cat('\ncopying', n_file, 'layers to', paste(names(regex), collapse=', '), '... ')
  for(v in seq_along(r_out)) {
    
    # assign timestamps then data values in a loop over layers (times)
    terra::time(r_out[[v]]) = grib_df[['posix_pred']][file_idx]
    names(r_out[[v]]) = terra::time(r_out[[v]])
    
    # assign values to raster cells
    terra::set.values(r_out[[v]], cells = idx_write, values = z_out[[v]][])
  }
  
  cat('done.\n')
  return(r_out) 
}
