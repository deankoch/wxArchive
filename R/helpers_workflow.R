#' Create or update a NetCDF version of a local RAP archive
#' 
#' Most GRIB files fetched by `archive_update` are around 20-30 MB compressed,
#' so loading tens of thousands at once is extremely slow. To make the time series
#' easier to work with, this function extracts the subset that we need and saves it in
#' a format that can be loaded more quickly.
#' 
#' The function opens the GRIB files in the time series created by `archive_update`
#' and merges the data into a single (multi-layer) NetCDF file for each of the variables
#' selected by `regex`. Only the sub-grid overlapping with `aoi` is copied.
#' 
#' The output files are given base names `names(regex)`, and extension '.nc', and
#' have layers named after the POSIXct time (as character, in GMT) at which the
#' forecast is valid. To see the paths of the output nc files, do:
#' 
#' `wx_file('nc', base_dir, output_nm, regex)`
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
#' The default `append=TRUE` only copies data from files whose time does not already
#' appear in the existing files listed in `output_nm`. With `append=FALSE`, the function
#' copies all times greater than or equal to `from`, overwriting any existing times in
#' the output nc files (any existing times before `from` will remain).
#' 
#' Within each of these sub-directories, another sub-directory, 'time', is created
#' to store JSON files, one per variable (ie one per nc file). These hold the indices
#' of NA layers, and the observed times, for quicker loading later on. To see their
#' output paths. do:
#' 
#' `wx_file('index', base_dir, output_nm, regex)`
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
#' @param aoi geometry object passed to `my_grib_idx` (area of interest)
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param output_nm list of character vectors, sub-directories in `base_dir` for the nc files
#' @param regex character vector passed to `my_grib_idx` (layer names)
#' @param n_chunk number of files to load before saving intermediate results to disk
#' @param memory_limit integer (GB) maximum memory usage passed to `my_grib_extract`
#' @param make_dummy logical, indicates to omit "pcp_total" for early years (see details)
#' @param from POSIXct time or vector of them, GRIB files for this and all earlier times are ignored
#'
#' @return nothing, but possibly writes to the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
my_update_nc = function(aoi,
                        base_dir,
                        output_nm = list(coarse=c('coarse'), fine=c('fine')),
                        regex = .rap_regex,
                        n_chunk = 5e3,
                        memory_limit = 8L,
                        make_dummy = FALSE,
                        from = NULL,
                        append = TRUE) {
  
  # pass lists to wx_file to loop over the variable/directory and set names properly
  var_nm = names(regex)
  grib_dir = wx_file('grib', base_dir)
  input_path = wx_file('nc', base_dir, as.list(output_nm), as.list(var_nm), make_dir=TRUE)
  
  # new data written to the first of the files listed in `output_nm`
  output_path = wx_file('nc', base_dir, lapply(output_nm, \(x) x[1]), as.list(names(regex)))
  output_json = wx_file('index', base_dir, lapply(output_nm, \(x) x[1]), as.list(names(regex)))
  
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
    
    cat('\n\nchecking', nm_res, 'grids')
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
      if( append ) {
        
        # of those, select all times not already processed
        t_done = my_nc_attributes(p, ch=TRUE)[['time']]
        is_eligible = is_eligible & !( grib_df[['posix_pred']] %in% t_done ) 
      } 
      is_eligible
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
        r_from_gribs = my_grib_extract(grib_df,
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
          my_nc_write(r = r_from_gribs[[nm]], 
                      p = output_path[[nm_res]][[nm]],
                      overwrite = TRUE,
                      append = append)
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
#' confusion downstream (list them with `wx_file('nc', base_dir, output_nm)`)
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
my_pcp_total = function(base_dir,
                        pcp_nm = 'pcp',
                        input_nm = list(coarse='coarse', fine='fine'),
                        output_nm = c('coarse', 'fine')) {
  
  # paths to expected inputs
  pcp_var_nm = list('pcp_total', 'pcp_large', 'pcp_small')
  input_path = wx_file('nc', base_dir, input_nm, var_nm=pcp_var_nm)

  # sanity check for pcp_nm
  nm_var_all = stats::setNames(nm=names(input_path[[1]]))
  if( pcp_nm %in% nm_var_all ) stop('pcp_nm cannot be an existing variable name')
  
  # new files to write
  output_nc_path = wx_file('nc', base_dir, as.list(output_nm), var_nm=pcp_nm)

  # find time coverage of each variable at both resolutions
  cat('\nscanning available times for', paste(nm_var_all, collapse=', '))
  var_info = lapply(input_path, \(r) lapply(r, \(p) my_nc_attributes(p, ch=TRUE)) )

  # loop over resolutions
  for( nm_res in output_nm ) {
    
    cat('\n\nchecking', nm_res, 'grids')
    t1 = proc.time()
    
    # check for available times by variable
    t_all = do.call(c, lapply(var_info[[nm_res]], \(v) v[['time']])) |> unique() |> sort()
    is_obs = var_info[[nm_res]] |> lapply(\(v) t_all %in% v[['time_obs']])
    is_relevant = is_obs[['pcp_large']] & is_obs[['pcp_small']] & !is_obs[['pcp_total']]
    time_relevant = t_all[is_relevant]
    
    # check for existing data (creating JSON as needed)
    p_out = output_nc_path[[nm_res]]
    time_done = character(0) |> as.POSIXct()
    if( file.exists(p_out) ) time_done = my_nc_attributes(p_out, overwrite=TRUE, lazy=TRUE)[['time']]
 
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
      r_large = input_path[[nm_res]][['pcp_large']] |> my_nc_layers(time_add) 
      r_small = input_path[[nm_res]][['pcp_small']] |> my_nc_layers(time_add)
      r_add = r_large + r_small
    }
    
    # handle NAs encoded as negatives
    cat('\nchecking for negatives... ')
    is_na = c(terra::global(r_add, 'min')[[1]]) < 0
    n_na = sum(is_na)
    if( n_na > 0 ) {

      r_add[seq(nrow(r_add)), seq(ncol(r_add)), which(is_na)] = NA
      cat('omitting', n_na, 'NA layers')
      
    } else { cat('none found') } 
    
    # append to existing data file (or create the file and write to it)
    r_add |> my_nc_write(p_out, overwrite=TRUE, append=TRUE)
    t2 = proc.time()
    cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
  }
}


#' Create/update a NetCDF archive of resampled coarse resolution grids
#' 
#' This resamples the coarse resolution time series to get data at fine resolution,
#' writing results to nc files in subdirectory `output_nm`. Only those times found in
#' the coarse series but not the fine series are included.
#' 
#' Specify the coarse and fine sub-directories as named elements in list `input_nm`. Each
#' should be a character vector of sub-directories to check for data on `var_nm`. When
#' multiple sub-directories contain data for a given time, the function extracts the data
#' from the first (according to the order in `input_nm`)
#' 
#' `var_nm` can be a list of character vectors, where each vector gives a set of variable
#' names to be treated as equivalent (eg 'pcp' and 'pcp_total', which are the same variable,
#' compiled from different sources). The output variable is assigned the corresponding name
#' from `var_nm`. If `var_nm` is unnamed, the function names it using the first element of
#' each vector.
#' 
#' The sub-directory `output_nm` will be created if it doesn't exist. If the file in
#' `output_nm` exists already, only those times available in `input_nm`, but not already
#' present in `output_nm`, are appended.
#' 
#' Coarse resolution grids are by default resampled by bilinear averaging (with GDAL).
#' Change this by setting named arguments in `...` to pass to `terra::project`.
#'
#' @param var_nm character vector or list of them, names of the variables to process
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm list of character vectors, namess of the source file subdirectories in `base_dir`
#' @param output_nm character name for the sub-directory to create.
#' @param SpatRaster optional target grid for resampling
#' @param ... additional arguments passed to `terra::project`
#'
#' @return nothing, but possibly modifies the nc and JSON files in `file.path(base_dir, output_nm)`
#' @export
my_resample = function(var_nm,
                       base_dir, 
                       input_nm = list(coarse='coarse', fine='fine'), 
                       output_nm = 'coarse_resampled',
                       r_fine = NULL,
                       from = NULL,
                       append = TRUE,
                       ...) {

  # paths to expected inputs and outputs 
  input_nc = wx_file('nc', base_dir, as.list(input_nm), as.list(var_nm))
  var_nm = var_nm |> stats::setNames(nm=names(input_nc[[1]]))
  output_nc = wx_file('nc', base_dir, output_nm, as.list(names(var_nm)), make_dir=TRUE)

  # find time coverage of each variable at both resolutions
  cat('\nscanning available times for', paste(names(var_nm), collapse=', '))
  var_info = lapply(input_nc, \(r) lapply(r, \(p) my_nc_attributes(p, ch=TRUE)) )
  
  # set default start times and check consistency of two arguments
  if( is.null(from) ) from = var_info[['coarse']] |> sapply(\(v) min(v[['time_obs']]) )
  if( length(from) == 1 ) from = rep(from, length(var_nm))
  if( length(from) != length(var_nm) ) stop('"from" must have the same length as "var_nm" (or 1)')
  names(from) = names(var_nm)
  
  # loop over variable names to check for existing resampled data
  time_add = stats::setNames(nm=names(var_nm)) |> lapply(\(nm) {

    # copy existing output times (or nothing if overwriting)
    time_done = if( !append | !file.exists(output_nc[[nm]]) ) as.POSIXct(integer(0)) else { 
      
      #open JSON or create it as needed
      my_nc_attributes(output_nc[[nm]], overwrite=TRUE, lazy=TRUE)[['time_obs']]
    }

    # find times in coarse series that haven't been processed yet
    time_coarse = var_info[['coarse']][[nm]][['time_obs']]
    time_fine = var_info[['fine']][[nm]][['time_obs']]
    time_pending = time_coarse[ !( time_coarse %in% c(time_fine, time_done) ) ]
    
    # filter to supplied start time (default includes all)
    time_pending = time_pending[ time_pending >= from[nm] ]
  })
  
  # finished if no updates to write
  needs_update = sapply(time_add, length) > 0
  if( !any(needs_update) ) {
    
    cat('\nall variables are up to date \U2713')
    return(invisible())
  }
 
  # a template SpatRaster at fine resolution from first file
  if( is.null(r_fine) ) r_fine = input_nc[['fine']][[1]][1] |> terra::rast(lyrs=1) |> terra::rast()
  
  # loop over variables
  for(nm in names(var_nm)) {
    
    # skip variables that are up to date
    if( needs_update[nm] ) {
      
      cat('\n\nprocessing', nm, '...')
      t1 = proc.time()
      
      # `project` calls GDAL
      cat('\nresampling', length(time_add[[nm]]), 'coarse resolution layers')
      r_resample = input_nc[['coarse']][[nm]] |> 
        my_nc_layers(time_add[[nm]]) |>
        terra::project(r_fine)
        #terra::project(r_fine, ...)
      
      # check for problem layers
      n_na = r_resample |> terra::global('isNA') |> as.matrix() |> as.numeric()
      if( any(n_na == terra::ncell(r_resample)) ) stop('one or more layers had all NA cells')
      
      # clean up NAs around borders with moving window mean
      is_na = n_na > 0
      sanity_i = 1
      while( any(is_na) & (sanity_i < 1e3) ) {
        
        r_resample[[which(is_na)]] = r_resample[[which(is_na)]] |> 
          terra::focal(fun=mean, na.rm=TRUE, na.policy='only')
        
        # update index of NA layers
        n_na[is_na] = r_resample[[is_na]] |> terra::global('isNA') |> as.matrix() |> as.numeric()
        is_na = n_na > 0
        sanity_i = sanity_i + 1
      }
      
      # append to existing data file (or create the file and write to it)
      r_resample |> my_nc_write(output_nc[[nm]], overwrite=TRUE, append=append)
      
      t2 = proc.time()
      cat('\nfinished in', round((t2-t1)['elapsed'] / 60, 2), 'minutes.\n')
      
    } else { cat(paste0('\n', nm), 'is up to date \U2713') }
  }
}



