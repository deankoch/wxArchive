#' Get or set fields in a JSON file on disk representing attributes of a time series
#' 
#' This uses `jsonlite` to write information about the observed times in the netCDF
#' files (`overwrite=TRUE`), or else reads it back to produce the following in a list:  
#' 
#' * "na": index of layers containing any number of NAs
#' * "time" : all times as POSIXct strings (in same order as layers)
#' * "time_obs" : elements of "time" whose layers are not listed "na"
#' * "time_na" : elements of "time" whose layers are listed "na"
#' 
#' All times are assumed to be in the UTC time zone.
#' 
#' In read mode (`overwrite=FALSE`), the function returns all fields in the JSON as
#' a list, along with "time_obs" and "time_na", which are computed from "na" and "time".
#' This happens without checking the contents of the nc file. If the file
#' doesn't exist in read mode, the function returns NULL.
#' 
#' In write mode (`overwrite=TRUE`), only the fields "na" and "time" are written to
#' the JSON. This should leave any other existing fields alone unless there is a name
#' collision, or the structure is too complicated for a return trip through `fromJSON`
#' and `toJSON`. If `lazy=TRUE` the function will only write the JSON if it doesn't
#' already exist, and otherwise does not modify any existing ones.
#' 
#' If a SpatRaster is passed to in `r`, the function will append the NA index and times
#' in `r` to the existing lists for the raster at `nc_path`. This is useful when
#' appending a small number of new times to a large existing time series, where it
#' would be slow and redundant to check for NAs in existing layers.
#' 
#' Pass a vector to `nc_path` in read mode to get results for each file back in a list.
#' If `ch=TRUE`, the function treats these files as chunks of a single time series, and
#' the list output is collapsed to into a single list for all of the files - 'time' lists
#' times appearing in any of the files; 'time_obs' lists the times at which a non-NA layer
#' is found in at least one of the files; 'time_na' lists the times at which every
#' appearance in every file is an NA layer; and 'na' gives the index of `time_na` in `time`.
#' 
#' By default the JSON is located in sub-directory 'time' and has the same name
#' as the nc (except for file extension). Change this with `json_path`.
#' 
#' @param nc_path character path to the (.nc) time series data file
#' @param json_path character path to the JSON (.json) attributes file
#' @param overwrite logical indicating to compute fields and write to disk
#' @param lazy logical indicating to not modify any existing JSONs
#' @param ch logical indicates to collapse results for chunked files into a single list
#' @param r SpatRaster to check for appended times (avoids checking `rast(nc_path)`)
#'
#' @return a list, the contents of the JSON (after updating, if `overwrite=TRUE`) 
#' @export
my_nc_attributes = function(nc_path, json_path=NULL, overwrite=FALSE, lazy=FALSE, ch=FALSE, r=NULL) {
  
  # by default we should have like-named JSONs in subdirectory "time"
  json_expected = nc_path |> tools::file_path_sans_ext() |> basename() |> paste0('.json')
  if( is.null(json_path) ) json_path = file.path(dirname(nc_path), 'time', json_expected)
  
  # length check
  if( length(json_path) != length(nc_path) ) stop('argument lengths did not match')
  is_added = !is.null(r) & file.exists(json_path)
  
  # loop for vectorized case
  list_out = seq_along(json_path) |> lapply(\(i) {
    
    # read mode
    if( !overwrite ) {
      
      if( !file.exists(json_path[i]) ) { NULL } else {
        
        # load associated JSON as list
        out_list = json_path[i] |> readLines() |> jsonlite::fromJSON()
        
        # create `POSIXct` times
        nm_posix = c('time', 'time_na', 'time_obs')
        is_na = seq_along(out_list[['time']]) %in% out_list[['na']] 
        out_list[['time_na']] = out_list[['time']][is_na]
        out_list[['time_obs']] = out_list[['time']][!is_na]
        out_list[nm_posix] = out_list[nm_posix] |> 
          sapply(\(t) { if( length(t) > 0 ) as.POSIXct(t, tz='UTC') else list() } )
        
        out_list
      }

    } else {
      
      # write mode: make the directory if necessary
      p = json_path[i]
      if( !file.exists(dirname(p)) ) dir.create(dirname(p))
      p_exists = file.exists(p)
      
      # create the JSON fields unless we're lazy and the file exists already
      if( !p_exists | !lazy | is_added[i] ) {
        
        # load existing JSON attributes as list (if any)
        out_list = list()
        if( p_exists ) out_list = my_nc_attributes(nc_path[i], overwrite=FALSE, lazy=FALSE)
        
        # if r not supplied as argument, open nc file and compute its NA index
        cat('\nindexing NA layers and writing to', p)
        if( !is_added[i] ) r = nc_path[i] |> terra::rast()
        time_r = terra::time(r)
        n_existing = length(out_list[['time']])
        n_na = r |> terra::global('isNA') |> as.matrix() |> as.numeric()
        idx_r_na = which( n_na > 0 )
        
        # build output NA index and times
        if( is_added[i] ) {
         
          # sanity check for append mode 
          is_dupe = out_list[['time']] %in% time_r
          msg_dupe = out_list[['time']][is_dupe] |> paste(collapse=', ')
          if( any(is_dupe) ) warning( paste('time(s) found in both r and the nc file:', msg_dupe) )
          
          # times for every layer and index of NAs
          out_list[['na']] = c(out_list[['na']], idx_r_na + n_existing)
          out_list[['time']] = c(out_list[['time']], time_r)
          
        } else {
          
          # if r was not supplied or JSON not found these are already complete
          out_list[['na']] = idx_r_na
          out_list[['time']] = time_r
        }

        # omit names added by this function
        out_list = out_list[ !( names(out_list) %in% c('time_na', 'time_obs') ) ]
        
        # export to JSON and write changes to disk
        out_list |> jsonlite::toJSON(pretty=TRUE) |> writeLines(p)
        cat(' \U2713')
      }
      
      # return read mode output for the existing JSON
      my_nc_attributes(nc_path[i], overwrite=FALSE, lazy=FALSE)
    }
    
  }) |> stats::setNames( names(nc_path) )
  
  # remove NULL results from missing files
  is_missing = sapply(list_out, is.null) 
  if( all(is_missing) ) { return(NULL) } else { list_out = list_out[!is_missing] } 
  
  # collapse lists representing a time series in chunks
  if( ch & length(list_out) > 1 ) {
    
    # merge and sort all indices
    time_all_list = lapply(list_out, \(x) x[['time']])
    time_obs_list = lapply(list_out, \(x) x[['time_obs']])
    time_all = do.call(c, time_all_list[sapply(time_all_list, length) > 0]) |> unique() |> sort()
    time_obs = do.call(c, time_obs_list[sapply(time_obs_list, length) > 0]) |> unique() |> sort()
    na_all = which(!(time_all %in% time_obs))
    
    # simplified list
    list_out = list(na = na_all,
                    time = time_all,
                    time_obs = time_obs,
                    time_na = time_all[na_all])
  }
  
  # collapse length-1 lists
  if( length(list_out) == 1 ) list_out = list_out[[1]]
  return(list_out)
}


#' Open a subset of times from one or more NetCDF files
#' 
#' Returns a SpatRaster with time-indexed layers in chronological order,
#' loaded from one or several NetCDF files. This wraps `terra::rast`, setting
#' the `lyrs` argument by matching `times` to times listed in the JSON file(s)
#' in subdirectory "/time".  If `times` is `NULL`, the function returns all
#' times found in the first file `p[1]`.
#' 
#'`p` can be a vector of paths, indicating to match `times` to the times in each
#' file. If a time is found in multiple files, the function returns the matching
#' layer from the first file (WRT to the order in `p`), and ignores the rest.
#' The returned SpatRaster will include one layer for all `times` that could be
#' matched, with different layers possibly coming from different files.
#' 
#' This functionality allows you to split time series into temporal chunks,
#' and store them in multiple files. Note that the function assumes all files
#' in `p` contain data for the same variable (variable names are ignored).
#' 
#' The function creates the JSON file(s) as needed, and assumes that
#' all of `times` can be found in the .nc file(s) at `p`. Available times
#' for the file at `p` can be listed with `terra::rast(p) |> terra::time()`,
#' and this should match `my_nc_attributes(p)[['time']]`, assuming the JSON
#' is up to date.
#' 
#' When `preload=TRUE` the function forces `terra` to copy all data values
#' into memory before returning the SpatRaster handle. Normally `terra` loads
#' values only when they are needed, but I have found that this can create
#' problems with file operations when there are a large number of times
#' (layers) in the file. 
#' 
#' When `na_rm=TRUE`, only layers with all grid points observed (non-NA) are
#' searched. In this case times that match a layers labelled as NA are reported
#' as "unmatched" (and no NA layers are returned).
#'
#' @param p character vector path to the nc file(s)
#' @param times vector of unique POSIXct times to match in the file(s)
#' @param preload logical indicating to load values into RAM
#' @param na_rm logical indicating to return only non-NA layers
#'
#' @return a SpatRaster with `length(t)` layers
#' @export
my_nc_layers = function(p, times=NULL, preload=TRUE, na_rm=FALSE) {
  
  # filter nonexistent files
  p_valid = file.exists(p)
  if( !any(p_valid) ) stop('file(s) not found: ', paste(p, collapse=', '))
  p = p[p_valid]
  
  # loop over vectorized input
  load_all = is.null(times)
  times = unique(times)
  r_out_list = seq_along(p) |> lapply(\(x) NULL)
  for(i in seq_along(p)) {
    
      # load/create attributes JSON and find matching times in this file
      attr_i = my_nc_attributes(nc_path=p[i], overwrite=TRUE, lazy=TRUE)
      time_available_i = if( na_rm ) attr_i[['time_obs']] else attr_i[['time']]
      if( load_all ) times = time_available_i
      is_i = times %in% time_available_i
      if( any(is_i) ) {
        
        # return the matching layers (+0 forces copying into memory)
        r_out_list[[i]] = terra::rast(p[i], lyrs=match(times[is_i], attr_i[['time']]))
        if( preload ) r_out_list[[i]] = r_out_list[[i]] + 0 
       
        # remove from to-do stack
        times = times[!is_i]
      }
  }
  
  # report any unmatched times
  msg_unmatched = as.character(times, tz='UTC') |> paste(collapse=', ')
  if( length(times) > 0 ) cat('\nunmatched times:', msg_unmatched)
  
  # remove unmatched files outputs and sort by time
  r_out = r_out_list[!sapply(r_out_list, is.null)] |> terra::rast()
  r_out = r_out[[ order(terra::time(r_out)) ]]
  return(r_out)
}

#' Modify a NetCDF file by adding new time layers
#' 
#' This creates or updates an existing NetCDF file by adding the layers (times)
#' in the supplied SpatRaster `r`.
#' 
#' The function does little checking to see if the request is valid. It assumes
#' that `p` is a valid path to a NetCDF file (.nc), and that raster grid in `r`
#' is compatible with the file at `p` in terms of dimensions etc.
#' 
#' `r` must have a time for each layer (check `terra::time(r)`). Only new times
#' (those not found in the file at `p`) are copied when `append=TRUE`. Set
#' `append=FALSE` to copy all times from `r`, replacing any existing ones. This
#' retains any times in the nc file that are not found in `r`.
#' 
#' If `overwrite=FALSE` (the default), the function overwrites nothing but
#' returns the times that would be added to the file.
#' 
#' When overwriting an existing file, the function initially writes to a temporary
#' file in the same directory. At the end of the function call the old file is
#' deleted and the temporary one is renamed to take its place. If something goes
#' wrong and the function halts, this temporary file can be safely deleted.
#' 
#' After writing to the .nc file, the function creates/updates the JSON attributes
#' file in the 'time' subdirectory. 
#'
#' @param r SpatRaster with POSIXct vector `terra::time(r)`, the data to write
#' @param p character path the (.nc) time series data file to write
#' @param overwrite logical indicating to write changes to the nc file on disk
#' @param append logical indicating to only add new times to the file 
#' 
#' @return vector of POSIXct times, the layers added to the file
#' @export
my_nc_write = function(r, p, overwrite=FALSE, append=TRUE) {

  # make attributes JSON for nc file at p if it's missing
  attr_nc = my_nc_attributes(p)
  is_replacement = file.exists(p)
  if( is_replacement & is.null(attr_nc) ) attr_nc = my_nc_attributes(p, overwrite=TRUE)
    
  # collect times from input raster
  t_new = terra::time(r)
  if( is.null(t_new) ) stop('terra::time(r) must return a time for each layer in r')
  
  # copy all times in the existing file
  t_nc = attr_nc[['time']]
  is_first = is.null(t_nc)
  
  # append mode omits duplicates in r, otherwise omit duplicates in nc file
  is_new =  if( append ) !(t_new %in% t_nc) else rep(TRUE, length(t_new))  
  t_fetch = as.POSIXct(t_nc[ !(t_nc %in% t_new[is_new]) ])
  
  # if there's nothing to add we are done
  if( !any(is_new) ) {
    
    # compute index JSON if it's not there
    cat('\nup to date')
    return( as.POSIXct(integer(0)) )
  } 
  
  # +0 forces terra to load data into RAM
  cat('\nloading and sorting', sum(is_new), 'input SpatRaster layer(s)')
  r_add = r[[ which(is_new) ]] + 0
  
  # merge new and old layers in memory
  if( is_first | ( length(t_fetch) == 0 ) ) { r_out = r_add } else {
   
    cat('\nmerging with', length(t_fetch), 'existing nc layer(s)')
    r_existing = my_nc_layers(p, t_fetch, preload=TRUE)
    r_out = c(r_existing, r_add) 
  }
  
  # sort and name output layers
  r_out = r_out[[ order(terra::time(r_out)) ]]
  names(r_out) = paste0('lyr_', seq(terra::nlyr(r_out)))
  
  # update the nc and JSON files
  if( overwrite ) {
    
    # a name for the dataset pulled from the file path
    nm = basename(p) |> tools::file_path_sans_ext()
    cat('\nwriting to', p)
    
    # foolproof overwrite via tempfile
    if( is_replacement ) {
      
      # a temporary file name for the existing file
      suffix_temp = paste0('_', basename(tempfile()), '.nc')
      p_temp = file.path(dirname(p), paste0(nm, suffix_temp))
      
    } else { p_temp = p }
    
    # write result
    terra::writeCDF(r_out, p_temp, varname=nm)
    if(is_replacement) {
      
      # remove old file and rename new one to replace it
      unlink(p)
      file.rename(from=p_temp, to=p)
      p_temp = p
    }
    
    # we pass `r_add` to my_nc_attributes for speed only when all new times appear at the end 
    if( any( t_nc >= min(terra::time(r_add)) ) ) r_add = NULL
    
    # update attributes
    cat(' \U2713')
    my_nc_attributes(p_temp, overwrite=TRUE, r=r_add)
  }
  return( t_new[is_new] )
}


