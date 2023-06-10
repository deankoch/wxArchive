#' Aggregate a gridded spatio-temporal dataset over user-specified regions
#'
#' This takes a set of input files (NetCDF or GeoTIFF) and applies the aggregation
#' function `fun` to its values at each time-step. If `poly_in` is supplied, the
#' aggregation happens separately over polygons.
#'
#' The function writes its output to CSV files in the output directory `output_nm`,
#' with each variable saved to its own sub-directory. A separate file is written for
#' each year, and results for different polygons are written different columns of the
#' output.
#'
#' For consistency with earlier datasets in the workflow, this function also
#' writes JSON files indexing the dates in each CSV file (in subdirectories named "time").
#' Note also that sub-directories of `output_nm` will use the suffix '.nc' (even though
#' we aren't dealing with NetCDF at this point)
#'
#' By default `from` is set to to 10 days before latest existing date in the output,
#' and `to` is set to latest available date in `input_nm`. If any part of a year is
#' included in this date range, then the whole year is processed (and the output file
#' is overwritten completely)
#'
#' Non-NULL `fun` is passed directly to `terra::extract` to specify the aggregation
#' function. Otherwise the function looks for variable names given in snake case (with
#' '_' separating parts), and assigns the last part to `fun`. For example if `fun=NULL`
#' then variables named 'tmp_max' or 'tmp_mean_max' would both prompt `max` for
#' aggregation. If a variable name cannot be split in this way, the function reverts to
#' the default `fun='mean'`.
#'
#' @param base_dir path to parent directory of GRIB storage subfolder
#' @param input_nm character, the sub-directory name of input to process
#' @param output_nm character, the sub-directory name for output
#' @param var_nm character vector, the variable names to process
#' @param poly_in polygon(s) to aggregate over, any object acceptable by `sf::st_geometry`
#' @param fun character, name of the function to apply for aggregation
#' @param from Date, the first date of the sequence to process
#' @param to Date, the last date of the sequence to process
#' @param file_ext character, the input file extension to load, either '.tif' or '.nc'
#'
#' @return returns nothing but writes CSV files to `output_nm`
#' @export
nc_aggregate_space = function(base_dir,
                              input_nm = .nm_down,
                              output_nm = .nm_export,
                              var_nm = .var_daily,
                              poly_in = NULL,
                              fun = NULL,
                              from = NULL,
                              to = NULL,
                              file_ext = '.tif') {

  # var_nm is list to ensure input/output paths are lists too
  input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
  output_nc = file_wx('nc', base_dir, output_nm, as.list(var_nm), make_dir=TRUE)

  # check for available input and output files
  input_var_info = lapply(input_nc, \(p) time_wx(p, file_ext=file_ext))
  output_var_info = lapply(output_nc, \(p) time_wx(p, file_ext='csv'))

  # TODO: loop over variables
  for(v in seq_along(var_nm)) {

    paste0('\n\nprocessing ', var_nm[[v]], '...\n') |> cat()

    # default function last part of variable name (or else mean)
    fun_v = fun
    if( is.null(fun_v) ) {

      nm_parts = var_nm[[v]] |> strsplit('_') |> unlist()
      fun_v = if( length(nm_parts) > 1 ) nm_parts[[ length(nm_parts) ]] else 'mean'
    }

    # user arguments (or NULL)
    from_v = from
    to_v = to

    # copy existing dates
    time_out = output_var_info[[v]][['time_obs']] |> as.Date()
    time_in = input_var_info[[v]][['time_obs']] |> as.Date()
    if( is.null(time_in) )  {

      cat('\nup to date')
      next
    }

    # set default starting/ending times
    is_initial = length(time_out) == 0
    if( is.null(to_v) ) to_v = max(time_in)
    if( is.null(from_v) ) {

      # on first call this writes everything
      from_v = min(time_in)

      # subsequently the default start time is 10 days before latest time
      if( !is_initial ) from_v = as.POSIXct(max(time_out)) - ( 10 * (60*60*24) )
    }

    # silently fix invalid start/end times
    from_v = as.Date(from_v)
    to_v = as.Date(to_v)
    if( from_v < min(time_in) ) from_v = min(time_in)
    if( to_v > max(time_in) ) to_v = max(time_in)

    # filter available times to requested range
    time_req = time_in[ ( time_in >= from_v ) & ( time_in <= to_v ) ]
    if( length(time_in) == 0 ) {

      cat('nothing to write\n')
      next
    }

    # find years available, and those in requested range
    yr_v = format(time_in, '%Y', tz=tz)
    yr_unique = format(time_req, '%Y', tz=tz) |> unique()

    # loop over years
    for( y in seq_along(yr_unique) ) {

      yr_y = yr_unique[y]
      paste0('\nloading year ', yr_y, '... ') |> cat()
      time_y = time_in[yr_v == yr_y]

      # load data for the year and CRS
      r_y = input_nc[[v]] |> nc_layers(times=time_y, file_ext='.tif')
      crs_r = terra::crs(r_y) |> sf::st_crs()

      # set default polygon or split input into list of polygons
      if( is.null(poly_in) ) poly_in = terra::ext(r_y) |> sf::st_bbox() |> sf::st_as_sfc()
      poly_split = poly_in |> sf::st_geometry() |> sf::st_transform(crs=crs_r)
      poly_nm = paste0('poly_', seq_along(poly_split))

      # split input raster into segments and clear source from RAM
      cat('cropping... ')
      r_split = poly_split |> lapply(\(p) terra::crop(r_y, p, snap='out'))
      rm(r_y)
      gc()

      # apply aggregation function over polygons
      cat('aggregating... ')
      r_stat_split = Map(\(x, p) terra::extract(x = x,
                                                y = as(p, 'SpatVector'),
                                                fun = fun_v,
                                                raw = TRUE,
                                                touches = TRUE,
                                                ID = FALSE), x=r_split, p=poly_split)

      # clear raster segments from RAM
      cat('done')
      rm(r_split)
      gc()

      # copy the series from each polygon (data.frames -> list of vectors -> data.frame)
      df_v = r_stat_split |> lapply(\(x) unlist(x, use.names=FALSE)) |> data.frame()
      df_v = data.frame(date=time_y) |> cbind(df_v) |> stats::setNames(c('date', poly_nm))

      # write to disk (overwrites any existing CSV by this name)
      if( !dir.exists(output_nc[[v]]) ) dir.create(output_nc[[v]])
      p_csv = output_nc[[v]] |> file.path( paste0(var_nm[[v]], '_', yr_y, '.csv') )
      cat('\nwriting CSV to', p_csv)
      df_v |> write.csv(p_csv, row.names=FALSE, quote=FALSE)

      # overwrite time index
      write_time_json(nc_path=p_csv, r=time_y)
      cat('\n')
    }
  }
}
