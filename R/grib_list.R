#' Create/update a CSV with info on GRIB files retrieved with `archive_get`.
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
#' Call `file_wx('csv')` to get the default `csv`; and `file_wx('grib', base_dir)`
#' to get the default `grib_dir` for project directory `base_dir`.
#'
#' @param grib_dir path to directory with GRIB files
#' @param csv character, file name for the CSV (set `NULL` to use default in `file_wx`)
#' @param dupe logical, indicating to keep only the first of every duplicate time
#' @param quiet logical, suppresses console output
#'
#' @return a tibble, the contents of the CSV
#' @export
grib_list = function(grib_dir, csv=NULL, dupe=NA, quiet=FALSE) {

  # default paths from another helper function
  if( is.null(csv) ) csv = file_wx('csv') |> basename()

  if( !quiet ) cat('scanning for GRIB files in', grib_dir)
  files_new = files_existing = list.files(grib_dir)

  # define output classes in data frame with 0 rows
  empty_df = data.frame(name = character(0),
                        grib2 = logical(0),
                        posix_pred = character(0) |> as.POSIXct(),
                        posix_rel = character(0) |> as.POSIXct(),
                        hour_pred = integer(0),
                        hour_rel = integer(0),
                        date_rel = character(0) |> as.Date(),
                        size = numeric(0) |> units::set_units('MB'),
                        file = character(0),
                        coarse = logical(0))

  # compare to existing data frame loaded from CSV file on disk (if any)
  csv_path = file.path(grib_dir, csv)
  has_csv = file.exists(csv_path)
  if( has_csv ) {

    # read the CSV, setting column classes
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

  # process any new GRIBs found in storage directory
  any_changes = length(files_new) > 0
  if( any_changes ) {

    cat('\nparsing', length(files_new), 'new file name(s)')
    file_df = data.frame(file = files_new) |>
      dplyr::filter(grepl('\\.grb2*$', file)) |>
      dplyr::mutate(grib2 = grepl('\\.grb2$', file)) |>
      dplyr::mutate(size = file.size(file.path(grib_dir, file))) |>
      dplyr::mutate(size = units::set_units(units::set_units(size, bytes), megabytes))

    if(nrow(file_df) == 0) {

      cat('\nno changes detected')
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
                    file)

    # match file name prefix with expected grid dimensions
    new_grib_df[['coarse']] = new_grib_df[['name']] |> strsplit('_') |>
      sapply(\(x) tail(x, 1L) != '130')

    # merge with any existing results
    grib_df = rbind(grib_df, new_grib_df) |> dplyr::arrange(posix_pred)

  } else {

    # either the CSV is up to date or there were no GRIBs in the directory
    if( has_csv ) {

      if( !quiet ) cat('\nno changes detected')

    } else { if( !quiet ) { message('\nno gribs found in ', grib_dir) } }
  }

  # ordered by prediction time, with ties sorted to show later release times first
  grib_df = grib_df |> dplyr::arrange(posix_pred, dplyr::desc(posix_rel))

  # write changes to csv on disk
  if( !is.na(csv) & any_changes ) {

    if( !quiet ) cat('\nupdating', csv)
    write.csv(grib_df, csv_path, row.names=FALSE)
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

  # append local path to result before returning
  grib_df = grib_df |> dplyr::mutate(path = file.path(grib_dir, file))
  return(grib_df)
}
