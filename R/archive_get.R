#' Download batches of GRIB files from the RAP archive or GFS
#'
#' This is a loosely based on `rNOMADS::ArchiveGribGrab`, but customized for our use case.
#' It downloads GRIB/GRIB2 file(s) from NCEI or NCEP, from a particular release date-time
#' (`date_rel`, `hour_rel`) and prediction hour (`hour_pred`), or a sequence of them.
#'
#' The available choices for `model` are:
#'
#' * 'rap_archive' : the archived RAP/RUC files at NCEI
#' * 'gfs_0p25' : the most current GFS release and previous few days
#'
#' Each of `date_rel`, `hour_rel`, and `hour_pred` can be vectors and all times are expected
#' in the UTC time zone. The function matches requests to files based on these three arguments,
#' and `overwrite=FALSE` will skip requests where a matching forecast file was found already
#' on disk.
#'
#' With the 'rap_archive' model there can be multiple matches for a given request: grb2
#' files are always preferred with .grb (version 1) only used as a fall-back. Similarly
#' rap_130 files are preferred over the (older model and coarse resolution) ruc_252 files.
#'
#' Note that with 'rap_archive' the only valid choices for `hour_pred` are 0 and 1, but
#' `hour_rel` can be anything in `seq(24)`. With 'gfs_0p25' the only valid choices for
#' `hour_rel` are 0, 6, 18, but `hour_pred` can be anything in `seq(120)`. To get a list
#' of all choices available on a given date, see `?archive_url`
#'
#' Requests for multiple files are processed in a loops that should continue in spite of
#' download errors. If an error occurs in the middle of a download, the function will delete
#' the (probably corrupt) file and continue on to the next. Multiple requests to NCEI/NCEP
#' servers are spaced out using `Sys.sleep` such that the request rate never exceeds
#' 60 / minute.
#'
#' @param model character, either 'rap_archive' or 'gfs_0p25'
#' @param date_rel Date vector, the forecast release date(s)
#' @param hour_rel integer vector, the forecast release hour
#' @param hour_pred integer vector, the prediction hours
#' @param aoi geometry object defining the area of interest, passed to `archive_url`
#' @param grib_dir character, path to storage directory for GRIB/GRIB2 files
#' @param quiet logical, prints progress information to console
#' @param overwrite logical, whether to overwrite existing files found in `grib_dir`
#'
#' @return data frame with info about the requested times, and results of download attempts
#' @export
archive_get = function(model,
                       date_rel,
                       hour_rel = seq(0, 23, by=2),
                       hour_pred = 1L,
                       aoi = NULL,
                       grib_dir = NULL,
                       quiet = FALSE,
                       overwrite = FALSE)
{
  # time the whole thing
  t_outer = last_ping = proc.time()

  # get a list of all existing files
  grib_df = grib_list(grib_dir, quiet=TRUE)

  # relevant keys in this data frame
  hour_key = c('hour_rel', 'hour_pred')

  # loop over supplied release dates
  n_date = length(date_rel)
  list_out = vector(mode='list', length=n_date)
  for(d in seq(n_date)) {

    msg_count = paste0('(date ', d, '/', n_date,')')
    if( !quiet ) cat(paste0('\n\n> ', format(date_rel[d], '%Y/%m/%d')), msg_count, ' ... ')

    # hours of interest on this date
    date_df = list(hour_rel = hour_rel,
                   hour_pred = hour_pred) |> expand.grid() |> dplyr::tibble()

    # flag existing files on this date by appending a file column
    if( !overwrite ) date_df = grib_df |>
      dplyr::filter(date_rel %in% !!date_rel[d]) |> #!! follows the label outside the data frame
      dplyr::right_join(date_df, by=hour_key)

    # subset of files and hours needed on this date
    request_df = date_df |> dplyr::filter( is.na(file) )
    if( nrow(request_df) == 0 ) {

      # we have all the necessary files already for this date
      cat('completed \U2713')
      next
    }

    # this is to avoid pinging the servers more than once a second
    t_since = (proc.time() - last_ping)['elapsed']
    if(t_since < 1) Sys.sleep(1 - t_since)
    last_ping = proc.time()

    # fetch preferred download URLs
    url_df = archive_url(date_rel[d], model=model, aoi=aoi, quiet=quiet) |>
      dplyr::filter(preferred) |>
      dplyr::mutate(downloaded = FALSE)

    # finished if there are no URLs available
    if( nrow(url_df) == 0 ) next

    # map URLs to requested times
    request_df = request_df |>
      dplyr::select( all_of(hour_key) ) |>
      dplyr::left_join(url_df, by=hour_key)

    # loop over file requests
    for( h in seq( nrow(request_df) ) ) {

      h_rel = request_df[['hour_rel']][h]
      h_pred = request_df[['hour_pred']][h]
      if( !quiet ) cat('\nhour', sprintf('%02d', h_rel), '+', sprintf('%03d', h_pred))

      # skip if we couldn't locate a URL for the file
      if( is.na( request_df[['url']][h] ) ) {

        # no file was found for this release hour
        if( !quiet ) cat(' : not found!')
        next

      } else {

        # prevents pinging the servers more than once a second
        t_since = (proc.time() - last_ping)['elapsed']
        if(t_since < 1) Sys.sleep(1 - t_since)
        last_ping = proc.time()

        # download the file and catch failures
        if( !quiet ) cat(' : file', request_df[['file']][h], 'downloading ... ')
        path_dl = file.path(grib_dir, request_df[['file']][h])
        download_failed = tryCatch({

          request_df[['url']][h] |> download.file(path_dl, mode='wb', quiet=TRUE)

        }, error = function(err) err) |> is('error')

        # report and clean up after failed downloads
        if(download_failed) {

          # file is probably corrupt if it exists
          if( file.exists(path_dl) ) unlink(path_dl)
          if( !quiet ) cat('failed!')
          next

        } else {

          request_df[['downloaded']][h] = TRUE
          cat('done')
        }
      }
    }

    # data frame detailing downloads and failed requests
    list_out[[d]] = request_df
  }

  # collapse list output
  cat('\n\ntotal time elapsed:', round((proc.time() - t_outer)['elapsed'] / 60, 2), 'minutes\n')
  do.call(rbind, list_out) |> dplyr::tibble()
}

