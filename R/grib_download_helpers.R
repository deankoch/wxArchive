#' Get a list of download URLs for GRIB files at NCEI and NCEP
#'
#' This fetches download links for files from the RAP/RUC archive at NCEI,
#' and the (current) GFS forecasts from NCEP. The date arguments specify the
#' release date (when the file was published) and the function returns links for
#' all forecast hours at all available release times.
#'
#' IMPORTANT NOTE: be polite and use the function in moderation! See below:
#'
#' Each function call causes `readLines` to request the HTML content of a webpage from
#' NCEP or NCEI. Users who make requests too frequently can have their IP addresses
#' banned by these servers! There is a 120 requests/minute limit, and the guidance from
#' NCEP is to aim for around 50 requests/minute when making requests in a loop. Enforce
#' this in your own code using `proc.time()` and `Sys.sleep()` as needed.
#'
#' When `model='gfs_0p25'`, the argument `aoi` is used to find appropriate bounding
#' box, and the resulting download links should produce a suitably cropped raster.
#' Only the first 120 prediction times (hours) are returned. GFS forecasts are released
#' 4x daily, so you can expect 1-4 links for each of the 120 prediction hours
#' depending on what time of day you make the call.
#'
#' `aoi` is ignored when `model='rap_archive'`, which doesn't support cropping.
#' All of its files cover all of North America. With this dataset, only the 0
#' and 1-hour prediction times are available, and releases happen hourly, so you can
#' typically expect 2 links for each of 24 release hours.
#'
#' with both datasets, the links are for files containing all variables at all levels.
#' GFS supports requests for a subset of variables, so maybe we can add functionality
#' for that later on. All times are in UTC.
#'
#' @param date_rel Date, the forecast release date of interest
#' @param model either 'rap_archive' or 'gfs_0p25'
#' @param aoi a geometry object, the area of interest (for cropping 'gfs_0p25')
#' @param quiet logical suppresses console output
#'
#' @return a tibble containing URLs and information about the files behind them
#' @export
my_download_url = function(date_rel, model='rap_archive', aoi=NULL, quiet=FALSE) {

  # classes and names in return data frame
  empty_df = data.frame(file = character(0),
                        url = character(0),
                        hour_rel = integer(0),
                        hour_pred = integer(0),
                        preferred = logical(0),
                        posix_rel = character(0) |> as.POSIXct(),
                        posix_pred = character(0) |> as.POSIXct())

  # Note that the most recent 2 years of RAP grids are not found in the archive URL returned by
  # `rNOMADS`. Currently, the cut-off date for archive files is in May 2020, but this is likely
  # to change in the future. If you are having issues with missing files near this date, you
  # should check if `cutoff_date` needs to be updated (using a web browser)

  # parent directory for download URL
  cutoff_date = as.Date('2020-05-15')
  rap_archive_old_URL = 'https://www.ncei.noaa.gov/data/rapid-refresh/access/historical/analysis'
  rap_archive_URL = 'https://www.ncei.noaa.gov/data/rapid-refresh/access/rap-130-13km/analysis'
  gfs_URL = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl'

  # useful strings for naming files and directories
  year = as.integer(format(date_rel, '%Y'))
  month = as.integer(format(date_rel, '%m'))
  day = as.integer(format(date_rel, '%d'))
  year_month_string = paste0(year, sprintf('%02d', month))
  date_string = paste0(year_month_string, sprintf('%02d', day))

  # RAP archive has directory structure for different dates
  if(model == 'rap_archive') {

    # different parent URL depending on requested date
    if( !quiet ) cat('fetching file names from NCEI')
    model_base_URL = if(date_rel < cutoff_date) rap_archive_old_URL else rap_archive_URL
    model_URL = file.path(model_base_URL, year_month_string, date_string)
  }

  # GFS file URLs specify parameters with '?=' syntax
  if(model == 'gfs_0p25') {

    if( !quiet ) cat('fetching file names from NCEP')
    model_URL = gfs_URL |> paste0('?dir=%2Fgfs.', date_string)
  }

  # check for available dates/times
  links = tryCatch({

    # this chunk is copied from rNOMADS::LinkExtractor
    html_text = readLines(model_URL, warn=FALSE)
    links_out = html_text[ which(grepl('href', html_text)) ] |>
      strsplit('</tr><tr><td>') |>
      unlist() |>
      stringr::str_extract('href=\".*\"') |>
      stringr::str_replace('href=', '') |>
      stringr::str_replace_all('\"', '') |>
      stringr::str_replace_all('http:', 'https:') |>
      stringr::str_replace_all(':80', ':443')

    # RAP archive links are copied directly from HTML text
    if( model == 'rap_archive' ) {

      # some extra cleaning of archive links
      char_last = regexpr('>', links_out) - 1L
      extra = char_last > 0
      if( any(extra) ) links_out[extra] = substr(links_out[extra], 1L, char_last[extra])

      # keep only GRIB file links and sort to get 130 < 252, rap < ruc
      links_out = links_out[ grepl('\\.grb2*$', links_out) ] |> sort()

      # parse hours from file names
      basename_len = links_out |> tools::file_path_sans_ext() |> nchar()
      hour_pred = links_out |> substr(basename_len-2L, basename_len) |> as.integer()
      hour_rel = links_out |> substr(basename_len-7L, basename_len-6L) |> as.integer()
      df_out = data.frame(file = links_out,
                          url = file.path(model_URL, links_out),
                          hour_pred = hour_pred,
                          hour_rel = hour_rel,
                          preferred = FALSE)
    }

    # gfs links are constructed based on release hours listed in HTML
    if( model == 'gfs_0p25' ) {

      # get lat/lon bounding box from aoi (or set defaults)
      if( is.null(aoi) ) {

        # this selects everything
        leftlon = 0
        rightlon = 360
        toplat = 90
        bottomlat = -90

      } else {

        # project to WGS84 and map [-180, 180] to [0, 360] longitude
        bbox = aoi |> sf::st_transform('EPSG:4326') |> sf::st_bbox()
        leftlon = floor( bbox['xmin'] ) #%% 360L )
        rightlon = ceiling( bbox['xmax'] ) #%% 360L )
        toplat = ceiling( bbox['ymax'] )
        bottomlat = floor( bbox['ymin'] )
      }

      # copy release hour strings
      n_rel = length(links_out)
      rel_hour_string = links_out |> substr(nchar(links_out)-1L, nchar(links_out))

      # define prediction hours to fetch (first 5 days worth of forecast hours)
      n_pred = 120L
      pred_hour_string = sprintf('%03d', seq(n_pred))

      # define expected file names
      files = paste0('gfs.t',  rel_hour_string, 'z.pgrb2.0p25.f') |>
        lapply(\(link) paste0(link, pred_hour_string) ) |>
        unlist()

      # set up arguments for lat/lon and variable names
      args_subset = paste0('&all_lev=on&all_var=on&subregion=') |>
        paste0('&leftlon=', leftlon) |>
        paste0('&rightlon=', rightlon) |>
        paste0('&toplat=', toplat) |>
        paste0('&bottomlat=', bottomlat)

      # build data frame of URLs, and file names to use locally (similar to RAP/RUC)
      urls = links_out |> rep(each=n_pred) |> paste0('%2Fatmos&file=', files, args_subset)
      df_out = data.frame(file = files, url = urls) |>
        dplyr::mutate( preferred = TRUE ) |>
        dplyr::mutate( hour_pred = rep(seq(n_pred), n_rel) ) |>
        dplyr::mutate( hour_rel = rep(as.integer(rel_hour_string), each=n_pred) ) |>
        dplyr::mutate( file = paste0(model, '_',
                                     date_string, '_',
                                     sprintf('%04d', 1e2 * hour_rel), '_',
                                     sprintf('%03d',  hour_pred), '.grb2') )
    }

    # a data frame with file names and URLs
    df_out

  }, error = function(err) err) |> suppressWarnings()

  # handle failed LinkExtractor calls
  is_error = is(links, 'error')
  if( is_error ) {

    msg_fail = ' \U2713 \nno GRIBs found on the requested date'
    if( !quiet ) cat(msg_fail, format(date_rel, '%Y/%m/%d'))
    return(empty_df)
  }
  cat(' \U2713')

  # add a column for date-time of prediction and publication
  df_out = df_out |>
    dplyr::mutate( posix_rel = as.POSIXct(date_rel, tz='UTC') + ( 60 * 60 * hour_rel ) ) |>
    dplyr::mutate( posix_pred = posix_rel + ( 60 * 60 * hour_pred ) ) |>
    dplyr::arrange( desc(preferred), file ) |>
    dplyr::tibble()

  # set preferred files for rap_archive
  if( !any( df_out[['preferred']] ) ) {

    # files are already ordered by preference, so we just label duplicates
    df_out[['preferred']] = !duplicated(df_out[['posix_pred']])
  }

  return(df_out)
}


#' Initialize or update a local archive of forecasts from RAP/RUC or GFS
#'
#' This creates the subfolder 'gribs' in your `base_dir` if it doesn't exist
#' already, and fills it with archived forecast GRIB files from the selected
#' model.
#'
#' Two models are currently supported: 'rap_archive' downloads the Rapid Refresh
#' (RAP) and Rapid Update Cycle (RUC, the predecessor of RAP) from NCEI;
#' and 'gfs_0p25' fetches the Global Forecast System (GFS) from NCEP. Available
#' release times for the RAP/RUC range from from 2005-01-01 until two days before
#' present; and for GFS there is a rolling window of availability from about 10
#' days before present until the present day.
#'
#' To download all available times from GFS, leave the `from`, `to`, and `hour_rel`
#' arguments `NULL`. Otherwise they specify the start and end dates, and the release
#' hours to look for.
#'
#' `from` and `to` can be NULL in RAP/RUC calls, in which case defaults are assigned
#' based on the date range of existing files in your archive. If there are no existing
#' files, `from` is set to the earliest available date (`.from_def`), otherwise it is
#' set to latest date among the existing files. The default for `to` is latest
#' available date.
#'
#' The time at which a forecast becomes/became valid is equal to the release time
#' `hour_rel` plus the prediction hour `hour_pred`. These can be vectors if you want
#' to select sequences on every date. See `?my_archive_getter` (which does all the work)
#' for details on availability, and how multiple versions of forecasts are handled.
#'
#' @param base_dir path for the GRIB storage subfolder
#' @param from Date, the first date in the sequence
#' @param to Date, the last date in the sequence
#' @param hour_rel integer vector, the release hours
#' @param hour_pred integer vector, the prediction hours
#' @param model character, either 'rap_archive' or 'gfs_0p25'
#'
#' @return a data frame of information about the files downloaded
#' @export
my_update_archive = function(base_dir,
                             from = NULL,
                             to = NULL,
                             hour_rel = 0L,
                             hour_pred = 1L,
                             model = 'rap_archive',
                             aoi = NULL) {

  # get the current date
  date_today = as.Date(Sys.time(), tz='UTC')

  # scan for existing files or create the directory as needed
  grib_dir = wx_file('grib', base_dir, make_dir=TRUE)
  grib_df = grib_dir |> my_archive_lister()
  date_rel = NULL

  # set up times to fetch from RAP
  if(model == 'rap_archive') {

    # fetch only the most recent additions by default
    end_existing = max(grib_df[['date_rel']])
    if( is.null(to) ) to = date_today - 2
    if( is.null(from) ) from = if( nrow(grib_df) == 0 ) .from_def else end_existing
  }

  # set up times to fetch from GFS
  if(model == 'gfs_0p25') {

    # typically the past 10 days are available
    if( is.null(to) ) to = date_today
    if( is.null(from) ) from = date_today - 9L
  }

  # return nothing if date range was invalid
  if( to < from ) return(NULL)
  date_rel = seq.Date(from, to, by='day')

  # download the archive files
  my_archive_getter(model,
                    date_rel = date_rel,
                    hour_rel = hour_rel,
                    hour_pred = hour_pred,
                    aoi = aoi,
                    grib_dir = grib_dir)
}

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
#' of all choices available on a given date, see `?my_download_url`
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
#' @param aoi
#' @param grib_dir character, path to storage directory for GRIB/GRIB2 files
#' @param quiet logical, prints progress information to console
#' @param overwrite logical, whether to overwrite existing files found in `grib_dir`
#'
#' @return data frame with info about the requested times, and results of download attempts
#' @export
my_archive_getter = function(model,
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
  grib_df = my_archive_lister(grib_dir, quiet=TRUE)

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
    url_df = my_download_url(date_rel[d], model=model, aoi=aoi, quiet=quiet) |>
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

