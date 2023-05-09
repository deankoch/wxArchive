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
archive_url = function(date_rel, model='rap_archive', aoi=NULL, quiet=FALSE) {

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

    msg_fail = ' \nno GRIBs found on the requested date'
    if( !quiet ) cat(msg_fail, format(date_rel, '%Y/%m/%d'))
    return(empty_df)
  }

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
