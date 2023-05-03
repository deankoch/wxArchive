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
#' to select sequences on every date. See `?archive_get` (which does all the work)
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
archive_update = function(base_dir,
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
  grib_df = grib_dir |> grib_list()
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
  archive_get(model,
              date_rel = date_rel,
              hour_rel = hour_rel,
              hour_pred = hour_pred,
              aoi = aoi,
              grib_dir = grib_dir)

}
