# Rscript
# Dean Koch, May 2023
#
# USAGE: Rscript run_wx.R [data_dir] [operation] [start_date] [end_date]
#
# [data_dir] directory containing sub-folders "rap" and gfs". All outputs go here
# [operation] usually either "list" or "update_all". See below for list of all valid operation names
# [start_date] start of the date range to process (set to non-Date to detect automatically)
# [end_date] end of the date range to process (set to non-Date to detect automatically)
#
# * fix posix_rel in grib_df, which is getting turned into integer somewhere (then rebuild CSVs)

library(wxArchive)

# put your own version of this file in data_dir to change the AOI
aoi_nm = 'aoi.geojson'
local_path = file.path('/home/wxarchive', aoi_nm)

# expected model files (if missing they are created automatically)
from_temporal_file = 'temporal_model.zip'
from_spatial_file = 'spatial_model.zip'

# valid operations
operation_valid =  c('list',
                     'update_all',
                     'update_rap',
                     'fit_temporal',
                     'impute_rap',
                     'update_gfs',
                     'daily',
                     'fit_spatial',
                     'downscale',
                     'extract')

# get user input
msg_args = 'Usage: Rscript run_wx [data_dir] [operation] [start_date] [end_date]'
user_args = commandArgs(trailingOnly=TRUE)
if( length(user_args) != 4 ) stop(msg_args)

# first argument should point to parent directory of "rap" and "gfs"
data_dir = user_args[1]
data_dir_exists = data_dir |> normalizePath(winslash='/') |> dir.exists()
if( ( length(data_dir) == 0 ) | !data_dir_exists ) stop('data directory "', data_dir, '" not found')

# second argument specifies the operation
operation = user_args[2]
operation_info = paste(operation_valid, collapse=', ')
if( !(operation %in% operation_valid) ) stop('first argument must be one of: ', operation_info)

# (optional) last two arguments specify the target date range for updates or output
check_date = function(d) tryCatch(as.Date(d), error = \(e) NULL)
start_date = user_args[3] |> check_date()
end_date = user_args[4] |> check_date()

# copy default AOI polygon file if it's not in the expected location
ext_path = file.path(data_dir, aoi_nm)
if( !file.exists(ext_path) ) {

  message('copying', aoi_nm, 'to', dirname(ext_path))
  if( !file.exists(local_path) ) stop('file ', local_path, ' not found')
  file.copy(local_path, ext_path)
}

# list stats about recognized files in the project directory
if( operation == 'list' ) data_dir |> wxArchive::workflow_list()

# update RAP/RUC archive
if( operation %in% c('update_rap', 'update_all') ) data_dir |>
  wxArchive::workflow_update_rap(from=start_date, to=end_date)

# fit temporal model to RAP/RUC archive
if( operation %in% c('fit_temporal', 'update_all') ) {

  # this unpacks a saved set of model parameter files unless in fit_temporal mode
  if(operation == 'fit_temporal') from_temporal_file = NULL
  data_dir |> wxArchive::workflow_fit_temporal(from_file=from_temporal_file)
}

# impute missing values in RAP/RUC archive
if( operation %in% c('impute_rap', 'update_all') ) {

  # expected location of temporal model directory
  model_dir = data_dir |> file.path(wxArchive:::.nm_temporal_model)

  # check if a temporal model has been fitted to the data yet
  msg_fail = paste('directory', model_dir, 'not found. Run operation "fit_temporal" first')
  if( !dir.exists(model_dir) ) stop(msg_fail)
  data_dir |> wxArchive::workflow_impute_rap()
}

# add wind speed variable
if( operation %in% c('update_all') ) data_dir |> wxArchive::workflow_wnd_rap()

# update GFS archive
if( operation %in% c('update_gfs', 'update_all') ) data_dir |> wxArchive::workflow_update_gfs()

# for each variable this writes a copy of the completed daily time series
if( operation %in% c('daily', 'update_all') ) data_dir |> wxArchive::workflow_daily()

# fit a spatial model to the daily data
if( operation %in% c('fit_spatial', 'update_all') ) {

  # this unpacks a saved set of model parameter files unless in fit_temporal mode
  if(operation == 'fit_spatial') from_spatial_file = NULL
  data_dir |> wxArchive::workflow_fit_spatial(from_file=from_spatial_file)
}

# downscale and aggregate (makes completed daily time series at two scales)
if( operation %in% c('downscale', 'update_all') ) data_dir |> wxArchive::workflow_downscale()
if( operation %in% c('extract', 'update_all') ) data_dir |> wxArchive::workflow_extract()


cat('\n')
message('all done')

