# usage: Rscript list_all [operation] [data_dir]
#
# [data_dir] a directory with sub-folders "rap" and gfs". All weather data is written here
# [operation] one of "list", "update_all", "update_rap", "fit_rap", "impute_rap", "update_gfs", "extract"

# put your own version of this file in data_dir to change the AOI
aoi_nm = 'aoi.geojson'
local_path = file.path('/home/wxarchive', aoi_nm)

library(wxArchive)
operation = commandArgs(trailingOnly=TRUE)[1]
data_dir = commandArgs(trailingOnly=TRUE)[2]
operation_valid =  c('list',
                     'update_all',
                     'update_rap',
                     'fit_rap',
                     'impute_rap',
                     'update_gfs',
                     'export')

# first argument specifies the operation
operation_info = paste(operation_valid, collapse=', ')
if( !(operation %in% operation_valid) ) stop('first argument must be one of: ', operation_info)

# optional second argument should point to parent directory of "rap" and "gfs"
data_dir_exists = data_dir |> normalizePath(winslash='/') |> dir.exists()
if( ( length(data_dir) == 0 ) | !data_dir_exists ) stop('data directory "', data_dir, '" not found')

# copy default AOI polygon file if it's not in the expected location
ext_path = file.path(data_dir, aoi_nm)
if( !file.exists(ext_path) ) {

  message('copying', aoi_nm, 'to', dirname(ext_path))
  if( !file.exists(local_path) ) stop('file ', local_path, ' not found')
  file.copy(local_path, ext_path)
}

# expected location of temporal model directory
model_dir = data_dir |> file.path('rap', wxArchive:::.nm_resample_rap[1], 'model')

# dates to download. NULL sets earliest/latest available
from = NULL
to = NULL

# list stats about recognized files in the project directory
if( operation == 'list' ) data_dir |> wxArchive::workflow_list()

# update RAP/RUC archive
if( operation %in% c('update_rap', 'update_all') ) data_dir |>
  wxArchive::workflow_update_rap(from=from, to=to)

# fit temporal model to RAP/RUC archive
if( operation %in% c('fit_rap', 'update_all') ) {

  # in "update_all" mode we fit the model only if it doesn't exist yet
  if( !dir.exists(model_dir) | (operation == 'fit_rap') ) data_dir |>
    wxArchive::workflow_fit_temporal()
}

# impute missing values in RAP/RUC archive
if( operation %in% c('impute_rap', 'update_all') ) {

  # check if a temporal model has been fitted to the data yet
  if( !dir.exists(model_dir) ) stop('directory ', model_dir, ' not found. Run operation "fit_rap" first')
  data_dir |> wxArchive::workflow_impute_rap()
}

# add wind speed variable
if( operation %in% c('update_all') ) data_dir |> wxArchive::workflow_wnd_rap()

# update GFS archive
if( operation %in% c('update_gfs', 'update_all') ) data_dir |> wxArchive::workflow_update_gfs()

# for each variable this will write a copy of the completed time series in one file
if( operation %in% c('export', 'update_all') ) data_dir |> wxArchive::wxArchive::workflow_export()

# a clue that you can close the bash terminal now
cat('\n')
message('all done')

