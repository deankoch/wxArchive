
#' Dean Koch
#' April 2023
#'

library(terra)
library(devtools)
load_all()
document()


project_dir = 'G:'
from = NULL
to = NULL
data_dir = project_dir

# # TODO: export the following:
# tmp_max
# tmp_min
# hum_mean
# pcp_mean
# wnd_mean



data_dir |> workflow_export()




data_dir |> wxArchive::workflow_list()
data_dir |> workflow_update_rap(from=from, to=to)
data_dir |> workflow_impute_rap()
data_dir |> workflow_wnd_rap()
data_dir |> workflow_update_gfs()

base_dir = data_dir
var_nm = 'hum'
output_nm = 'export'
write_csv = TRUE
fun = 'mean'
tz = 'MST'
origin_hour = 0L


###
project_dir |> nc_list()

project_dir |> workflow_list()
project_dir |> workflow_wnd_rap()


file.path(project_dir, 'rap') |> wnd_update()

###
base_dir = file.path(project_dir, 'rap')
wnd_nm = 'wnd'
uv_nm = c('wnd_u', 'wnd_v')
input_nm = .nm_complete_rap
output_nm = wnd_nm





input_path = file_wx('nc', base_dir, input_nm, var_nm=as.list(uv_nm))
p = input_path[[1]]
vinfo = time_wx(p)
vinfo[['time_obs']]
vinfo[['time_na']]


r = p |> nc_layers(vinfo[['time_obs']], na_rm=TRUE)
x = r[]
xx = x |> apply(2, anyNA)
xx |> which()



times = time(r)[xx]

p
times





###

#base_dir |> nc_export()



##





##

#project_dir |> workflow_list()
p_all = project_dir |> workflow_list(quiet=TRUE)
p = p_all[[2]]
r = nc_aggregate(p)

r[[1]] |> plot()





project_dir |> workflow_update_rap()
project_dir |> workflow_impute_rap()
project_dir |> workflow_update_gfs()


message('\nupdating NetCDF files')
aoi_path = project_dir |> file.path('aoi.geojson')
aoi = sf::st_read(aoi_path)
base_dir_rap = project_dir |> file.path('rap')
nc_update(aoi = aoi,
          base_dir = base_dir_rap,
          output_nm = .nm_src_rap,
          regex = .rap_regex) |> invisible()

# copy precip from components (applies to early years)
message('\nprocessing precipitation layers')
pcp_update(base_dir = base_dir_rap,
           pcp_nm = .var_pcp,
           input_nm = .nm_src_rap,
           output_nm = .nm_rap) |> invisible()

# resampled coarse to fine
message('\nresampling')
nc_resample(var_nm = .nm_output_var,
            base_dir = base_dir_rap,
            input_nm = .nm_src_rap,
            output_nm = .nm_resample) |> invisible()

hour_rel
hour_pred


xx = c(NA, 2, NA, NA)
for(i in 1:4) {

  if( is.na(xx[i]) ) next
  cat('test i =', i)


}




# project directories, an area of interest in North America and a fine-resolution DEM covering it
base_dir_rap = 'G:/rap'
base_dir_gfs = 'G:/gfs'

dem_path = 'G:/weather_db/NAM_elv_msk.tif'
aoi_path = 'G:/weather_db/aoi.geojson'

# build an AOI polygon
aoi = sf::st_read(aoi_path)

# release hours to fetch from each model
hour_rel_rap = seq(0, 23, by=2)
hour_rel_gfs = c(6L, 18L)

# prediction hours of interest
hour_pred_rap = 1L
hour_pred_gfs = seq(1, 120, by=2)

# sub-directories for different resolutions and time periods of RAP/RUC
nm_rap = c('coarse', 'fine')
nm_old_rap = nm_rap |> paste0('_lts')

# sub-directories for synthesized/transformed layers
nm_resample = 'coarse_resampled'
nm_complete = 'completed'

#  sets of subdirectories with semi-completed series
nm_src_rap = cbind(nm_rap, nm_old_rap) |> apply(1, identity, simplify=FALSE) |> stats::setNames(nm_rap)
nm_resample_rap = c(nm_src_rap[['fine']], nm_resample)
nm_complete_rap = c(nm_resample_rap, nm_complete)

# which to use for model fitting
nm_temporal = nm_resample
nm_spatial = 'fine'

# variables to fetch in GRIBs from each model
regex_rap = .rap_regex
regex_gfs = .gfs_regex
var_pcp = 'pcp'
var_pcp_old = c('pcp_large', 'pcp_small', 'pcp_total')
var_other = names(.rap_regex)[ !( names(.rap_regex) %in% var_pcp_old ) ]

# list grouping variable names that are considered equivalent in each step
nm_gfs_var = names(regex_gfs)
nm_output_var = c(list(c(var_pcp, 'pcp_total')), as.list(var_other))

#
##
###
##
#

# library(terra)
#
rap_nc_path = file_wx('nc', base_dir_rap, nm_complete_rap, nm_output_var)
gfs_nc_path = file_wx('nc', base_dir_rap, nm_resample, as.list(nm_gfs_var))

nc_path = rap_nc_path[['wnd_u']][['fine']]
r = terra::rast(nc_path)
rt = terra::time(r)
length(rt)
length(unique(rt))

xx = time_wx(nc_path)
rtx = xx$time
length(rtx)
length(unique(rtx))
all(rt == rtx)



# max(rt)
#
#
# nc_path = rap_nc_path[['tmp']]['fine']
# r = rast(nc_path)
# rt = time(r)
# length(rt)
# length(unique(rt))
# min(rt)
#
#
#
#
# time_wx(rap_nc_path[['tmp']])[['time']] |> range()
#


#
##
###
##
#

base_dir = base_dir_rap
output_nm = nm_src_rap
regex = regex_rap
n_chunk = 5e3
memory_limit = 8L
make_dummy = FALSE
from = NULL
append = TRUE


# # part 2: export all to netCDF long term storage
nc_update(aoi = aoi,
          base_dir = base_dir_rap,
          output_nm = nm_src_rap,
          regex = regex_rap) |> invisible()


#
##
###
##
#


base_dir = base_dir_rap
pcp_nm = var_pcp
input_nm = nm_src_rap
output_nm = nm_rap


# part 3: compute pcp_total from large + small
pcp_update(base_dir = base_dir_rap,
           pcp_nm = var_pcp,
           input_nm = nm_src_rap,
           output_nm = nm_rap) |> invisible()


#
##
###
##
#

var_nm = nm_output_var
base_dir = base_dir_rap
input_nm = nm_src_rap
output_nm = nm_resample
r_fine = NULL
from = NULL


# part 4: impute fine resolution grids from coarse by spatial resampling
nc_resample(var_nm = nm_output_var,
            base_dir = base_dir_rap,
            input_nm = nm_src_rap,
            output_nm = nm_resample) |> invisible()


#
##
###
##
#


var_nm = nm_output_var
base_dir = base_dir_rap
input_nm = nm_resample_rap
output_nm = nm_complete
model_nm = input_nm[1]
n_max = NULL
until = NULL
quiet = FALSE

# part 8: impute missing times in fine grid series (run time_fit first)
time_impute(var_nm = nm_output_var,
            base_dir = base_dir_rap,
            input_nm = nm_resample_rap,
            output_nm = nm_complete) |> invisible()




nc_path_rap = file_wx('nc', base_dir_rap, 'coarse', 'pcp')
time_wx(nc_path_rap)[['time_obs']] |> range()




nc_path_rap = file_wx('nc', base_dir_rap, nm_resample_rap, nm_output_var)
nc_path = nc_path_rap[[3]]
nc_path

times = time_wx(nc_path)


