
#' Dean Koch
#' April 2023
#'

library(terra)
library(devtools)
load_all()
document()

###
##
# TESTING NEW TIME AGGREGATION UPDATER

data_dir = project_dir = 'G:'
from=NULL
to=NULL
tz='MST'

k = 2
var_nm = .var_daily_pairs[[k]][['var']]
fun = .var_daily_pairs[[k]][['fun']]
base_dir = project_dir
output_nm = .nm_daily
origin_hour = 0L

nc_aggregate_time(base_dir,
                  var_nm,
                  output_nm,
                  fun ,
                  tz,
                  origin_hour)

#
# *
## *   \
### * *  \
#### \ \ / /
# --> --> - -=-o
#### / / \ \
### * *  /
## *   /
# *
#


###
##
#### SPATIAL MODELLING
####

# path to area of interest polygon
data_dir = project_dir = 'G:'
aoi_path = project_dir |> file.path('aoi.geojson')
aoi = sf::st_read(aoi_path)

# all the output files go here
#base_dir_rap = project_dir |> file.path('rap')

var_nm = .var_daily

# testing spatial fit
#var_nm = 'tmp_max'
base_dir = project_dir |> file.path('test')
dem_path = project_dir |> file.path('elev_m.tif')
input_nm = .nm_daily
model_nm = .nm_model
n_max = 5e2
pos = NULL
time_fit = NULL

#spatial_fit(var_nm, base_dir, dem_path, input_nm, model_nm, n_max, pos)

# testing down-scale of entire grid - select subsets with poly_out
poly_path = project_dir |> file.path('export.geojson')
poly_in = poly_path |> sf::st_read()


# TESTING
xx = file.path('D:/UYRW_data', 'data/prepared/usgs_catchments.rds') |> readRDS()
poly_out = xx$boundary

# poly_out should lie entirely within AOI
plot(aoi)
poly_out |> sf::st_geometry() |> sf::st_transform(sf::st_crs(aoi)) |> plot(add=TRUE)

## TODO: put this into a new function `nc_downscale()`

poly_out = poly_out
dem = terra::rast(dem_path)
output_nm = output_nm = .nm_export
var_nm = .var_daily
base_dir
input_nm
model_nm
down = 2 # downscaling factor
edge_buffer = NULL
dates = NULL
write_nc = TRUE

nc_downscale(base_dir,
             dem,
             down = 2,
             input_nm = .nm_daily,
             model_nm = .nm_model,
             output_nm = .nm_export,
             var_nm = .var_daily,
             poly_out = NULL,
             edge_buffer = NULL,
             fun = NULL,
             dates = NULL,
             write_nc = TRUE)




# data input/output paths (var_nm is list to ensure output paths are in list)
input_nc = file_wx('nc', base_dir, input_nm, as.list(var_nm))
var_nm = var_nm |> stats::setNames(nm=names(input_nc))
var_nm_list = names(var_nm) |> as.list()
output_nc = file_wx('nc', base_dir, output_nm, var_nm_list, make_dir=TRUE)

# read fitted parameter values of existing models for the variable(s)
pars_json = file_wx('spatial', base_dir, model_nm, var_nm_list, make_dir=TRUE)
pars_json = pars_json[ sapply(pars_json, file.exists) ]
if( length(pars_json) == 0 ) stop('model file(s) not found in ', model_nm)
pars_all =  pars_json |> lapply(jsonlite::fromJSON)

# get input grid info from first nc file
r_grid_in = nc_chunk(input_nc[[1]][1])[1] |> terra::rast(lyrs=1)
crs_in = terra::crs(r_grid_in)
aoi_in = r_grid_in |> terra::ext() |> sf::st_bbox() |> sf::st_as_sfc(crs=crs_in)
poly_in = poly_in |> sf::st_geometry() |> sf::st_transform(crs=crs_in)

# set default buffer size based on source resolution
if( is.null(edge_buffer) ) edge_buffer = sum(terra::res(r_grid_in)^2) |> sqrt()

# transform output polygons to input projection and find their bounding box
poly_out = poly_in |> sf::st_transform(crs_in)
bbox_out = sf::st_bbox(poly_out) |> sf::st_as_sfc()
bbox_out_big = bbox_out |> sf::st_buffer(dist=edge_buffer)

plot(bbox_out_big)
plot(bbox_out, add=T)
plot(poly_out, add=T)

# set cell values to pixel key
r_grid_in[] = terra::ncell(r_grid_in) |> seq()

# get down-scaled (no data) version of input via snapKrig then crop to bbox
r_grid_out = r_grid_in |>
  snapKrig::sk() |>
  snapKrig::sk_rescale(down=down) |>
  snapKrig::sk_export() |>
  terra::crop(bbox_out_big, snap='out')

# snapKrig version and mapping from uncropped SpatRaster
g_input = r_grid_out |> snapKrig::sk()
is_obs = !is.na(g_input)
idx_obs = g_input[is_obs]

# select a variable
i = 4

# use the latest parameter fit
pars_i = pars_all[[i]][[1]]

# bilinear averaging to get DEM points on same grid as r_grid_out
X_out = r_grid_out |> space_X(dem,
                              dem_knots = pars_i[['knots']],
                              X_center = pars_i[['center']],
                              X_scale = pars_i[['scale']],
                              intercept = FALSE)

g = g_input
pars = pars_i[['pars']]
X = X_out
poly_list = poly_out |> split( seq_along(poly_out) )

# load all data
r_obs = input_nc[[i]] |> nc_layers()

# slow loading of all required data points
z_obs_all = r_obs |> lapply(\(x) x[][idx_obs])
z_obs = z_obs_all[1:10]
#  r_obs can be deleted at this point?


xx = .nc_downscale(g, z_obs, pars, X)
xx |> plot()



split.screen(c(1,2))

screen(1)
r_obs[[2]] |> terra::crop(bbox_out_big) |> plot(reset=F)
plot(bbox_out_big, add=T)
plot(bbox_out, add=T)
plot(poly_out, add=T)

screen(2)
xx |> plot()
plot(bbox_out_big, add=T)
plot(bbox_out, add=T)
plot(poly_out, add=T)





i = 0

i = i + 1
xx[[i]][[1]] |> plot(reset=FALSE)
poly_list[[i]] |> plot(add=TRUE)




# loop over layers
j = 0

j = j + 1


#copy data to snapKrig grid
g_input[['gval']][is_obs] = c(r_obs[[j]][][idx_obs])

# universal kriging
g_output = g_input |> snapKrig::sk_cmean(pars=pars_i[['pars']], X=X_out)
plot(g_output, reset=F)
poly_out |> plot(add=TRUE)






#



poly_out[3,] |> plot(add=TRUE)




# TODO: variance?
#g_var = g_input |> snapKrig::sk_cmean(pars=pars_i[['pars']], X=X_out, what='v')
#plot(g_var)









pts_out = r_grid_out |> snapKrig::sk() |> snapKrig::sk_coords(out='sf')

plot(dem_out, reset=F)
plot(pts_out, add=TRUE)
plot(poly_out, add=TRUE)
plot(bbox_out, add=TRUE)















# bounding box of input raster
poly_bbox


if( is.null(poly_out) ) poly_out = poly_bbox
























start_date = NULL
end_date = NULL


data_dir |> wxArchive::workflow_update_rap(from=start_date, to=end_date)
data_dir |> wxArchive::workflow_impute_rap()
data_dir |> wxArchive::workflow_wnd_rap()
data_dir |> wxArchive::workflow_update_gfs()
data_dir |> wxArchive::workflow_list()


# download RUC/RAP files
# cat('\n')
# message('updating RAP/RUC GRIB archive')
# archive_update(base_dir = base_dir_rap,
#                hour_rel = .hour_rel_rap,
#                model = 'rap_archive') |> invisible()

# r = 'G:/rap/fine/wnd_v.nc/wnd_v_2007.nc'
# xx = time_nc(r)
# yy = time_json(r)
#
# str(xx)
# str(yy)



# # export to nc
# cat('\n')
# message('updating NetCDF files')
# nc_update(aoi = aoi,
#           base_dir = base_dir_rap,
#           output_nm = .nm_src_rap,
#           regex = .rap_regex) |> invisible()


# TODO: fix grib_list posix_rel (converted to integer)



project_dir = 'G:/test'
check_date = function(d) tryCatch(as.Date(d), error = \(e) NULL)
start_date= '2014-09-01' |> check_date()
end_date = '2014-09-01' |> check_date()
data_dir = project_dir

base_dir_rap = project_dir |> file.path('rap')
file_wx('grib', base_dir_rap)

# path to area of interest polygon
aoi_path = project_dir |> file.path('aoi.geojson')
aoi = sf::st_read(aoi_path)

# data_dir |> wxArchive::workflow_list()
data_dir |> wxArchive::workflow_update_rap(from=start_date, to=end_date)

base_dir = base_dir_rap

archive_update(base_dir = base_dir_rap,
               hour_rel = .hour_rel_rap,
               from = start_date,
               to = end_date,
               model = 'rap_archive') |> invisible()

wxArchive::nc_update(aoi = aoi,
                     base_dir = base_dir_rap,
                     output_nm = .nm_src_rap,
                     regex = .rap_regex,
                     from=start_date,
                     to=end_date)

# resampled coarse to fine
var_nm = .nm_output_var
base_dir = base_dir_rap
input_nm = .nm_src_rap
output_nm = .nm_resample
r_fine = NULL
from = NULL






cat('\n')
message('resampling')
nc_resample(var_nm = .nm_output_var,
            base_dir = base_dir_rap,
            input_nm = .nm_src_rap,
            output_nm = .nm_resample) |> invisible()






output_nm = .nm_src_rap
regex = .rap_regex
from=start_date
to=end_date
grib_dir = file_wx('grib', base_dir)
wxArchive::nc_update(aoi = aoi,
                     base_dir = base_dir_rap,
                     output_nm = .nm_src_rap,
                     regex = .rap_regex,
                     from=start_date,
                     to=end_date)

wxArchive::nc_update(aoi = aoi,
                     base_dir = base_dir_rap,
                     output_nm = .nm_src_rap,
                     regex = .rap_regex)



# export to nc
cat('\n')
message('updating NetCDF files')
base_dir_rap = project_dir |> file.path('rap')
base_dir = base_dir_rap
output_nm = .nm_src_rap
regex = .rap_regex
from = NULL
to = NULL
grib_dir = file_wx('grib', base_dir)

nc_update(aoi = aoi,
          base_dir = base_dir_rap,
          output_nm = .nm_src_rap,
          regex = .rap_regex) |> invisible()


# path to area of interest polygon
aoi_path = project_dir |> file.path('aoi.geojson')
aoi = sf::st_read(aoi_path)

p = project_dir |> file.path('rap/fine/tmp.nc/tmp_2022.nc')
x = time_wx(p)

str(x)


r = nc_layers(p)
xx = r |> time_nc()

str(xx)







x = rast(p)
x[][5,] |> plot(type='l')


nc_path = p



project_dir |> workflow_list()

x = .var_daily_pairs[[1]]

base_dir = project_dir
var_nm = x['var']
output_nm = .nm_daily
fun = x['fun']
tz = 'MST'
origin_hour = 0


project_dir |> workflow_update_rap()
project_dir |> workflow_daily()


# copy main nc files to test directory
base_dir_rap = project_dir |> file.path('rap')
out_dir = 'G:/test/rap'

# res_nm = 'fine'
# nc_in = file_wx('nc', base_dir_rap, .nm_src_rap[[res_nm]], as.list(names(.rap_regex)))
# for(i in seq_along(nc_in)) {
# f_out = names(nc_in)[[i]] |> paste0('.nc')
# p_out = out_dir |> file.path(res_nm, f_out)
# nc_layers(nc_in[[i]]) |> nc_write_chunk(p_out)
# }
####
##
#

# GFS

project_dir |> workflow_list()

# path to area of interest polygon
aoi_path = project_dir |> file.path('aoi.geojson')
aoi = sf::st_read(aoi_path, quiet=TRUE)

# all the output files go here
base_dir_gfs = project_dir |> file.path('gfs')

## GFS

# RAP archive end time determines start time of the GFS output
base_dir_rap = project_dir |> file.path('rap')
rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
rap_time = lapply(rap_nc_path, \(p) time_wx(p))
from = do.call(c, lapply(rap_time, \(x) max(x[['time_obs']]) ))

# download GFS files
cat('\n')
message('updating GFS GRIB archive')
gfs_result = archive_update(base_dir = base_dir_gfs,
                            hour_pred = .hour_pred_gfs,
                            hour_rel = .hour_rel_gfs,
                            aoi = aoi,
                            model = 'gfs_0p25',
                            alternate = FALSE)

# delete the old GFS NetCDF directories
base_dir_gfs |> file.path('coarse') |> unlink(recursive=TRUE)
base_dir_gfs |> file.path(.nm_resample) |> unlink(recursive=TRUE)
base_dir_gfs |> file.path(.var_wnd) |> unlink(recursive=TRUE)

# export latest GFS data to nc (creates "coarse" subdirectory)
cat('\n')
nc_update(aoi = aoi,
          base_dir = base_dir_gfs,
          output_nm = list(coarse=.nm_gfs),
          regex = .gfs_regex,
          from = from) |> invisible()

# load an example grid at fine resolution (second rast call drops cell values)
r_fine = nc_chunk(file_wx('nc', base_dir_rap, 'fine', names(.rap_regex)[[1]]))[1] |>
  terra::rast() |> terra::rast()

# resample to match RAP grid
cat('\n')
message('resampling')
nc_resample(var_nm = .nm_gfs_var,
            base_dir = base_dir_gfs,
            input_nm = list(coarse=.nm_gfs),
            output_nm = .nm_resample,
            r_fine = r_fine) |> invisible()

# predictions will run n_ahead days past the end of the GFS data
gfs_nc_path = file_wx('nc', base_dir_gfs, .nm_resample, .nm_gfs_var)
gfs_time = lapply(gfs_nc_path, \(p) time_wx(p))
gfs_end = do.call(c, lapply(gfs_time, \(x) max(x[['time_obs']]) )) |> max()
until = gfs_end + ( n_ahead * 24 * 60 * 60 )

# prediction using model results from RAP analysis
cat('\n')
message('extending forecasts by ', n_ahead, ' day(s)')
time_impute(var_nm = .nm_gfs_var,
            base_dir = base_dir_gfs,
            until = until,
            model_dir = base_dir_rap,
            model_nm = .nm_model,
            input_nm = .nm_resample,
            output_nm = .nm_complete) |> invisible()

# create wind speed layers
cat('\n')
message('computing wind speed from u/v components')
base_dir_gfs |> wnd_update(wnd_nm = .var_wnd,
                           uv_nm = .var_wnd_uv,
                           input_nm = .nm_complete_gfs) |> invisible()





##
###
####

# download RUC/RAP files
cat('\n')
message('updating RAP/RUC GRIB archive')
archive_update(base_dir = base_dir_rap,
               hour_rel = .hour_rel_rap,
               model = 'rap_archive') |> invisible()

# export to nc
cat('\n')
message('updating NetCDF files')
nc_update(aoi = aoi,
          base_dir = out_dir,
          output_nm = list(coarse='coarse', fine='fine'),
          regex = .rap_regex,
          grib_dir = file.path(base_dir_rap, 'grib')) |> invisible()

# copy precip from components (applies to early years)
cat('\n')
message('processing precipitation layers')
pcp_update(base_dir = out_dir,
           pcp_nm = .var_pcp,
           input_nm = .nm_src_rap,
           output_nm = .nm_rap) |> invisible()

# resampled coarse to fine
cat('\n')
message('resampling')
nc_resample(var_nm = .nm_output_var,
            base_dir = out_dir,
            input_nm = .nm_src_rap,
            output_nm = .nm_resample) |> invisible()

# part 6: fit temporal model to fine grid (include all layers)
time_fit(var_nm = .nm_output_var,
         base_dir = out_dir,
         model_nm = 'model',
         input_nm = .nm_resample_rap) |> invisible()

# fill gaps
cat('\n')
message('imputing missing times')
time_impute(var_nm = .nm_output_var,
            base_dir = out_dir,
            model_nm = .nm_model,
            input_nm = .nm_resample_rap,
            output_nm = .nm_complete) |> invisible()

# create/update the wind file
cat('\n')
message('computing wind speed from u/v components')
wnd_update(base_dir = out_dir,
           wnd_nm = .var_wnd,
           uv_nm = .var_wnd_uv,
           input_nm = .nm_complete_rap) |> invisible()










p = out_dir |> file.path('tmp.nc')
tt = time_wx(p)
r = nc_layers(p)


#
p_src =
nc_chunk()
vname = 'tmp'
psrc = 'G:/rap/coarse_lts/tmp.nc'
r = rast(psrc)
p = 'G:/test/tmp.nc'
r |> nc_write_chunk(p)

p = 'G:/test/tmp.nc'
nc_chunk(p)

time_wx(p)
r = nc_layers(p)




p = 'G:/rap/completed'
nc_chunk(p, check_ext=F)
nc_chunk(p)






nc_chunk(p) |> time_wx()



data_dir |> workflow_update_gfs()
project_dir |> workflow_list()

project_dir |> workflow_daily()
project_dir |> workflow_list()



# lake = sf::st_read('D:/rswat_testing/data/nhd_lake.geojson') |> sf::st_geometry()
# river = sf::st_read('D:/rswat_testing/data/nhd_flow.geojson') |> sf::st_geometry()
# catch_df = sf::st_read('D:/rswat_testing/data/nhd_catchment.geojson')
# catch = catch_df |> sf::st_geometry()
#
# i = catch_df$FEATUREID |> order() |> tail(5)
#
# #fname = 'example_areas.png'
# #png(fname, width=25e2, height=45e2, units='px')
#   #plot(catch[i], border=NA)
#   plot(catch, add=F, border='white', col='grey90')
#   plot(river, add=TRUE, col=adjustcolor('blue', alpha.f=0.2))
#   plot(lake, add=TRUE, border=NA, col='blue')
#   plot(catch[i], add=TRUE, border='white', col='red')
# #dev.off()
#
# catch[i] |> sf::st_transform(4326) |> sf::st_write('G:/export.geojson')
#



# create wind speed layers
base_dir_gfs = project_dir |> file.path('gfs')
base_dir_gfs |> file.path(.var_wnd) |> unlink(recursive=TRUE)
cat('\n')
message('computing wind speed from u/v components')
base_dir_gfs |> wnd_update(wnd_nm = 'wnd',
                           uv_nm = c('wnd_u', 'wnd_v'),
                           input_nm = .nm_resample) |> invisible()









# check for spatial model
nm_spatial = 'spatial_daily'
dir_spatial = project_dir |> file.path(nm_spatial)
dir.exists(dir_spatial)


# fit spatial model
var_nm = 'pcp'
base_dir = data_dir
#dem_path =
input_nm = 'fine'
model_nm = input_nm[[1]]
n_max = 5e2
positive = NULL







# starting on down-scaling function
# inputs: from, to,







# # TODO: export the following:
# tmp_max
# tmp_min
# hum_mean
# pcp_mean
# wnd_mean

# p = file.path(project_dir, 'gfs', 'coarse_resampled', 'wnd_v.nc')
# p |> rast() |> crs() |> cat()
#
# library(ncdf4)
# nc = nc_open(p)
# nc |> str()

# DEBUGGING:
library(profvis)

base_dir = data_dir
var_nm = 'pcp'
output_nm = 'export'
write_csv = F
fun='mean'
tz='MST'
origin_hour = 0L

#nc_export(data_dir, var_nm, )

profvis({

  data_dir |> workflow_aggregate(write_csv=FALSE)

  })




# set max RAM available to terra in GB
# terraOptions(memmax=2)
data_dir |> workflow_aggregate()








data_dir |> wxArchive::workflow_list()
data_dir |> workflow_update_rap(from=from, to=to)
data_dir |> workflow_impute_rap()
data_dir |> workflow_wnd_rap()
data_dir |> workflow_update_gfs()
data_dir |> workflow_aggregate()

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
r = nc_aggregate_time(p)

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


