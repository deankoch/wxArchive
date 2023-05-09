
#' Dean Koch
#' April 2023
#'
#'
#'

library(devtools)
#install_github('deankoch/snapKrig')
load_all()
document()

project_dir = 'G:'
grib_dir = file_wx('grib', file.path(project_dir, 'rap'))
grib_df = grib_dir |> grib_list()
project_dir |> workflow_list()

# is_zero = grib_df$hour_pred == 0
# sum(is_zero)
# sum( grib_df[['date_rel']][is_zero] > as.Date('2006-10-26') )

grib_df[['date_rel']][is_zero] |> hist()



project_dir |> workflow_update_rap()
# project_dir |> workflow_impute_rap()
# project_dir |> workflow_update_gfs()

# set me to TRUE on first run. This adds about 1-2 hours processing time
model_fitting_run = FALSE

# project directories, an area of interest in North America and a fine-resolution DEM covering it
# base_dir_rap = 'G:/weather_db/rap'
# base_dir_gfs = 'G:/weather_db/gfs'
base_dir_rap = 'G:/rap'
base_dir_gfs = 'G:/gfs'

dem_path = 'G:/weather_db/NAM_elv_msk.tif'
aoi_path = 'G:/weather_db/aoi.geojson'

# build an AOI polygon
aoi = sf::st_read(aoi_path)

#
##
###
##
#

# part 1: download GRIBs
archive_update(base_dir = base_dir_rap,
               hour_rel = .hour_rel_rap,
               model = 'rap_archive') |> invisible()

# To speed things on the initial run, set `dummy_total=TRUE` to skip
# loading NA pcp_total layers in early years

# part 2: export all to netCDF long term storage
nc_update(aoi = aoi,
          base_dir = base_dir_rap,
          output_nm = .nm_src_rap,
          regex = .rap_regex) |> invisible()


# part 3: compute pcp_total from large + small
pcp_update(base_dir = base_dir_rap,
             pcp_nm = .var_pcp,
             input_nm = .nm_src_rap,
             output_nm = .nm_rap) |> invisible()

# To speed things up after the initial run, copy all the nc files and the
# time subfolder from "coarse" to "coarse_lts" (and same for fine).
# Subsequent updates will be saved to a smaller new file in the "coarse"
# subdirectory (making updates much faster), and the older data will
# still be available for all downstream steps.

# part 4: impute fine resolution grids from coarse by spatial resampling
nc_resample(var_nm = .nm_output_var,
            base_dir = base_dir_rap,
            input_nm = .nm_src_rap,
            output_nm = .nm_resample) |> invisible()

# These steps only have to be run once initially to establish a model
# for imputation and interpolation. Run them again occasionally to update
# the fitted parameters.

if(model_fitting_run) {

# run these in any order

  var_nm = .nm_output_var
  base_dir = base_dir_rap
  dem_path = dem_path
  input_nm = .nm_src_rap[['fine']]
  n_max = 1e3
  model_nm = input_nm[[1]]
  positive = NULL


# part 5: fit spatial model to fine grid (don't use resampled layers)
spatial_fit(var_nm = .nm_output_var,
            base_dir = base_dir_rap,
            dem_path = dem_path,
            input_nm = .nm_src_rap[['fine']],
            n_max = 1e3) |> invisible()

# part 6: fit temporal model to fine grid (include all layers)
time_fit(var_nm = .nm_output_var,
         base_dir = base_dir_rap,
         input_nm = .nm_resample_rap) |> invisible()

}


# part 8: impute missing times in fine grid series (run time_fit first)
time_impute(var_nm = .nm_output_var,
            base_dir = base_dir_rap,
            input_nm = .nm_resample_rap,
            output_nm = .nm_complete) |> invisible()


# after updating the temporal model, delete the contents of `.nm_complete` and
# re-run part 8 to get new imputed values for missing times in archive

#
##
###
##
#
# GFS

# part 7: grab new GFS forecast data (up to 10 days worth of releases)
gfs_result = archive_update(base_dir = base_dir_gfs,
                            hour_pred = .hour_pred_gfs,
                            hour_rel = .hour_rel_gfs,
                            aoi = aoi,
                            model = 'gfs_0p25')


# find the earliest of the forecast times downloaded above
# time_added_gfs = gfs_result |> dplyr::filter(downloaded) |> dplyr::pull(posix_pred)
# t_gfs = if( length(time_added_gfs) == 0 ) NULL else min(time_added_gfs)

# time info from the completed RAP/RUC archive
rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
rap_time = lapply(rap_nc_path, \(p) time_wx(p))

# find the latest observed times in the RAP/RUC archive
from = do.call(c, lapply(rap_time, \(x) max(x[['time_obs']]) ))

# delete the old GFS NetCDF directories
base_dir_gfs |> file.path('coarse') |> unlink(recursive=TRUE)
base_dir_gfs |> file.path(.nm_resample) |> unlink(recursive=TRUE)

# part 8: export latest GFS data to nc (creates "coarse" subdirectory)
nc_update(aoi = aoi,
          base_dir = base_dir_gfs,
          output_nm = list(coarse='coarse'),
          regex = .gfs_regex,
          from = from) |> invisible()

# load an example grid at fine resolution (second rast call drops cell values)
r_fine = file_wx('nc', base_dir_rap, .nm_spatial, names(.rap_regex)[[1]]) |>
  terra::rast() |> terra::rast()

# part 9: resample to match fine
nc_resample(var_nm = .nm_gfs_var,
            base_dir = base_dir_gfs,
            input_nm = list(coarse='coarse'),
            output_nm = .nm_resample,
            r_fine = r_fine) |> invisible()



#
##
###
##
# TODO: put this stuff into helper functions
#
# nc_get()
# nc_daily()
#

# merge datasets and prefer RAP archive over GFS
rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
gfs_nc_path = file_wx('nc', base_dir_gfs, .nm_resample, as.list(.nm_gfs_var))
p_all = Map(\(rap, gfs) c(rap, gfs), rap = rap_nc_path, gfs = gfs_nc_path)

# load example variable
var_i = 2
.nm_gfs_var[var_i] |> print()
p = p_all[[var_i]]
p_attr = time_wx(p)
range(p_attr[['time_obs']])

# load all layers into RAM (~5GB)
r = nc_layers(p, times=p_attr[['time_obs']], na_rm=TRUE)
t_obs = terra::time(r)

# xx = r |> terra::global('mean') |> as.matrix() |> as.numeric()
# plot(xx ~ t_obs, type='l')

# check for gaps
ts_df = data.frame(posix_pred=t_obs) |> archive_pad()


# step size in the data
step_data = time_step(ts_df)

# desired starting hour and time zone for alignment
origin_hour = 0L
tz = 'MST'

# starting hour of the data in target time zone
start_hour = as.integer( format(min(t_obs), '%H', tz=tz) )
start_idx = ( (origin_hour - start_hour) %% 24 ) / step_data
start_date = t_obs[start_idx] |> as.Date(tz=tz)

# list of indices to aggregate in each step and corresponding dates
n_per = 24L / step_data
n_out = floor( ( length(t_obs) - start_idx + 1 ) / n_per )
list_idx = seq(n_out) |> lapply(\(i) seq(n_per) + (i-1)*n_per )
date_out = seq.Date(start_date, by='day', length.out=n_out)

# aggregate
library(terra)
r_min = do.call(c, lapply(list_idx, \(j) min(r[[j]])) )
r_mean = do.call(c, lapply(list_idx, \(j) mean(r[[j]])) )
r_max = do.call(c, lapply(list_idx, \(j) max(r[[j]])) )
#terra::time(r_mean) = date_out

# pick a grid point at random
k_test = sample(terra::ncell(r_mean), 1)

# extract time series
y_low = r_min[][k_test,] |> c()
y = r_mean[][k_test,] |> c()
y_high = r_max[][k_test,] |> c()

all_len = length(y)
show_len = 5e2
idx_plot = all_len + seq(show_len) - 1

idx_plot = idx_plot - show_len
plot(y[idx_plot]~date_out[idx_plot], type='l', ylim=range(c(y_low, y_high)))
lines(y_low[idx_plot]~date_out[idx_plot], col='grey')
lines(y_high[idx_plot]~date_out[idx_plot], col='grey')



#
##
###
##
#

# merge datasets and prefer RAP archive over GFS
rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
gfs_nc_path = file_wx('nc', base_dir_gfs, .nm_resample, as.list(.nm_gfs_var))
p_all = Map(\(rap, gfs) c(rap, gfs), rap = rap_nc_path, gfs = gfs_nc_path)

# load example variable
var_i = 2
.nm_gfs_var[var_i] |> print()
p = p_all[[var_i]]
p_attr = time_wx(p)
range(p_attr[['time_obs']])

# load all layers into RAM (~5GB)
r = nc_layers(p, times=p_attr[['time_obs']], na_rm=TRUE)
t_obs = terra::time(r)

# xx = r |> terra::global('mean') |> as.matrix() |> as.numeric()
# plot(xx ~ t_obs, type='l')

# check for gaps
ts_df = data.frame(posix_pred=t_obs) |> archive_pad()

# identify times that were resampled from RAP (violet)
p_res = file_wx('nc', base_dir_rap, .nm_resample, names(p_all)[[var_i]])
t_res = time_wx(p_res)[['time_obs']]

# identify times that were imputed (red)
p_impute = file_wx('nc', base_dir_rap, .nm_complete, names(p_all)[[var_i]])
t_imp = time_wx(p_impute)[['time_obs']]

# identify times from GFS (blue)
p_gfs = file_wx('nc', base_dir_gfs, .nm_resample, names(p_all)[[var_i]])
t_gfs = time_wx(p_gfs)[['time_obs']]

# pick a random grid point and plot its time series
k_test = sample(terra::ncell(r), 1)
y_black = y_red = y_blue = y_violet = r[][k_test,] |> c()
y_black[t_obs %in% c(t_gfs, t_imp, t_res)] = NA
#y_red[t_obs %in% c(t_gfs, t_res)] = NA
y_blue[t_obs %in% c(t_imp, t_res)] = NA
y_violet[t_obs %in% c(t_gfs, t_imp)] = NA

all_len = length(y_black)
show_len = 1e3
idx_plot = all_len + seq(show_len) - 1

idx_plot = idx_plot - show_len
plot(y_red[idx_plot] ~ t_obs[idx_plot], type='l', col='red')
lines(y_blue[idx_plot] ~ t_obs[idx_plot], col='blue')
lines(y_violet[idx_plot] ~ t_obs[idx_plot], col='violet')
lines(y_black[idx_plot] ~ t_obs[idx_plot])




#
# t_obs |> length()
# t_obs |> unique() |> length()
#
# # consistency check
# all(p_attr[['time_obs']] == t_obs)
#
# #p_attr[['time_obs']][p_attr[['time_obs']] - t_obs]
# ts_df |> dplyr::filter(is.na(ts_hours))




