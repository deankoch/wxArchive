#' rap_workflow.R
#' Dean Koch
#' March 2023
#' 
#' Create and maintain a 2-hour archive of NOAA's Rapid Refresh (RAP) model
#' and its predecessor, Rapid Update Cycle (RUC), from 2005-present
#' 
#' Instructions: run this entire script daily to update the weather database
#' for SWAT forecast simulations.

library(rNOMADS)
library(terra)
library(sf)
library(forecast)
library(dplyr)
library(jsonlite) 
library(snapKrig)

source('D:/rapid_refresh/helpers_misc.R')

# set me to TRUE on first run. This adds about 1-2 hours processing time
model_fitting_run = FALSE 

# project directory
base_dir_rap = 'G:/rap'
base_dir_gfs = 'G:/gfs'

# an area of interest in North America and a fine-resolution DEM covering it
dem_path = 'L:/spatial_data/NAM_elv_msk.tif'
aoi_path = 'D:/UYRW_data/data/prepared/uyrw_nhd_boundary.rds'

# build an AOI polygon
aoi = readRDS(aoi_path) |> st_buffer(units::set_units(250, km))

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

# part 1: download GRIBs
archive_update(base_dir = base_dir_rap,
                  hour_rel = hour_rel_rap,
                  model = 'rap_archive') |> invisible()

# To speed things on the initial run, set `dummy_total=TRUE` to skip
# loading NA pcp_total layers in early years

# part 2: export all to netCDF long term storage
nc_update(aoi = aoi,
             base_dir = base_dir_rap,
             output_nm = nm_src_rap,
             regex = regex_rap) |> invisible()

# part 3: compute pcp_total from large + small
my_pcp_total(base_dir = base_dir_rap,
             pcp_nm = var_pcp,
             input_nm = nm_src_rap,
             output_nm = nm_rap) |> invisible()

# To speed things up after the initial run, copy all the nc files and the
# time subfolder from "coarse" to "coarse_lts" (and same for fine).
# Subsequent updates will be saved to a smaller new file in the "coarse"
# subdirectory (making updates much faster), and the older data will
# still be available for all downstream steps.

# part 4: impute fine resolution grids from coarse by spatial resampling
my_resample(var_nm = nm_output_var,
            base_dir = base_dir_rap,
            input_nm = nm_src_rap,
            output_nm = nm_resample) |> invisible()

# These steps only have to be run once initially to establish a model
# for imputation and interpolation. Run them again occasionally to update
# the fitted parameters.

if(model_fitting_run) {
  
# run these in any order

# part 5: fit spatial model to fine grid (don't use resampled layers)
my_fit_spatial(var_nm = nm_output_var,
               base_dir = base_dir_rap,
               dem_path = dem_path,
               input_nm = nm_src_rap[['fine']],
               n_max = 1e3,
               append = TRUE) |> invisible()

# part 6: fit temporal model to fine grid (include all layers)
time_fit(var_nm = nm_output_var,
                base_dir = base_dir_rap,
                input_nm = nm_resample_rap,
                n_max = NA) |> invisible()

}

# part 8: impute missing times in fine grid series (run time_fit first)
time_impute(var_nm = nm_output_var,
                   base_dir = base_dir_rap,
                   input_nm = nm_resample_rap,
                   output_nm = nm_complete) |> invisible()


# after updating the temporal model, delete the contents of `nm_complete` and
# re-run part 8 to get new imputed values for missing times in archive

#
##
### 
##
#
# GFS

# part 7: grab new GFS forecast data (up to 10 days worth of releases)
gfs_result = archive_update(base_dir = base_dir_gfs,
                               hour_pred = hour_pred_gfs,
                               hour_rel = hour_rel_gfs,
                               aoi = aoi,
                               model = 'gfs_0p25')

# find the earliest of the forecast times downloaded above
time_added_gfs = gfs_result |> dplyr::filter(downloaded) |> dplyr::pull(posix_pred)
t_gfs = if( length(time_added_gfs) == 0 ) NULL else min(time_added_gfs) 

# part 8: export latest GFS data to nc
nc_update(aoi = aoi,
             base_dir = base_dir_gfs,
             output_nm = list(coarse='coarse'),
             regex = regex_gfs,
             from = t_gfs,
             append = is.null(t_gfs)) |> invisible()


# resample done fresh each update? Test this

# find latest observed time with all variables present in RAP/RUC archive
nc_path_rap = file_wx('nc', base_dir_rap, nm_resample_rap, nm_output_var)
nc_tmax = do.call(c, lapply(nc_path_rap, \(x) max( my_nc_attributes(x, ch=TRUE)[['time_obs']] ) ))

# load an example grid at fine resolution (second rast call drops cell values)
r_fine = file_wx('nc', base_dir_rap, nm_spatial, names(regex_rap)[[1]]) |> 
  terra::rast() |> terra::rast()

# part 9: resample to match fine
my_resample(var_nm = nm_gfs_var,
            base_dir = base_dir_gfs,
            input_nm = list(coarse='coarse'),
            output_nm = nm_resample,
            r_fine = r_fine,
            from = nc_tmax,
            append = FALSE) |> invisible()

#
##
### 
##
#
# testing

# merge datasets and prefer RAP archive over GFS
p_all = Map(\(rap, gfs) c(rap, gfs),
            gfs = file_wx('nc', base_dir_gfs, nm_resample, as.list(nm_gfs_var)),
            rap = file_wx('nc', base_dir_rap, nm_complete_rap, nm_output_var))
             
# load example variable
var_i = 2
nm_gfs_var[var_i] |> print()
p = p_all[[var_i]]
p_attr = my_nc_attributes(p, ch=TRUE)

# load all layers into RAM (~5GB)
r = nc_layers(p, times=p_attr[['time_obs']], na_rm=TRUE)
t_obs = time(r)

# consistency check
all(p_attr[['time_obs']] == t_obs)

# check for gaps
ts_df = data.frame(posix_pred=t_obs) |> archive_pad()

# identify times that were imputed
p_impute = file_wx('nc', base_dir_rap, nm_complete, names(p_all)[[var_i]])
t_imp = my_nc_attributes(p_impute, ch=TRUE)[['time_obs']]

# identify times from GFS
p_gfs = file_wx('nc', base_dir_gfs, nm_resample, names(p_all)[[var_i]])
t_gfs = my_nc_attributes(p_gfs)[['time_obs']]

# pick a random grid point and plot its time series
k_test = sample(ncell(r), 1)
y_black = y_red = y_blue = r[][k_test,] |> c()
y_black[t_obs %in% c(t_gfs, t_imp)] = NA
y_red[t_obs %in% t_gfs] = NA

all_len = length(y_black)
show_len = 2e3
idx_plot = all_len + seq(show_len) - 1

idx_plot = idx_plot - show_len
plot(y_blue[idx_plot] ~ t_obs[idx_plot], type='l', col='blue')
lines(y_red[idx_plot] ~ t_obs[idx_plot], col='red')
lines(y_black[idx_plot] ~ t_obs[idx_plot])


lines(y_lm[idx_plot] ~ t_obs[idx_plot], col=adjustcolor('green', alpha.f = 0.5))

# add linear model lines
pp = 'G:/rap/fine/model/temporal/tmp_2023-04-22.20.49.18.nc'
ppp = 'G:/rap/fine/model/tmp_temporal.json'
model_info = jsonlite::fromJSON(readLines(ppp))[[1]]
pars_r = terra::rast(pp)
betas = pars_r[][,-seq(2)][k_test,]
alphas = pars_r[][,seq(2)][k_test,]
knots = model_info[['knots']]
#model_info
# betas = rast(pp)[][,-seq(2)][k_test,]
# alphas = rast(pp)[][,seq(2)][k_test,]
X = time_X(t_obs, knots_t=knots)
y_lm = X %*% betas




idx = idx_plot[seq(1e2) + 3300]
plot(y_blue[idx] ~ t_obs[idx], type='l', col='blue')
lines(y_red[idx] ~ t_obs[idx], col='red')
lines(y_black[idx] ~ t_obs[idx])
lines(y_lm[idx] ~ t_obs[idx], col=adjustcolor('green', alpha.f = 0.5))


z_res = y_black[idx] - y_lm[idx]


# recurse forward to compute expected values, filling residuals matrix
for( j in which(is.na(z_res))-2 ) {
  
  # sum of lag 1 and 2 terms 
  z_res[j+2] = ( alphas[1] * z_res[j+1] ) + ( alphas[2] * z_res[j] )
}

z_out = z_res + y_lm[idx]
plot(z_out~t_obs[idx], type='l')
lines(y_lm[idx] ~ t_obs[idx], col=adjustcolor('green', alpha.f = 0.5))
#lines(y_red[idx] ~ t_obs[idx], col='red')









idx = -1
idx = idx + 1
idx_plot = seq(show_len) + idx * show_len






# 
# ttest = '2023-04-18 19:00:00' |> as.POSIXct(tz='UTC')
# r = nc_layers(p[1], times=ttest)
# plot(r)
# ttest
# time(r)
# 
# 
# is_na = r[] |> apply(2, anyNA)
# t_obs[is_na]

# tt = as.POSIXct('2007-08-14 21:00:00', tz='UTC')
# tt %in% my_nc_attributes(p[5])$time_obs
# r2 = nc_layers(p[5], tt, na_rm=TRUE)
# r2

# # check for gaps
# ts_df = data.frame(posix_pred=t_obs) |> archive_pad()
# ts_df[is.na(ts_df$ts_hours),]

#g = r |> sk()








show_len = 2e3
# idx = ceiling(length(y_plot)/show_len)
# 
# idx = idx - 1
# idx_plot = seq(show_len) + idx * show_len 
idx = -1

idx = idx + 1
idx_plot = seq(show_len) + idx * show_len

# compute linear predictor

plot(y_blue[idx_plot] ~ t_obs[idx_plot], type='l', col='blue')
lines(y_red[idx_plot] ~ t_obs[idx_plot], col='red')
lines(y_black[idx_plot] ~ t_obs[idx_plot])
lines(y_lm[idx_plot] ~ t_obs[idx_plot], col='green')

# impute missing values?

my_nc_attributes(nm_all)
nc_layers()


# part 9/10: aggregate to daily



#' 
#' # DETAILS
#' 
#' #' ## IMPORTANT FUNCTIONS
#' 
#' 1. `archive_update` initializes the archive and downloads new files #
#' 2. `nc_update` exports the (many) GRIB files to (a few) monolithic nc files 
#' 3. `my_pcp_total` fills gaps in "pcp_total" using component precipitation
#' 3. `my_fine_from_coarse` fills gaps in fine resolution series using coarse grid data
#' 3. `my_fit_spatial` fits a model of spatial covariance to assist with gap-filling
#' 4. `time_fit` fits an AR(2) model to residuals to assist with gap-filling
#' 5. `my_fill_gaps` fills missing time points with imputed data
#' 6. `my_aggregate` produces daily aggregate variables suitable for input to SWAT+
#' 
#' ## LAYER NAMES
#' 
#' The GRIB files we are download here have named layers identifying the variables.
#' Some of these names have changed slightly (at least once) in the history of RAP/RUC.
#' To get around this we use regular expressions (regex) instead of exact matching.
#' 
#' Select a particular set of variables by modifying the regular expressions in the
#' `regex = .rap_regex` at the beginning of the workflow. This can be changed after
#' step 1 (all variables are included in all downloads), but not after step 2, which
#' extracts the requested subset of variables and discards the rest.
#' 
#' `names(regex)` gives a short form name for each variable to use as a file name.
#' The default in `.rap_regex` selects the following named variables
#' 
#' 1. "pcp_total" = 01 hr Total precipitation (kg/m2)
#' 2. "pcp_small" = Large scale precipitation at ground or water surface (kg/m2)
#' 3. "pcp_large" = Convective precipitation at ground or water surface (kg/m2)
#' 4. "tmp" = Temperature at 2m above ground (C)
#' 5. "hum" = Relative humidity at 2m above ground (%)
#' 6. "wnd_u" = wind u-component (roughly West->East) at 10m above ground (m/s)
#' 7. "wnd_v" = wind v-component (roughly South->North) at 10m above ground (m/s)
#' 
#' There are often multiple options for height above ground, so make sure you
#' are selecting the one you want with your `regex`. For example the regex for
#' precipitation total is "^SFC.*hr Total precipitation", which selects names that
#' begin with "SFC" (short for "surface", or 0 metres vertical height).
#' 
#' 
#' ## COMPUTATION TIME
#' 
#' After Part 1, the workflow script spends the vast majority of its compute time either
#' reading layer data from disk, or writing it to the nc file. This is very slow, and I
#' think a faster physical storage device would help a lot (I'm using am external disk via USB).
#' 
#' There are some compression settings to tweak here as well but I'm not sure if they're
#' accessible from `terra`, or if their potential impact would be worth the fuss. 
#' 
#' 
#' ## PART 1 `archive_update` (populates "/grib")
#' 
#' This part downloads archived forecast grids from NCEI servers, using `rNOMADS` to
#' get the https download links. These servers are sometimes down for maintenance. If
#' you encounter download errors due to network issues, just halt the function and call
#' it again when the issues are resolved. The script will skip downloads of forecast
#' times found among the existing files in your `rap_base_dir`.
#' 
#' Forecast files are posted to the online archive on a delay of about two days, so the
#' default call fetched all times up to the present (`Sys.Date() |> as.Date()`) minus two
#' days. Once some of the downloads are completed, you can run rap_gaps.R to make a chart
#' showing availability of the two resolutions in your local archive, by time and date.
#' 
#' Note that the default request is for over 70,000 files (at the time of writing) - it
#' will likely take more than a week to download everything, so I recommend setting the
#' arguments `from` and `to` in `archive_update` to select smaller chunks at a time for
#' downloading. Once the bulk of it is downloaded, smaller updates (dozens of files) should
#' complete in under a minute.
#' 
#' 
#' ## PART 2 `nc_update` (creates or updates "/fine" and "/coarse")
#' 
#' We are interested in a small number of variables (about 8 of the 300+ in each
#' GRIB), and a relatively small sub-grid (the Yellowstone/Tetons and surrounding
#' 250km bubble) in each of the GRIB files. Part 2 repackages this data into netCDF
#' files of size 50-1000 MB, each containing the entire observed time series for one
#' of the variables. These nc files are must faster to load and update (seconds to
#' minutes, instead of half a day).
#' 
#' Once the nc files are initialized, small updates of a few days worth of forecasts
#' should complete in around 10 minutes.
#' 
#' Initializing the nc files for the whole archive takes a long time - about 12 hours. 
#' To make the script more robust to interruptions, `n_chunk` sets a maximum number
#' of GRIB files to load before writing results to disk (default 5000). Previously
#' loaded files are always skipped, so you can interrupt the "loading" phase at any time
#' and restart the script at a later time to resume from the current chunk.
#' 
#' The output from this script is two sub-folders in your `rap_base_dir`, "coarse" and
#' "fine", for the nc files at two resolutions. Each of these directories has a
#' sub-folder called "time" containing JSONs with names matching the nc files.
#' The JSON files are for keeping track of NAs and timestamps.
#' 
#' 
#' ## PART 3 `my_pcp_total` (creates/updates "pcp.nc" in "/fine" and "/coarse")
#' 
#' This part defines a new variable "pcp" and creates/updates the "pcp.nc" file on
#' disk at each resolution. "pcp" is just a copy of "pcp_total" where missing time
#' layers are imputed using the sum of "pcp_large" and "pcp_small". I have found that
#' this closely approximates "pcp_total", and it is available for long time period
#' early on (pre-2007) where "pcp_total" would otherwise be missing
#' 
#' This part is relatively fast, taking about 5 minutes to complete. Like in the
#' last part, existing times are skipped when updating existing files, and a matching
#' JSON is created/updated in the 'time' subfolders.
#' 
#' 
#' ### PART 4 `my_fine_from_coarse` (creates/updates "/fine_from_coarse")
#' 
#' This creates a copy of the "/fine" files (excluding all precipitation variables
#' except "pcp"), then fills temporal gaps by resampling grids from "/coarse". It should
#' complete in a about 5-10 minutes.
#' 
#' In the early years of the archive - and in numerous small gaps later on - forecasts
#' are available at coarse resolution only. This part converts those coarse grids
#' to fine resolution using bilinear averaging. Resampling is handled by GDAL so it can be
#' done very quickly on the whole stack.
#' 
#' The results are merged with the existing fine resolution data and saved to like named
#' files in "/fine_from_coarse", along with the associated JSON files in the "time"
#' sub-folder. Layers with (any number of) negatives are set to NA here.
#' 
#' 
#' ### PART 5 `my_fit_spatial` (creates or updates "/fine/model")
#' 
#' It will be useful later on to have a simple model for spatial trends and covariance.
#' This part fits such a model to a random subsample of layers from each variable, saving
#' model parameters to a JSON file in the "model" subdirectory.
#' 
#' We do this for the fine resolution grids only, and use only the original data (not the
#' resampled layers). Layer on we will use these model parameters to downscale grids from
#' the "fine" resolution here to a much finer one. 
#' 
#' Model fitting is done with snapKrig, which can be very fast when using a relatively small
#' number of subsample layers (`n_max` in the hundreds). Model parameters, along with an
#' index of the training set, are saved in a list in the JSON output file. The `append=TRUE`
#' argument causes new results to be appended to any existing list, so we can execute this
#' part repeatedly to grow a list of model fits.
#' 
#' ...
#' 
#' 
#' ## REFERENCES
#' 
#' https://nomads.ncep.noaa.gov/
#' https://www.r-bloggers.com/2015/04/how-to-install-rnomads-with-grib-file-support-on-windows/
