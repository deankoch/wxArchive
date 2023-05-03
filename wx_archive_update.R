
#' Dean Koch
#' April 2023
#'

library(devtools)
load_all()
document()

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
my_update_nc(aoi = aoi,
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
my_update_nc(aoi = aoi,
             base_dir = base_dir_gfs,
             output_nm = list(coarse='coarse'),
             regex = regex_gfs,
             from = t_gfs,
             append = is.null(t_gfs)) |> invisible()


# resample done fresh each update? Test this

# find latest observed time with all variables present in RAP/RUC archive
nc_path_rap = wx_file('nc', base_dir_rap, nm_resample_rap, nm_output_var)
nc_tmax = do.call(c, lapply(nc_path_rap, \(x) max( my_nc_attributes(x, ch=TRUE)[['time_obs']] ) ))

# load an example grid at fine resolution (second rast call drops cell values)
r_fine = wx_file('nc', base_dir_rap, nm_spatial, names(regex_rap)[[1]]) |>
  terra::rast() |> terra::rast()

# part 9: resample to match fine
my_resample(var_nm = nm_gfs_var,
            base_dir = base_dir_gfs,
            input_nm = list(coarse='coarse'),
            output_nm = nm_resample,
            r_fine = r_fine,
            from = nc_tmax,
            append = FALSE) |> invisible()
