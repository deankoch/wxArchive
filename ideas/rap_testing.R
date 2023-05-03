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
base_dir_permanent = 'G:/rap'
base_dir_rap = 'G:/rap_test'
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


grib_dir = my_file_path('grib', base_dir_rap, make_dir=TRUE)

# copy small subset of the archive GRIBS containing both coarse and fine grids and NAs
copy_df = my_file_path('grib', base_dir_permanent) |>
  my_archive_lister(dupe=FALSE) |>
  filter(date_rel < as.Date('2008-12-24')) |>
  filter(date_rel > as.Date('2008-11-30'))
copy_dest_path = grib_dir |> file.path(copy_df[['file']])
file.copy(copy_df[['path']], copy_dest_path, overwrite=FALSE)


# part 2: export all to netCDF long term storage
my_update_nc(aoi = aoi,
             base_dir = base_dir_rap,
             output_nm = nm_src_rap,
             regex = regex_rap) |> invisible()


# add some more GRIBS (with gap) and do an update to establish near term file
copy_df = my_file_path('grib', base_dir_permanent) |>
  my_archive_lister(dupe=FALSE) |>
  filter(date_rel < as.Date('2009-01-02')) |>
  filter(date_rel > as.Date('2008-12-26'))
copy_dest_path = grib_dir |> file.path(copy_df[['file']])
file.copy(copy_df[['path']], copy_dest_path, overwrite=FALSE)

# part 2: export new file to netCDF near term storage
my_update_nc(aoi = aoi,
             base_dir = base_dir_rap,
             output_nm = nm_src_rap,
             regex = regex_rap) |> invisible()

# part 3: compute pcp_total from large + small
my_pcp_total(base_dir = base_dir_rap,
             pcp_nm = var_pcp,
             input_nm = nm_src_rap,
             output_nm = nm_rap) |> invisible()

# part 4: impute fine resolution grids from coarse by spatial resampling
my_resample(var_nm = nm_output_var,
            base_dir = base_dir_rap,
            input_nm = nm_src_rap,
            output_nm = nm_resample) |> invisible()
                    
# part 5: fit spatial model to fine grid (don't use resampled layers)
my_fit_spatial(var_nm = nm_output_var,
               base_dir = base_dir,
               dem_path = dem_path,
               input_nm = nm_src_rap[['fine']],
               n_max = 5e2,
               append = TRUE) |> invisible()

# part 6: fit temporal model to fine grid (include all layers)
my_fit_temporal(var_nm = nm_output_var,
                base_dir = base_dir,
                input_nm = nm_resample_rap) |> invisible()

# part 7: impute missing times in time series
my_impute_temporal(var_nm = nm_output_var,
                   base_dir = base_dir,
                   input_nm = nm_resample_rap,
                   output_nm = 'fine_complete',
                   model_nm = nm_resample_rap[1]) |> invisible()
                



var_nm = nm_output_var
input_nm = c(nm_all_rap[['fine']], rap_resample=resampled_nm)
p_na = NULL
n_max = NULL



base_dir
pcp_nm = 'pcp'
nc_nm = rap_nm

nc_path = input_path$fine$pcp_large
my_nc_attributes(nc_path, overwrite=TRUE, lazy=TRUE, ch=TRUE)

#grib_df = my_file_path('grib', base_dir) |> my_archive_lister()

#
##
###
##
#
aoi
base_dir
nc_nm = lapply(rap_chunk_nm, rev)
regex = regex_rap
n_chunk = 5e3
memory_limit = 8L
dummy_total = TRUE
append = TRUE


# part 2: export to netCDF
my_update_nc(aoi = aoi,
             base_dir = base_dir,
             nc_nm = rap_nm,
             regex = regex_rap) |> invisible()





# part 5: fit spatial model to fine grid (don't use resampled layers)
my_fit_spatial(var_nm = var_nm,
               base_dir = base_dir_rap,
               dem_path = dem_path,
               train_nm = spatial_train_nm,
               append = TRUE) |> invisible()



# testing 
grib_dir = my_file_path('grib', base_dir, make_dir=T)
list_out = my_file_path('nc', base_dir, 'coarse') |> my_nc_attributes()
rap_nm_extra = list(coarse=c('coarse', 'foobar', 'coarse_archive'),
                    fine=c('fine', 'fine_archive'))

var_nm = names(regex_rap) |> as.list()
var_nm[[1]] = c('pcp_total', 'pcp')



my_file_path('nc', base_dir, var_nm=var_nm, sub_dir=rap_nm, collapse=F) 



out_list = my_file_path('nc', base_dir, var_nm=names(regex_rap), sub_dir=rap_nm) 
out_list

data.frame(out_list) |> t() |> as.data.frame() |> as.list()


#out_list |> lapply(\(x) stats::setNames(x, nm=sapply(x, \(y) names(y)[1]) )



my_file_path('nc', base_dir, sub_dir=c('test'), var_nm=var_nm[1:2], collapse=T, simplify=T) |> str()






##








#
##
###
##
#


# copy test file
p_src = my_file_path('nc', base_dir_rap, 'coarse', 'wnd_u')
p_test = my_file_path('nc', base_dir, 'chunks', 'wnd_u', make_dir=T)
if(!file.exists(p_test)) file.copy(p_src, p_test)
p_attr = my_nc_attributes(p_test, overwrite=TRUE)
t_src = p_attr[['time']]

# first chunk is the start
p_big = my_file_path('nc', base_dir, 'chunks', 'wnd_big')
t_big = t_src |> head(-4e3)

# second chunk is the end (much smaller). Include some overlap 
p_small = my_file_path('nc', base_dir, 'chunks', 'wnd_small')
t_small = t_src |> tail(4e3 + 10)

# write both to disk
my_nc_layers(p_test, t_small) |> my_nc_write(p=p_small, overwrite=TRUE, append=FALSE)
my_nc_layers(p_test, t_big) |> my_nc_write(p=p_big, overwrite=TRUE, append=FALSE)

my_nc_attributes(c(p_small, p_big), ch=TRUE)


# select a set of times with duplicate matches and one file is not used
times = c( tail(t_big, 5), head(t_small, 6)) |> unique()
r = my_nc_layers(p_test, times)

# check a random layer
r_test = r[[sample(nlyr(r), 1)]]
r_validate = my_nc_layers(p_src, time(r_test))
max(r_test[] - r_validate[])

# select a set of times with duplicate matches and both files are used
times = c( tail(t_big, 25), head(t_small, 26)) |> unique()
r = my_nc_layers(p_test, times)

# check a random layer
r_test = r[[sample(nlyr(r), 1)]]
r_validate = my_nc_layers(p_test, time(r_test))
max(r_test[] - r_validate[])

my_attributes()



