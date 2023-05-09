library(devtools)
load_all()

# get this from user input
project_dir = 'G:'

# base directories for storing all files from RAP/RUC and GFS
base_dir_rap = project_dir |> file.path('rap')
base_dir_gfs = project_dir |> file.path('gfs')

# make a list of all datasets with preference for RAP archive over GFS
rap_nc_path = file_wx('nc', base_dir_rap, .nm_complete_rap, .nm_output_var)
gfs_nc_path = file_wx('nc', base_dir_gfs, .nm_resample, as.list(.nm_gfs_var))
p_all = Map(\(rap, gfs) c(rap, gfs), rap = rap_nc_path, gfs = gfs_nc_path)

#


p_all
