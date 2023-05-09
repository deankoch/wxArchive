library(wxArchive)

# first argument should point to parent directory of "rap" and "gfs"
project_dir = commandArgs(TRUE)[1]
if( length(project_dir) == 0 ) stop('provide the project directory as first argument')

# path to elevation raster in metres
dem_path = project_dir |> file.path('elev_m.tif')

# all the output files go here
base_dir_rap = project_dir |> file.path('rap')

# part 5: fit spatial model to fine grid (don't use resampled layers)
spatial_fit(var_nm = .nm_output_var,
            base_dir = base_dir_rap,
            dem_path = dem_path,
            input_nm = .nm_src_rap[['fine']],
            n_max = 1e3) |> invisible()
