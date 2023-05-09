library(sf)
library(devtools)
load_all()

# get this from user input
project_dir = 'G:'
project_dir |> workflow_update_rap()
project_dir |> workflow_impute_rap()
project_dir |> workflow_update_gfs()


