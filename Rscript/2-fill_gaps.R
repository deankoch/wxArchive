library(wxArchive)

# first argument should point to parent directory of "rap" and "gfs"
project_dir = commandArgs(TRUE)[1]
if( length(project_dir) == 0 ) stop('provide the project directory as first argument')
project_dir |> workflow_impute_rap()


