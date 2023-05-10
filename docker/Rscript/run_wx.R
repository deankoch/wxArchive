library(wxArchive)

# usage: Rscript list_all [command] [data_dir]
command = commandArgs(trailingOnly=TRUE)[1]
data_dir = commandArgs(trailingOnly=TRUE)[2]

# first argument specifies the operation
command_valid =  c("list", "update_all", "update_rap", "impute_rap", "update_gfs", "extract")
command_info = paste(command_valid, collapse=', ')
if( !(command %in% command_valid) ) stop('first argument must be one of: ', command_info)

# optional second argument should point to parent directory of "rap" and "gfs"
if( length(data_dir) == 0 ) data_dir = 'home/wxarchive/data'
if( !dir.exists(data_dir) ) stop('data directory ', data_dir, ' not found')

# list stats about recognized files in the project directory
if( command == 'list' ) data_dir |> wxArchive::workflow_list()

# update RAP/RUC archive (temporarily set this to )
if( command %in% c('update_rap', 'update_all') ) data_dir |> workflow_update_rap(to=as.Date('2005-02-01'))

# impute missing values in RAP/RUC archive
if( command %in% c('impute_rap', 'update_all') ) data_dir |> workflow_impute_rap()

# update GFS archive
if( command %in% c('update_gfs', 'update_all') ) data_dir |> workflow_update_gfs()

# update GFS archive
if( command %in% c('extract', 'update_all') ) stop('not yet implemented')
