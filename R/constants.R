## TIME

# earliest available start date for RAP/RUC downloads
.from_def = as.Date('2005-01-01')

# release hours to fetch from each model
.hour_rel_rap = seq(0L, 23L, by=2)
.hour_rel_gfs = c(6L, 18L)

# prediction hours of interest
.hour_pred_rap = 1L
.hour_pred_gfs = seq(1L, 120L, by=2L)


## VARIABLES

# table of regex for my variables of interest in RAP/RUC
.rap_regex = c(pcp_total = '^SFC.*hr Total precipitation', # total over prev hour
               pcp_large = '^SFC.*Large scale precipitation', # total attributed to non-convective
               pcp_small = '^SFC.*Convective precipitation', # total attributed to convective
               tmp = '^2\\[m\\].*Temperature',
               hum = '^2\\[m\\].*Relative humidity',
               wnd_u = '^10\\[m\\].*u-component of wind',
               wnd_v = '^10\\[m\\].*v-component of wind')

# same for GFS (only pcp is different)
.gfs_regex = c(pcp = '^SFC.*hr Total precipitation', # same as RAP but with short name
               tmp = '^2\\[m\\].*Temperature',
               hum = '^2\\[m\\].*Relative humidity',
               wnd_u = '^10\\[m\\].*u-component of wind',
               wnd_v = '^10\\[m\\].*v-component of wind')

# output variable names from GFS
.nm_gfs_var = names(.gfs_regex)

# output names for precipitation variables split from the others
.var_pcp = 'pcp'
.var_pcp_old = c('pcp_large', 'pcp_small', 'pcp_total')
.var_not_pcp = names(.rap_regex)[ !( names(.rap_regex) %in% .var_pcp_old ) ]

# list of names to consider equivalent
.nm_output_var = c(list(c(.var_pcp, 'pcp_total')), as.list(.var_not_pcp))


## DIRECTORIES

# sub-directory names for NetCDF files at two resolutions
.nm_gfs = 'coarse'
.nm_rap = c('coarse', 'fine')
.nm_old_rap = .nm_rap |> paste0('_lts')

# list with old and new at both resolutions
.nm_src_rap = cbind(.nm_rap, .nm_old_rap) |>
  apply(1, identity, simplify=FALSE) |>
  stats::setNames(.nm_rap)

# sub-directory names for transformed layers
.nm_resample = 'coarse_resampled'
.nm_complete = 'completed'

# set of sub-directories to use for fitting temporal model and spatial model
.nm_resample_rap = c(.nm_src_rap[['fine']], .nm_resample)
.nm_spatial = 'fine'

# set of sub-directories forming completed time series
.nm_complete_rap = c(.nm_resample_rap, .nm_complete)

