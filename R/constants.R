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

# same for GFS (but use short name pcp for pcp_total)
.gfs_regex = c(pcp = .rap_regex[['pcp_total']], .rap_regex[c('tmp', 'hum', 'wnd_u', 'wnd_v')] )

# output variable names from GFS
.nm_gfs_var = names(.gfs_regex)

# names for precipitation (pcp_large + pcp_small = pcp_total = pcp)
.var_pcp = 'pcp'
.var_pcp_old = c('pcp_large', 'pcp_small', 'pcp_total')
.var_not_pcp = names(.rap_regex) |> setdiff(.var_pcp_old)

# list of names to consider equivalent (eg pcp and pcp_total refer to same variable)
.nm_output_var = c(list(c(.var_pcp, 'pcp_total')), as.list(.var_not_pcp))

# list of variables available for export (wind speed is added last)
.var_wnd = 'wnd'
.var_wnd_uv = c('wnd_u', 'wnd_v')
.var_rap_export = .nm_output_var |> c( list(.var_wnd) )
.var_gfs_export = .nm_gfs_var |> c( list(.var_wnd) )

# name and aggregation function for export
.var_daily_pairs = list(c(var='tmp', fun='max'),
                        c(var='tmp', fun='min'),
                        c(var='hum', fun='mean'),
                        c(var='pcp', fun='mean'),
                        c(var='wnd', fun='mean'))


# file names for the aggregate series
.var_daily = .var_daily_pairs |> sapply(\(x) paste(x, collapse='_'))


## DIRECTORIES

# sub-directory names for NetCDF files at two resolutions
.nm_gfs = 'coarse'
.nm_rap = c('coarse', 'fine')
.nm_src_rap = as.list(.nm_rap) |> stats::setNames(.nm_rap)

# sub-directory names for transformed layers
.nm_resample = 'coarse_resampled'
.nm_complete = 'completed'
.nm_daily = 'daily'
.nm_down = 'daily_down'
.nm_export = 'daily_aggregate'

# set of sub-directories to use for fitting temporal model and spatial model
.nm_resample_rap = c(.nm_src_rap[['fine']], .nm_resample)
.nm_model = 'model'

# set of sub-directories forming completed time series
.nm_complete_rap = c(.nm_resample_rap, .nm_complete)
.nm_complete_gfs = c(.nm_resample, .nm_complete)

# set of sub-directories to use for building exports
.nm_rap_export = .nm_complete_rap |> c(.var_wnd)
.nm_gfs_export = .nm_complete_gfs |> c(.var_wnd)

