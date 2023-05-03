# earliest available start date for RUC downloads
.from_def = as.Date('2005-01-01')

# table of regex for my variables of interest in RAP/RUC
.rap_regex = c(pcp_total = '^SFC.*hr Total precipitation', # total over prev hour
               pcp_large = '^SFC.*Large scale precipitation', # total attributed to non-convective
               pcp_small = '^SFC.*Convective precipitation', # total attributed to convective
               tmp = '^2\\[m\\].*Temperature',
               hum = '^2\\[m\\].*Relative humidity',
               wnd_u = '^10\\[m\\].*u-component of wind',
               wnd_v = '^10\\[m\\].*v-component of wind')

# same for GFS
.gfs_regex = c(pcp = '^SFC.*hr Total precipitation', # same as RAP but with short name
               .rap_regex[c('tmp', 'hum', 'wnd_u', 'wnd_v')])
