#' rap_show_gaps.R
#' Dean Koch
#' March 2023
#' 
#' Plots a chart showing presence/absence of forecast files at two resolutions in the
#' 2-hourly RAP archive. This should be run after downloading the files with rap_get.R
#' 
#' My notes on RAP/RUC (see also https://rapidrefresh.noaa.gov/ and https://ruc.noaa.gov/ruc/)
#' 
#' > RUC becomes available 2005-01-01 at 25km resolution
#' > additional 13km grids offered from 2007-04-01 onward (but missing for entire year of 2008)
#' > largest gaps in RUC happen in 2008, where the entire months of January, March, April are missing
#' > RUC replaced by RAP on 2012-05-01 (with no significant gap in the transition)
#' > largest gap in RAP is entire month of January 2014
#' > several other large gaps in 13km series can be filled by 25km files
#' > overall about 5% of times are missing (at both resolutions)
#' 
#' I think this is all much easier to see with a chart. This script makes the chart by scanning
#' for downloaded GRIB files, and parsing dates, times, and other information from their names. 

# this graphic written to disk by the script
png_path = 'D:/rapid_refresh/gaps.png'

# helper function for archive data
source('D:/rapid_refresh/helpers_misc.R')

# create a data frame of GRIB files and missing times, split by resolution
base_dir = 'L:/spatial_data/rap'
grib_df = base_dir |> 
  my_file_path(what='grib') |> 
  my_archive_lister() |> 
  my_archive_padder() |> 
  dplyr::mutate( gap_days = gap/24 ) |>
  dplyr::mutate( hour_of_day = as.numeric(format(posix_pred, '%H', tz='UTC')))

# identify resolutions in use
label_res = grib_df |> split(~coarse) |>
  stats::setNames(c('fine', 'coarse')) |> 
  lapply(\(x) as.numeric(x[1L, c('nx', 'ny')])) |>
  sapply(\(x) paste0('(', paste(x, collapse=' x '), ')'))

# set up colors for background and files from the two resolutions
my_cols = list(bg = adjustcolor('grey70', alpha.f=0.2),
               fine = adjustcolor('orange', alpha.f=0.8),
               coarse = adjustcolor('violet', alpha.f=0.8))

# set up a custom x axis with monthly ticks
from = grib_df[['posix_pred']] |> min() 
to = grib_df[['posix_pred']] |> max() 
xtick = seq(from, to, by='month')
xtick_label = format(xtick, '%b', tz='UTC')
xtick_label[xtick_label == 'Jan'] = format(xtick[xtick_label == 'Jan'], '%Y') |> paste('Jan')

# write an image to disk as .png file
png_height_px = 3200
png(png_path, width=10*png_height_px, height=png_height_px, units='px', pointsize=108)
{
  # set up 2-panel stacked plot
  split.screen(c(2,1))
  
  # upper panel (shares x axis with lower)
  screen(1)
  par(mar=c(0.1, 4.1, 4.1, 2.1))
  
  # initialize plot for gap length
  title_extra = paste0('(', as.Date(from, tz='UTC'), ' to ', as.Date(to, tz='UTC'), ' UTC)')
  plot(gap_days ~ posix_pred,
       data = grib_df,
       pch = NA,
       xaxt = 'n',
       xlab = NA,
       ylab = 'gap length (days)',
       main = paste('RAP/RUC archive gaps and coverage', title_extra),
       frame = FALSE)
  
  # add neutral background color, month and hour lines, then gap lines
  usr = par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col=my_cols[['bg']], border=NA)
  abline(v=xtick, col='white')
  abline(h=seq(round(usr[4L])), col='white')
  lines(gap_days ~ posix_pred, data=grib_df)
  
  # legend of colors for high and low res results in lower plot
  leg_label = names(label_res) |> sapply(\(nm) paste(nm, label_res[[nm]]))
  legend('topright',
         legend = leg_label,
         fill = adjustcolor(my_cols[names(label_res)], alpha.f=0.5),
         title = 'resolution')
  
  # lower panel shows hours observed coloured by resolution
  screen(2)
  par(mar=c(6.1, 4.1, 0.1, 2.1))
  
  # plot observed hours
  plot(hour_of_day ~ posix_pred,
       data = grib_df,
       pch = NA,
       ylab = 'prediction hour',
       xlab = NA,
       xaxt = 'n',
       frame = FALSE)
  
  # add neutral background color, month lines, custom x axis
  usr = par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col=my_cols[['bg']], border=NA)
  abline(v=xtick, col='white')
  axis(side=1, at=xtick, labels=xtick_label, las=2)
  
  # add vertical lines to indicate files, look over resolutions
  for( nm_res in names(label_res) ) {
    
    # filter drops the NA rows so we add them back again below
    df_plot = grib_df |> dplyr::filter(coarse == (nm_res == 'coarse'))
    if(nrow(df_plot) > 0) lines(hour_of_day ~ posix_pred,
                                data =  my_archive_padder(df_plot, quiet=TRUE),
                                col = my_cols[[nm_res]] ) 
  }

  close.screen(all.screens=TRUE)
}
dev.off()
