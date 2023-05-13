# deankoch/wxarchive docker image
Dean Koch
May 10 2023

### INTRODUCTION

This uses the wxArchive R package to download GRIB files for weather forecasts and convert them to
NetCDF. It is currently set up to download RAP data from the past week but it can be configured to
go as far back as 2005.

All output files are written to the directory /home/wxarchive/data in the container. To make this
persistent you will need to map a local volume to this location. In the examples below I map the
empty directory G:/test/ on my local machine to /home/wxarchive/data with the -v argument.

The image contains a default area of interest polygon (the Yellowstone area). Change this by putting
a file named "aoi.geojson" in the directory mapped above (eg in G:/test/aoi.geojson).


### OPERATIONS

Set the "WX_OPERATION" environmental variable to select a task:

* "list"           : (the default) lists all available times in the archive at fine resolution
* "update_all"     : runs the five steps below in sequence

* "update_rap"     : downloads from RAP/RUC and does some transformation
* "fit_rap"        : fits a temporal model to the data at fine resolution
* "impute_rap"     : fills missing time points at fine resolution
* "update_gfs"     : downloads from GFS and does some transformation
* "export"         : writes a copy of daily aggregate data to NetCDF and CSV


### NOTES

"update_rap" will download data at coarse resolution when it can't find the fine resolution
version for a requested time. After all requested times are downloaded, the function resamples
the coarse data to match the fine resolution grids.

"fit_rap" operates on the fine resolution grids and expects at least one year's worth of data

"fit_rap" must be run before "impute_rap" (or else you must copy model files over manually)

You must download at least one fine resolution grid before running "update_gfs"

"export" does some further transformation (max, min, mean) to produce daily output, creating
five variables and saving them (as .csv and .nc) to your "export" subdirectory:

* tmp_daily_max
* tmp_daily_min
* pcp_daily_mean
* hum_daily_mean
* wnd_daily_mean


### EXAMPLES

docker run --rm -e WX_OPERATION=list -v G:/test/:/home/wxarchive/data deankoch/wxarchive

docker run --rm -e WX_OPERATION=export -v G:/test/:/home/wxarchive/data deankoch/wxarchive

docker run --rm -e WX_OPERATION=update_all -v G:/test/:/home/wxarchive/data deankoch/wxarchive
