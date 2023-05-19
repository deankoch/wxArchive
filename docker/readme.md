# deankoch/wxarchive docker image
Dean Koch
May 10 2023

### INTRODUCTION

This uses the wxArchive R package to download GRIB files for weather forecasts and convert them to
NetCDF. It is currently set up to download RAP/RUC data as far back as 2005, as well as current
GFS forecasts up to 5 days ahead. These are stitched together and gap-filled to produce a complete
time series.

All output files are written to the directory /home/wxarchive/data in the container. To make this
persistent you will need to map a local volume (some empty directory) to this location. In the
examples below I map the directory G:/ on my local machine using the -v argument.

The docker image contains a default area of interest polygon (the Yellowstone area). Change this
by putting a file named "aoi.geojson" in the directory mapped above (eg in G:/aoi.geojson).

### SET-UP

To build the image you will need "dockerfile", "aoi.geojson", and "Rscript/run_wx.R" in your
dockerfile directory (the same as they are laid out on github). Then, from this directory call

docker build . --no-cache -t deankoch/wxarchive

The --no-cache argument causes the latest version of the `wxArchive` R package to be
downloaded from github if you are building the image repeatedly. 

### OPERATIONS

To run a task in the container, set the WX_OPERATION environmental variable in your docker run call
by appending "-e WX_OPERATION=name", where name is one of:

* "list"           : (the default) lists all available times in the archive at fine resolution
* "update_all"     : runs the seven steps below in sequence

* "update_rap"     : downloads from RAP/RUC and does some transformation
* "fit_rap"        : fits a temporal model to the data at fine resolution
* "impute_rap"     : fills missing time points at fine resolution
* "update_gfs"     : downloads from GFS and does some transformation

* "daily"          : writes a copy of daily aggregate data to NetCDF
* "fit_daily"      : (not yet implemented) fits a spatial model to the daily data at fine resolution
* "export"         : (not yet implemented) downscale and spatially aggregate the daily data


### NOTES

"update_rap" will download data at coarse resolution when it can't find the fine resolution
version for a requested time. After all requested times are downloaded, the function resamples
the coarse data to match the fine resolution grids.

"fit_rap" operates on the fine resolution grids and expects at least one year's worth of data

"fit_rap" must be run before "impute_rap" (or else you must copy model files over manually)

You must download at least one fine resolution grid before running "update_gfs"

"daily" does some further transformation (max, min, mean) to produce daily output, creating
five variables and saving them to your "export" subdirectory

* tmp_max
* tmp_min
* pcp_mean
* hum_mean
* wnd_mean

The last step, "export", summarizes the updated daily times series over the the polygons
found in "export.geojson". Specify your own polygons by modifying this file. There can be any
number of polygons but they must all lie within the bounding box of "aoi.geojson".
Outputs are numbered in the same order as the polygons in "export.geojson".

"fit_daily" must be run at least once (to save a set of model parameter files) before "export".
If a set of parameter files exist already, "update_all" will skip this step and export using
the existing parameters.


### EXAMPLES

docker run --rm -e WX_OPERATION=list -v G:/test/:/home/wxarchive/data deankoch/wxarchive

docker run --rm -e WX_OPERATION=export -v G:/test/:/home/wxarchive/data deankoch/wxarchive

docker run --rm -e WX_OPERATION=update_all -v G:/test/:/home/wxarchive/data deankoch/wxarchive
