# wxArchive: a database of weather forecasts from RAP/RUC/GFS 
Dean Koch
**updated May 19 2023**

## Overview

This R package collects 2-hourly data from the following NOAA forecast models: 

* Rapid Update Cycle (RUC) from 2005 to 2007 (26km resolution)
* Rapid Refresh (RAP) from 2007 to present less 2 days (13km resolution)
* Global Forecast System (GFS) from 10 days before present to 5 days ahead (28km resolution)

The package handles downloads, model-fitting, and prediction automatically and produces
outputs in convenient formats like netCDF, CSV and geoJSON. It can be scheduled to run daily
in a docker container to produce an up-to-date, long-running, gap-less times series of your
variables of interest.

## Background

The RAP/RUC GRIB source files cover all of North America, so the package crops them to an 
an area of interest (AOI) supplied by the user before saving a copy. GFS forecasts can be
requested for specific bounding boxes, so the GFS GRIB source files cover only the AOI.

Source files are first saved to disk in their original GRIB/GRIB2 format before being
converted to NetCDF v4 for easier access. Each source file contains about 310 variables
(eg temperature at various altitudes) but only a small number of relevant ones are copied
to NetCDF, and only the sub-grid covering the AOI is included. 

Data points from RAP/RUC are comprised of 1-hour-ahead forecasts for odd-numbered hours,
whenever they are available. When these are unavailable, the corresponding 0-hour forecast
are used instead. Data points from GFS comprise the n-hour-ahead ahead forecasts for all
odd-numbered n from 1 to 119 (inclusive). These are drawn from the daily releases at
06:00 and 18:00, with a preference for the latest release. 

These various data sources are joined together and missing times are imputed automatically.
The result is a 2-hourly gap-less time series extending from 2005 until several days into
the future (see below). Optionally the package will aggregate to daily variables, and
efficiently extrapolate (using `snapKrig`) and/or spatially aggregate outputs over
user-supplied polygons

## Variables

For our purposes we need the following variables

* `pcp_total` = `'^SFC.*hr Total precipitation'`
* `pcp_large` = `'^SFC.*Large scale precipitation'`
* `pcp_small` = `'^SFC.*Convective precipitation')`
* `tmp` = `'\^2\\[m\\].*Temperature')`
* `hum` = `'\^2\\[m\\].*Relative humidity')`
* `wnd_u` = `'\^10\\[m\\].*u-component of wind')`
* `wnd_v` = `'\^10\\[m\\].*v-component of wind')`

The first string (eg `pcp_total`) is my name for the variable, whereas the one in parentheses is a regular expression to match layer names in the GRIB/GRIB2 files. The wildcards accommodate a small degree of variation in names among the three models, and among the different forecast times published on a given day. Notice that a vertical level has to be specified (eg "SFC" stands for "surface") and square brackets in the name should be escaped (eg `\\[` instead of `[`).

The precipitation type we are interested in is `pcp_total`. Variables `pcp_large` and `pcp_small` are included so that we can derive an estimate of `pcp_total` for time layers where it is missing.


## Grids

All outputs use coordinates from the LCC projection, matching the 13km RAP grid. The projection WKT is as follows:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJCRS["unnamed",
    BASEGEOGCRS["Coordinate System imported from GRIB file",
        DATUM["unnamed",
            ELLIPSOID["Sphere",6371200,0,
                LENGTHUNIT["metre",1,
                    ID["EPSG",9001]]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433,
                ID["EPSG",9122]]]],
    CONVERSION["Lambert Conic Conformal (2SP)",
        METHOD["Lambert Conic Conformal (2SP)",
            ID["EPSG",9802]],
        PARAMETER["Latitude of false origin",25,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8821]],
        PARAMETER["Longitude of false origin",-95,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8822]],
        PARAMETER["Latitude of 1st standard parallel",25,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8823]],
        PARAMETER["Latitude of 2nd standard parallel",25,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8824]],
        PARAMETER["Easting at false origin",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8826]],
        PARAMETER["Northing at false origin",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8827]]],
    CS[Cartesian,2],
        AXIS["easting",east,
            ORDER[1],
            LENGTHUNIT["metre",1,
                ID["EPSG",9001]]],
        AXIS["northing",north,
            ORDER[2],
            LENGTHUNIT["metre",1,
                ID["EPSG",9001]]]]

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To harmonize the spatial layout of the three grids, the 26km and 28km datasets are automatically resampled to 13km resolution using bilinear averaging with GDAL.


## Missing Data

About 10% of the RUC/RAP time points are missing (as of 2023-05-02). This can happen for a number of different reasons: a forecast might never get been published due to operational issues; or a file might simply not have been made available on the NCEI/NCEP archive server when it was last checked. Network issues can also lead to failed downloads or corrupted GRIB files, and this is not unlikely given the large number of file transfers (around 80,000).

Time points missing for any of these reasons are automatically imputed by one of three methods, attempted in this order:

1. derive the missing variable by combining other variables (eg precip from its components)
2. downscale the RUC layer (if present) to replace a missing RAP layer
3. Use a fitted AR2 model with seasonal trends to estimate the remaining missing points

## Updates

The update script can be run as often as you like. Running it once daily should be sufficient. It should complete a one-day update in well under 30 minutes. An update does the following:

* check existing and available RAP/RUC files and download any (applicable) new ones
* clean the data and and append to existing NetCDF files
* imputation: compute total precipitation as sum of component precipitation (as needed)
* imputation: resample to 13km (as needed)
* imputation: predict missing times using AR2 and seasonal trend (as needed)

The first two steps are then repeated for GFS. Then the same gap-filling routine is used to impute any missing/future unseen times.


## Files

The total data volume is about 1.3 TB (as of 2023-05-02). Almost all of this space is used up by source files in GRIB2 format. The output files in NetCDF format currently make up around 6 GB or about 0.5% of the total. Each daily update introduces about 250 MB of new source files.

Each NetCDF file contains data for a single variable, but can include numerous time points. To help speed up access in files with many times, every NetCDF file has an associated JSON file that indexes times and missing data points. This
JSON can be found in the "time" subfolder of the directory containing the netCDF file (eg. if the file is "/\<foo\>.nc" then the JSON is "/time/\<foo\>.json").

The NetCDF data for a given variable \<var\> is spread out over several files, using the following directory structure:

* rap/coarse/\<var\>.nc stores layers from RUC/RAP (26km)
* rap/fine/\<var\>.nc stores layers from RAP (13km)
* rap/coarse_resampled/\<var\>.nc stores layers from RUC resampled to match RAP
* rap/completed/\<var\>.nc stores layers imputed using temporal model trained on RAP
* rap/wnd/wnd.nc stores wind speed layer computed from wnd_u and wnd_v

Each of these items is a directory containing NetCDF files split by year. For example "rap/coarse/tmp.nc" is a directory
containing files of the form "rap/coarse/tmp.nc/tmp_\<year\>.nc"

The GFS model only appears in one resolution (which I am naming "coarse"), so it has a slightly simpler structure:

* gfs/coarse/\<var\>.nc stores layers at their original resolution
* gfs/coarse_resampled/\<var\>.nc stores a version resampled to match RAP
* gfs/completed/\<var\>.nc stores layers imputed using RAP temporal model
* gfs/wnd/wnd.nc stores wind speed layer computed from wnd_u and wnd_v

Three additional top-level folders store output files

* daily/<out_var>.nc stores daily aggregate values; like mean, max, min over time
* spatial/<out_var>.nc stores spatial aggregate values of daily; like mean, max, min over space
* model stores model parameter files (either fitted to the data or loaded from a previous fit)


## Reading/Exporting

*this needs updating*

Our data access functions (in R) are designed to operate on the the whole grid at a specified subset of times, with support for chunking time across different files. This means you can have different sets of times in different files, and the access functions will sort out which files/layers are needed for given query. 

On my workstation it takes about 5 seconds to load the entire time series for a given variable - this includes all grid points from 2005 to present + 5 days. The times are seamlessly drawn from multiple files (eg. for `pcp_total` there are 9 different files, representing different steps in the extraction/imputation process described above), and the result is complete, containing no missing times or NA values whatsoever.

It should be easy to export this completed series to a single monolithic NetCDF file, or CSV, or whatever you like. One exception is GeoTIFF, which has limit of about 65,000 layers (times) because it uses a 16 bit layer index. 
