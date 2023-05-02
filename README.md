RAP/RUC/GFS database
Dean Koch, May 02 2023
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This collects 2-hourly data from the following (NOAA) forecast models: 

* Rapid Update Cycle (RUC) from 2005-2007 (26km resolution)
* Rapid Refresh (RAP) from 2007-present less 2 days (13km resolution)
* Global Forecast System (GFS) from 10 days before present to 5 days ahead (28km resolution)

For past times, a grid of data covering all of North America is downloaded. For future times, only the subgrid covering our area of interest (AOI) is downloaded. All past time points use the 1-hour-ahead forecast (if available), and otherwise revert to the 0-hour forecast.

Source files are first saved to disk in their original GRIB/GRIB2 format, before being coverted to NetCDF v4 for easier access. Each source file contains about 310 variables (eg temperature at various altitudes), but only a small number of relevant ones are copied to NetCDF, and only the sub-grid covering the AOI is included.

Grids

All resampled data use coordinates from the LCC projection, matching the 13km RAP grid. The projection WKT is as follows:

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

To harmonize the spatial layout of these three grids, the 26km and 28km datasets are automatically resampled to 13km resolution using bilinear averaging in GDAL.


Missing Data

About 10% of the time points are missing (as of 2023-05-02). This can happen for a number of different reasons: a forecast might never get been published due to operational issues; or a file might simply not have been posted to the NCEI/NCEP archive server due to resource limitations or human error.

Network issues can also lead to failed downloads or corrupted GRIB files, and this is not unlikely given the large number of file transfers (around 80,000).

Time points missing for any of these reasons are automatically imputed by one of three methods, attempted in this order:

1. derive the missing variable by combining other variables (eg precip from its components)
2. downscale the RUC layer (if present) to replace a missing RAP layer
3. Use a fitted AR2 model with seasonal trends to estimate the remaining missing points


Updates

The update script can be run as often as you like. Running it once daily should be sufficient. It will complete a one-day update within about 30 minutes. The script does the following:

* check existing and available RAP/RUC files and download any (applicable) new ones
* clean the data and and append to existing NetCDF files
* imputation: compute total precipitation as sum of component precipitation (as needed)
* imputation: resample to 13km (as needed)
* imputation: predict missing times using AR2 and seasonal trend (as needed)

The first two steps are then repeated for GFS. So far I have not encountered any missing data issues with GFS/ If we do, it would be straightforward to implement the same gap-filling routine as with RAP/RUC.


Files

The total data volume is about 1.3 TB (as of 2023-05-02). The vast majority of this space is used up by source files in GRIB2 format. The output files in NetCDF format currently make up around 6 GB or about 0.5% of the total. Each daily update introduces about 250 MB of new source files.

Each NetCDF file contains data for a single variable, but can include numerous time points. To help speed up access in files with many times, every NetCDF file has an associated JSON file that indexes times and missing data points. 

The NetCDF data for a given variable <var> is spread out over several files, using the following directory structure:

* rap/coarse/<var>.nc stores layers found in RUC
* rap/coarse_resampled/<var>.nc stores layers from RUC resampled to 13km resolution
* rap/fine/<var>.nc stores layers found in RAP (13km resolution)
* rap/completed/<var>.nc stores layers imputed using the AR(2) model (at 13km resolution)

To speed up access and updates, non-imputed times prior to 2023-03 are held in these "long term storage" files:

* rap/coarse_lts/<var>.nc
* rap/fine_lts/<var>.nc

The GFS component doesn't involve imputation, so it has a simpler structure:

* gfs/coarse/<var>.nc stores layers at their original resolution
* gfs/coarse_resampled/<var>.nc stores layers resampled to 13km resolution


Reading/Exporting

Data access and write functions (in R) are designed to operate on the the whole grid at a specified subset of times, with support for chunking time across different files. This means you can have different sets of times in different files, and the access functions will sort out which files/layers are needed for given query.

On my workstation it takes about 5 seconds to load the entire time series (2005 to present + 5 days) for a given variable (all grid points). It would be very easy to export this completed series to a single monolithic NetCDF file, or CSV, or whatever you like. One exception is GeoTIFF, which has limit of about 65,000 layers (times) because it uses a 16 bit layer index. 
