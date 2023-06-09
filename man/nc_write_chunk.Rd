% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nc_chunk.R
\name{nc_write_chunk}
\alias{nc_write_chunk}
\title{Write time series data to NetCDF files that are chunked by year}
\usage{
nc_write_chunk(r, p, path_only = FALSE, insert = FALSE)
}
\arguments{
\item{r}{SpatRaster to write}

\item{p}{path to the output directory}

\item{path_only}{logical, if TRUE the function creates the directory but writes nothing to it}

\item{insert}{logical, enables replacement of existing times (passed to \code{nc_write})}
}
\value{
a list of times (the result of \code{nc_write} for each year in \code{r})
}
\description{
A wrapper for \code{nc_write} with multiple output files. This splits the data
from \code{r} by year and writes a separate file for each year in directory \code{p}
according to the naming scheme in \code{?nc_chunk}.
}
\details{
\code{p} must be a directory and its name must end with ".nc"

Set \code{name_only} to return the file names that would be written, but not
actually write anything.

Note that the year for a given time can depend on the time zone. Years are
delineated by this function in UTC.
}
