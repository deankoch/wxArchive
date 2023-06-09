% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nc_aggregate_time.R
\name{nc_aggregate_time}
\alias{nc_aggregate_time}
\title{Aggregate the completed 2-hourly time series to daily}
\usage{
nc_aggregate_time(
  base_dir,
  var_nm = NULL,
  output_nm = .nm_daily,
  fun = "mean",
  tz = "UTC",
  origin_hour = 0L,
  from = NULL,
  to = NULL
)
}
\arguments{
\item{base_dir}{path to parent directory of GRIB storage subfolder}

\item{var_nm}{character vector, the variable(s) to export (NULL to export all)}

\item{output_nm}{character, the sub-directory name for output}

\item{fun}{function, a function name like "mean", "min", or "max" (see \code{?.nc_aggregate_time})}

\item{tz}{character time zone for \code{origin_hour}}

\item{origin_hour}{integer, steps are aligned to start at this hour of the day}

\item{from}{POSIXct, the start of the time range to process}

\item{to}{POSIXct, the end of the time range to process}
}
\value{
vector of file paths written to the \code{output_nm} directory
}
\description{
This creates a (single) NetCDF file in \code{output_nm} for each of the variables named in
\code{var_nm} by applying the function named in \code{fun} to aggregate by day. \code{var_nm} should
be a subset of the names returned by \code{nc_list(base_dir)}.
}
\details{
If \code{fun=NULL}, the function does no aggregation and writes all times to the output file.

If \code{fun} is not \code{NULL} then the last three arguments are passed to \code{.nc_aggregate_time} to
control the alignment of the aggregation window. \code{fun} specifies the function to use for
combining times within the window (see \code{?.nc_aggregate_time}).  Set \code{tz} to the desired
output time zone, and leave \code{origin_hour=0} to have each day begin at 12AM (in time zone
\code{tz}).

Use \code{from} and \code{to} to specify a date range to update (inclusive), or leave them \code{NULL}
to use a default range. The default range is meant to all layers originating from GFS
(allowing them to be can be replaced by newly added RAP layers, or more recently released
GFS forecasts).

The default for \code{to} is always the latest available date in the input. The default for
\code{from} is 10 days before the latest date found in the existing output files. If there
are no existing outputs, the default is set to the earliest available date in the input.

File names for aggregate data are given the suffix \verb{_daily_<fun>} - eg with \code{fun='mean'},
the output "tmp.nc" becomes "tmp_daily_mean.nc".
}
