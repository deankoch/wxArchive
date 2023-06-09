% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow.R
\name{workflow_daily}
\alias{workflow_daily}
\title{Exported completed time series to daily aggregate values and write to disk}
\usage{
workflow_daily(project_dir, from = NULL, to = NULL, tz = "MST")
}
\arguments{
\item{project_dir}{character path to the project root directory}

\item{from}{POSIXct start of time range (default is earliest available)}

\item{to}{POSIXct end of time range (default is latest available)}
}
\value{
returns nothing but possible writes to \code{project_dir}
}
\description{
Wrapper for \code{nc_aggregate_time}
}
\details{
This prepares five output variables aggregated to daily average, or maximum or minimum.
The pairing of variable names and aggregation functions is set up in the global
constant \code{.var_daily_pairs} and the name of the output folder is \code{.nm_daily}.

By default the function overwrites layers starting from 10 days before the latest
date in the existing output files. This ensures that any/all GFS forecasts are
updated with fresh estimates, or else replaced by the more precise RAP archive
estimates. Specify a different time range to process using \code{from} and \code{to}.

The time zone string \code{tz} controls the alignment of the output frames, so that each
"day" starts (and ends) on the hour 12AM in the time zone \code{tz}. Set this to the local
time zone for your AOI.
}
