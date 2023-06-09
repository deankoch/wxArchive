% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_wx.R
\name{file_wx}
\alias{file_wx}
\title{Build file paths for writing outputs}
\usage{
file_wx(
  what,
  base_dir = "",
  sub_dir = list("coarse", "fine"),
  var_nm = as.list(names(.rap_regex)),
  make_dir = FALSE
)
}
\arguments{
\item{what}{one of 'grib', 'csv', 'nc', 'index', 'spatial', 'temporal' (see details)}

\item{base_dir}{character path the parent directory for RAP file storage}

\item{sub_dir}{character vector or list of them, names of sub-directories in \code{base_dir}}

\item{var_nm}{character vector of variable names, or a list of them}

\item{make_dir}{logical indicating to create directories where needed}
}
\value{
a named character vector, or a named list of them (depending on \code{what} and \code{sub_dir})
}
\description{
This organizes the input/output files for various steps in the RAP workflow. When
\code{make_dir=TRUE}, the function creates the parent directory containing the requested
file name (if it doesn't already exist).
}
\details{
The argument \code{sub_dir} defines the directory structure. eg by default we have:
\itemize{
\item .../grib/ stores the many source files (.grb and .grb2)
\item .../coarse/ stores NetCDF versions (.nc) for each variable at coarse resolution
\item .../fine/ stores NetCDF versions (.nc) for each variable at fine resolution
}

Note that a "file" may refer to a directory with the suffix ".nc" in its name. Such
directories should contain a set of yearly NetCDF files named as defined by \code{nc_chunk}.

Argument \code{what} specifies the path(s) to return:
\enumerate{
\item 'grib' : character, the source GRIB file directory
\item 'nc' : a list of character vectors, paths to the .nc file
\item 'index' : a list of character vectors, paths to JSON files storing time index
\item 'spatial' : a list of character vectors, paths to JSON files storing fitted spatial parameters
\item 'temporal_index' : character, path to a JSON file storing information about temporal fit
\item 'temporal_nc' : character, path to a new .nc file for storing fitted temporal parameters
}

With options (1) and (5-6) the function returns a character path, whereas with (2-4) it
returns a list of character path vectors - one list entry per resolution, and one vector
entry per variable name. If \code{sub_dir} has length 1 (one resolution only) then the result is
automatically unlisted and the function returns a character vector.

With option (6), the file name will contain the current system time - use this only when
creating a new file, then append that file name to the JSON from (5) for later use.

When \code{collapse=TRUE}, the top level of the list result is collapsed by concatenating
its contents using the names from the first list entry. This has no effect when
the output would otherwise have length 1, or when \code{sub_dir} has length 1.

Both \code{sub_dir} and \code{var_nm} can be vectors or a lists of vectors. Lists are used to group
like subdirectories (eg different pieces of the same time series) or variable names that
should be viewed as equivalent. The function loops over lists, passing each element to itself
in a recursive call and collapsing the results. This makes a (possibly nested) list of
vectors, the result of looping over \code{sub_dir} first (outer), then \code{var_nm} (inner).
}
