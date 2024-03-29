% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fingerprint_tools.R
\name{calculate_fingerprints_from_inchi}
\alias{calculate_fingerprints_from_inchi}
\title{Calculate molecular fingerprints from inchi}
\usage{
calculate_fingerprints_from_inchi(inchi, fp_type = "maccs",
  fp_mode = "bit", verbose = FALSE, ...)
}
\arguments{
\item{inchi}{list of strings (1 x n_samples)}

\item{fp_type}{string, fingerprint definition}

\item{fp_mode}{string, which fingerprint type should be returned:
\itemize{
  \item "bit": binary fingerprint-vector
  \item "raw": (key, value)-pairs
  \itemize{
    \item key: the ftp definition, e.g. SMARTS string
    \item value: the number of occurances of fpt
  }
  \item "count": counting fingerprint-vector
}}

\item{verbose}{boolean, should the function be verbose}

\item{...}{parameters passed to the fingerprint calculation function}
}
\value{
list of \code{\link[fingerprint]{fingerprint-class}} objects for the provided
    InChIs. \code{is.null(list[i])} is \code{TRUE}, if the fingerprint for
    the corresponding InChI could not be computed.

put an example for counting maccs fingerprints here.
}
\description{
This function takes a set of inchis and calculates for those the desired
set of fingerprints. For a list of valid fingerprint definitions see
\code{\link[rcdk]{get.fingerpint}}.
}
