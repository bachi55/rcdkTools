% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptor_tools.R
\name{descriptors_to_matrix}
\alias{descriptors_to_matrix}
\title{Construct a descriptor matrix from \code{calculate_descriptor_from_*} function}
\usage{
descriptors_to_matrix(desc)
}
\arguments{
\item{desc}{Output of \code{calculate_descriptor_from_*} functions, i.e. the
descriptor values for a set of (n_samples, ) molecules.}
}
\value{
fingerprint matrix (n_samples x n_desc_vals)
}
\description{
Given a list of descriptor values, i.e., the values of a _single_ descriptor
for a _set_ of molecules, this function will return a matrix object containing
all values.
}
