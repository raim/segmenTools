% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmenTools.R
\name{segmentMerge}
\alias{segmentMerge}
\title{merge genomic coordinates by type.}
\usage{
segmentMerge(x, type, verb = 1)
}
\arguments{
\item{x}{genomic interval (segments) table.}

\item{type}{column name of the types which should be merged
separately.}

\item{verb}{verbosity level, 0 for silent}
}
\description{
NOTE: merging segments requires bedtools!
TODO: allow merging without bedtools.
}
