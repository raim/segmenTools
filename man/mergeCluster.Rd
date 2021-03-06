% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterTools.R
\name{mergeCluster}
\alias{mergeCluster}
\title{using \pkg{flowMerge} to merge clusterings}
\usage{
mergeCluster(tset, cset, selected)
}
\arguments{
\item{tset}{the `timeseries' object from segmenTier's
\code{\link[segmenTier:processTimeseries]{processTimeseries}} used
for initial clustering}

\item{cset}{the `clustering' object from segmenTier's 
\code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}}

\item{selected}{optional pre-selection of clustering to merge; if missing
the pre-selected clustering (usually max. BIC)  in \code{cset} will be used}
}
\description{
using \pkg{flowMerge} to merge clusterings
}
