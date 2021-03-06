% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterTools.R
\name{reCluster}
\alias{reCluster}
\title{re-cluster clustering by \code{\link[stats:kmeans]{kmeans}}}
\usage{
reCluster(tset, cset, k, select = TRUE, ...)
}
\arguments{
\item{tset}{the `timeseries' object from segmenTier's
\code{\link[segmenTier:processTimeseries]{processTimeseries}} used
for initial clustering}

\item{cset}{the `clustering' object from segmenTier's 
\code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}}

\item{k}{colum name or index of the clustering that should be
re-clustered; defaults to the pre-selected clustering if missing}

\item{select}{use the re-clustered clustering as the new pre-selected
clustering}

\item{...}{parameters to \code{\link[stats:kmeans]{kmeans}}}
}
\description{
Use cluster centers from an initial clustering to initialize
\code{\link[stats:kmeans]{kmeans}}. This is still experimental,
and used to re-associated data rows to cluster centers from
a best clustering found by
\code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}.
While the latter clustering works best to extract specific time-courses
from the data set, it often comes with a high fraction of badly
associated individual data sets. Re-clustering with
\code{\link[stats:kmeans]{kmeans}} seems to clean this up, e.g., the
phase distributions of re-clustered clusterings are often tighter.
TODO: allow to generate cluster centers from novel data, to
account for different/more data then during clustering!
}
