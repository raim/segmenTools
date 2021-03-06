% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterTools.R
\name{clusterAverages}
\alias{clusterAverages}
\title{calculates cluster averages}
\usage{
clusterAverages(ts, cls, cls.srt, avg = "median", q = 0.9, rm.inf = TRUE)
}
\arguments{
\item{ts}{a matrix of time series, with time points in columns}

\item{cls}{a clustering of the time series, \code{length(cls)}
must equal \code{nrow(cls)}}

\item{cls.srt}{optional sorting of the clusters}

\item{avg}{a function (or the name of a function as a string)
for calculating an `average' value for each cluster; default is
the \code{median}}

\item{q}{either numeric 0-1, the fraction of data for which high
and low data cut-offs are calculated, or a function name for
calculating variance (eg. "sd", "var"), which will be added and
subtracted from the average (argument \code{avg}) for high
and low data cut-offs}

\item{rm.inf}{remove infinite values (e.g. from log transformations)}
}
\description{
calculates average values and distributions for each cluster
and time point of a time series
}
