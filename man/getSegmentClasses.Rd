% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmenTools.R
\name{getSegmentClasses}
\alias{getSegmentClasses}
\title{splits segmenTier segment class strings into classes}
\usage{
getSegmentClasses(sgtypes, sep = "_", gsep = ":")
}
\arguments{
\item{sgtypes}{a list of strings, that is converted to a table
of classes based on a string separator (\code{sep})}

\item{sep}{the separator for classes in the string}

\item{gsep}{currently not used, a separator within classes
to be used as classification ID (column header of the produces
class table)}
}
\description{
specifically tailored to segment strings in segmenTier
}
