% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterTools.R
\name{add_alphas}
\alias{add_alphas}
\title{replace alpha values of an RGB string color vector}
\usage{
add_alphas(col, alpha = rep(1, length(col)))
}
\arguments{
\item{col}{a string vector of RGB color specifiers with ("#FF0000AA ")
or without ("#FF0000") existing alpha values}

\item{alpha}{a numeric vector of equal length as argument \code{col}
providing new alpha values from 0 to 1}
}
\description{
adds (or replaces existing) alpha values (from 0 to 1) to
RGB string colors; extended to vectors from \url{http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html}
}
