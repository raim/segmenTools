% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotUtils.R
\name{num2col}
\alias{num2col}
\title{convert numeric values to color range:}
\usage{
num2col(x, limits, q, pal, colf = viridis::viridis, n = 100)
}
\arguments{
\item{x}{a numeric vector}

\item{limits}{optional data limits for min/max color, every x
lower/higher will get the extreme colors}

\item{q}{auto-select limits by these two quantiles (argument
\code{prob} in function \code{\link[stats:quantile]{quantile}}}

\item{pal}{color palette, alternatively \code{colf} and \code{n}
can be supplied}

\item{colf}{color palette function}

\item{n}{number of different colors}
}
\description{
convert numeric values to color range:
}
