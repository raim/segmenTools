% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osciTools.R
\name{get_fft}
\alias{get_fft}
\title{Discrete Fourier Transformation}
\usage{
get_fft(x)
}
\arguments{
\item{x}{data to be transformed}
}
\description{
A simple wrapper around \code{\link[stats:fft]{mvfft}} to perform
Discrete Fourier Transformation, discard redundant components
(for real numbers only, the second half of the transform),
and name columns by `DC` (direct current, a term from electrical
engineering) and the numbers of cycles in the data set reflected
by each component. The last column is the Nyquist frequency.
}
