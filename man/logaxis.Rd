% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotUtils.R
\name{logaxis}
\alias{logaxis}
\title{Logarithmic axis ticks}
\usage{
logaxis(side, lat = -10:10, base = 10, labels, ...)
}
\arguments{
\item{side}{an integer specifying which side of the plot the axis is
    to be drawn on.  The axis is placed as follows: 1=below,
    2=left, 3=above and 4=right.}

\item{lat}{as argument \code{at} of \link[graphics:axis]{axis} but
as the exponents of \code{log(y, base=base)}}

\item{base}{base of the logarithm}

\item{labels}{this can either be a logical value specifying whether
    (numerical) annotations are to be made at the tickmarks, or a
    character or expression vector of labels to be placed at the
    tickpoints.  (Other objects are coerced by \code{\link{as.graphicsAnnot}}.)
    If this is not logical, \code{at} should also be supplied and of the
    same length.  If \code{labels} is of length zero after coercion,
    it has the same effect as supplying \code{TRUE}.}

\item{...}{other \link[graphics]{graphical parameters} may also be passed as
    arguments to this function, particularly, \code{cex.axis}, \code{col.axis}
    and \code{font.axis} for axis annotation, i.e. tick labels, \code{mgp}
    and \code{xaxp} or \code{yaxp} for positioning, \code{tck} or
    \code{tcl} for tick mark length and direction, \code{las} for
    vertical/horizontal label orientation, or \code{fg} instead of
    \code{col}, and \code{xpd} for clipping.  See \code{\link[graphics]{par}} on these.

    Parameters \code{xaxt} (sides 1 and 3) and \code{yaxt} (sides 2 and
    4) control if the axis is plotted at all.

    Note that \code{lab} will partial match to argument
    \code{labels} unless the latter is also supplied.  (Since the
    default axes have already been set up by \code{\link[graphics]{plot.window}},
    \code{lab} will not be acted on by \code{axis}.)}
}
\description{
Logarithmic axis ticks
}
