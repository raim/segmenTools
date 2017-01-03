% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmenTools.R
\name{presegment}
\alias{presegment}
\title{pre-segmentation of whole-genome data into chunks that can
be handle by segmenTier.}
\usage{
presegment(ts, chrS, avg = 1000, favg = 100, minrd = 8, minds = 250,
  map2chrom = FALSE, seg.path, fig.path, fig.type = "png", verb = 1)
}
\arguments{
\item{ts}{the time-series of readcounts for the complete chromosome,
rows are chromosomal positions and columns are time-points; reverse
strand rows at the bottom of the matrix. Option \code{chrS} can be
used to handle chromosome ends and to optionally (\code{map2chrom})
map the resulting primary segment coordinates to chromosome coordinates.}

\item{chrS}{a chromosome index, indicating at wich positions
chromosomes start; this is required for handling chromosome ends
and forward and reverse strand values, but can be omitted.}

\item{avg}{the broad moving average of read-count presence
(number of time-points with >0 reads) for a first broad segmentation}

\item{favg}{as \code{avg}, but a narrower moving average used in
end scanning that can result in fusing back segments w/o good separation}

\item{minrd}{the minimal number of time-points with reads in the broad
moving average used as cutoff between segments}

\item{minds}{minimum distance between two segments (will be fused otherwise)}

\item{map2chrom}{if true, argument \code{chrS} is required to map
the segment coordinates to chromosomal coordinates}

\item{seg.path}{a directory path where individual segments' data will
be written to as tab-delimited .csv files; no files will be written if
\code{seg.path} is not provided.}

\item{fig.path}{a directory path for plots of the segment end scanning;
no figures will be plotted if \code{fig.path} is not provided.}

\item{fig.type}{image type, "png" or "pdf"}
}
\description{
pre-segmentation of whole-genome data into chunks that can
be handle by segmenTier.
}
