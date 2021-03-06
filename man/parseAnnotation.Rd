% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterTools.R
\name{parseAnnotation}
\alias{parseAnnotation}
\title{Parse an annotation file (a bidirectional map)}
\usage{
parseAnnotation(got, idcol = 1, keycol = 6, termcol, rm.empty = TRUE)
}
\arguments{
\item{got}{input table, e.g. a GO annotation}

\item{idcol}{column where feature IDs can be found}

\item{keycol}{column where annotation terms are found}

\item{termcol}{optional column where a readable description of annotation
terms is found; will only be used if keycol and termcol are a constant map,
ie. each key is always associated with the same term; TODO: terminology?}

\item{rm.empty}{rm all empty annotations; TODO: why does this occur?}
}
\description{
Parses a bidirectional map of feature IDs vs. annotation terms, e.g.
the GO annotation file at \url{ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt.gz} and returns a TRUE/FALSE table that
can be used as input to \code{\link{clusterAnnotation}}.
}
