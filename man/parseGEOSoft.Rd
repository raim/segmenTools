% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.R
\name{parseGEOSoft}
\alias{parseGEOSoft}
\title{parse GO Soft archives}
\usage{
parseGEOSoft(
  file,
  idcol,
  exp = 1,
  only.info = FALSE,
  only.data = TRUE,
  valcol = "VALUE",
  title = TRUE,
  desc = TRUE
)
}
\arguments{
\item{file}{a GEO Soft archive file (eg. GSE18902_family.soft.gz)}

\item{idcol}{columns to retrieve from the platform table}

\item{exp}{number of the platform to retrieve from a GEO file
with multiple platforms. TODO: allow to parse
multiple platforms.}

\item{only.info}{only parse information on probes and experiments, but
skip the actual data}

\item{only.data}{skip all samples without data}

\item{valcol}{value column (currently only 1 is allowed) in the sample tables}

\item{title}{if TRUE a sample title will be retrieved}

\item{desc}{if TRUE, value descriptions in field "#<valcol>" will be retrieved}
}
\description{
Parses GEO Soft microarray experiment archives into a
data matrix and an ID-probe mapping table. Probes can be further
summarized for features with \code{\link{summarizeGEOSoft}}.
}
\seealso{
\code{\link{summarizeGEOSoft}}
}
