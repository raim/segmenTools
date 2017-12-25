#!/usr/bin/Rscript

## script to convert a csv table to an RData file
## to allow faster data loading

verb <- FALSE
header <- TRUE
na2zero <- FALSE

args  <- R.utils::commandArgs(excludeReserved=TRUE, asValues=TRUE)
for ( arg in names(args)[2:length(args)] ) {
  ## fix(?) since R 3.0.0 ?
  if ( is.null(args[[arg]]) ) args[arg] <- TRUE
  assign(arg, as.character(args[arg]))
}

if ( !exists("o", mode="character") ) o <- sub(".csv$",".RData",i)
if ( !exists("sep", mode="character") ) sep <- "\t"

verb <- as.logical(verb)
header <- as.logical(header)
na2zero <- as.logical(na2zero)

if ( verb ) cat(paste("read table:", i, "\n"))
dat <- read.table(i, sep=sep, header=header)

if ( verb ) cat(paste("process table\n"))
dat <- as.matrix(dat)
if ( na2zero ) dat[is.na(dat)] <- 0

if ( verb ) cat(paste("store RData:", o, "\n"))
save('dat', file=o)
