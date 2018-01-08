#!/usr/bin/env Rscript

## ANNOTATING SEGMENTS BY GENOMIC FEATURES THEY COVER

library(segmenTools)
##segtools <- "~/programs/segmenTools/"
##source(file.path(segtools,"R/segmenTools.R")) # for segment analysis
##source(file.path(segtools,"R/coor2index.R")) # coor2index

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=msgfile)

## required options:
## 1: input segments 'allsegs.csv'
## 2: annotation, use a keyword from genomeData 'annotation' or load table/bed;
##    select types (ORF, etc.)
## 3: upstream and downstream ranges of targets?
## output:
## input segments with added columns, giving semicolon-separated lists
## of covered features; mutual coverages and
## relative position (covers,left,right,inside)

## output feeds into test:
## calculate cumulative dist function (rcdf and rrcdf); and get recovery
## number


suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("--chrfile"), type="character", default="",
              help="chromosome index file, providing a sorted list of chromosomes and their lengths in column 3 [default %default]"),
  ## QUERY OPTIONS
  make_option(c("-q", "--query"), type="character", default="", 
              help="query set of chromosomal segments"),    
  make_option(c("--qclass"), type="character", default="", 
              help="query classes to test"),
  ## TARGET OPTIONS
  make_option(c("-t", "--target"), type="character", default="", 
              help="target set of chromosomal segments, stdin is used if missing, allowing for command line pipes"),    
  make_option(c("--tclass"), type="character", default="", 
              help="target classes to test"),
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search target matches on reverse strand (if target is empty; search will be done for sense query vs. antisense query!"),
  ## OUTPUT
  make_option(c("-o", "--outfile"), type="character", default="", 
              help="file name to write annotated target list"),
  make_option(c("-v", "--verb"), type="integer", default=1, 
              help="0: silent, 1: main messages, 2: warning messages"))

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## process comma-separated list arguments
lst.args <- c()
for ( i in 1:length(lst.args) ) {
    idx <- which(names(opt)==names(lst.args)[i])
    opt[[idx]] <- unlist(strsplit(opt[[idx]], ","))
    for ( j in 1:length(opt[[idx]]) ) {
        tmp <- strsplit(opt[[idx]][j], ":")
    }
    mode(opt[[idx]]) <- lst.args[i]
}
## promote options to main environment 
for ( i in 1:length(opt) ) {
    arg <- names(opt)[i]
    assign(arg, opt[[arg]])
}

## select file pointer for output and messages
if ( outfile=="" ) {
    outfile <- stdout()
    msgfile <- stderr()
} else {
    msgfile <- stdout()
}


## print out arguments
if ( verb>0 )
    msg(paste("SETTINGS:\n"))
for ( i in 1:length(opt) ) {
    if ( verb>0 )
        msg(paste(names(opt)[i], "\t", 
                  paste(opt[[i]],collapse=", "), "\n",sep=""))
}
if ( verb>0 )
    msg(paste("\n"))

if ( verb>0 )
    msg(paste("LOADING DATA FILES\t",time(),"\n",sep=""))

## load chromosome index - DOESNT WORK WITHOUT
if ( verb>0 )
    msg(paste("Loading chromosome index file:", chrfile, "\t\n"))
cf <- read.table(chrfile,sep="\t",header=FALSE)
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos

## READ SEGMENTS TO BE TESTED 
if ( verb>0 ) msg(paste("Loading query:", query, "\t\n"))
query <- read.delim(query, stringsAsFactors=FALSE)

if ( verb>0 ) msg(paste("Loading target:", target, "\t\n"))

## target = antisense of query?
samesame <- FALSE # only do forward or reverse strand if testing against self?
if ( target=="" & antisense ) {
    target <- query
    samesame <- TRUE
} else {
    target <- read.delim(target, stringsAsFactors=FALSE)
}

if ( verb>0 )
    msg(paste("TARGETS\t", nrow(target), "\n",
              "QUERIES\t", nrow(query), "\n",sep=""))

if ( nrow(query)==0 | nrow(target)==0 )
    stop("Empty query (",nrow(query),") or target (", nrow(target), ")")

## TODO: convert to function in R/segmenTools from here:

## converting both to continuous index
query <- coor2index(query, chrS)
target <- coor2index(target, chrS)

## search on other strand 
if ( antisense ) 
    target <- switchStrand(target, chrS)

if ( verb>0 )
    msg(paste("CALCULATE OVERLAPS\t",time(),"\n",sep=""))

## TODO: calculate J_real>J_random
## for all query classes vs. all target classes
## allowing for sense:antisense test WITHIN same segments (query=target)

## query classes
qcls <- as.character(query[,qclass])
qcls.srt <- unique(qcls)
qN <- length(qcls.srt)

## target classes
tcls <- as.character(target[,tclass])
tcls.srt <- unique(tcls)
tN <- length(tcls.srt)

## get full ranges for all query classes
qcls.rng <- rep(list(NA), qN)
names(qcls.rng) <- qcls.srt
for ( i in 1:qN ) {
    cl <- qcls.srt[i]
    idx <- qcls==cl
    rng <- apply(as.matrix(query[idx,c("start","end")]), 1,
                 function(x) x["start"]:x["end"])
    names(rng) <- NULL
    qcls.rng[[cl]] <- unlist(rng)
}

## get full ranges for all target classes
tcls.rng <- rep(list(NA), tN)
names(tcls.rng) <- tcls.srt
for ( i in 1:tN ) {
    cl <- tcls.srt[i]
    idx <- tcls==cl
    rng <- apply(as.matrix(target[idx,c("start","end")]), 1,
                 function(x) x["start"]:x["end"])
    names(rng) <- NULL
    tcls.rng[[cl]] <- unlist(rng)
}

## get intersect/union of all query:target class pairs
jc <- matrix(NA, nrow=qN, ncol=tN)
for ( i in 1:qN ) {
    for ( j in 1:tN ) { ## asymmetric only for target = antisense of query?
        is <- length(intersect(qcls.rng[[i]],tcls.rng[[j]]))
        un <- length(union(qcls.rng[[i]],tcls.rng[[j]]))
        jc[i,j] <- is/un
    }
}

## randomize queries

## get intersegment distributions

## TODO: same as above but after permutation of all segments
## (or only queries)?

