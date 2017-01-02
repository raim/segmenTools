#!/usr/bin/env Rscript

## ANNOTATING SEGMENTS BY GENOMIC FEATURES THEY COVER

## genome annotation is loaded via the genomeBrowser/genomeData project
browser.path <- sub("GENBRO=","",system("env|grep GENBRO",intern=TRUE))
browser.data <- sub("GENDAT=","",system("env|grep GENDAT",intern=TRUE))
browser.data <- file.path(browser.data,"yeast")
source(file.path(browser.path,"src/genomeBrowser_utils.R")) # for coor2index
source(file.path(browser.path,"src/segmenTools.R")) # for segment analysis

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
  make_option(c("-q", "--query"), type="character", default="", 
              help="query set of chromosomal segments"),    
  make_option(c("--qtypes"), type="character", default="", 
              help="sub-sets of query to use for annotation"),
  make_option(c("--qtypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
  make_option(c("--qcol"), type="character", default="", 
              help="columns in query to copy to best matching target"),
  make_option(c("--prefix"), type="character", default="", 
              help="optional name prefix for copied columns"),
  make_option(c("--details"), action="store_true", default=FALSE,
              help="add details of the overlap"),
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search matches on reverse strand"),
  make_option(c("-t", "--target"), type="character", default="", 
              help="target set of chromosomal segments"),    
  make_option(c("--ttypes"), type="character", default="", 
              help="sub-sets of testset to annotate"),
  make_option(c("--ttypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
  make_option(c("--tcol"), type="character",
              default=c("query,intersect,qlen,qpos"), 
              help="overlap statistics to copy to best matching target"),
  make_option(c("-o", "--outfile"), type="character", default="", 
              help="file name to write annotated target list"),
  make_option(c("-v", "--verb"), type="integer", default=1, 
              help="0: silent, 1: main messages, 2: warning messages"))

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## process comma-separated list arguments
lst.args <- c(qtypes="character",
              ttypes="character",
              qcol="character",
              tcol="character")
for ( i in 1:length(lst.args) ) {
    idx <- which(names(opt)==names(lst.args)[i])
    opt[[idx]] <- unlist(strsplit(opt[[idx]], ","))
    for ( j in 1:length(opt[[idx]]) ) {
        tmp <- strsplit(opt[[idx]][j], ":")
    }
    mode(opt[[idx]]) <- lst.args[i]
}
## promote options to main environment and print all arguments
if ( opt$verb>0 )
    cat(paste("SETTINGS:\n"))
for ( i in 1:length(opt) ) {
    if ( opt$verb>0 )
        cat(paste("\t",names(opt)[i], ":", #typeof(opt[[i]]),
                  paste(opt[[i]],collapse=", "), "\n"))
    arg <- names(opt)[i]
    assign(arg, opt[[arg]])
}
if ( verb>0 )
    cat(paste("\n"))

## load chromosome index - DOESNT WORK WITHOUT
if ( verb>0 )
    cat(paste("Loading chromosome index file:", chrfile, "\n"))
cf <- read.table(chrfile,sep="\t",header=FALSE)
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos

## READ SEGMENTS TO BE TESTED 
if ( verb>0 ) cat(paste("Loading query:", query, "\n"))
query <- read.table(query,sep="\t",header=TRUE, stringsAsFactors=FALSE)

if ( verb>0 ) cat(paste("Loading target:", target, "\n"))
target <- read.table(target,sep="\t",header=TRUE, stringsAsFactors=FALSE)

## filter by type
if ( length(ttypes)>0 )
  query <- query[query[,qtypcol]%in%qtypes,]
if ( length(ttypes)>0 )
  target <- target[target[,ttypcol]%in%ttypes,]

## converting both to continuous index
query <- coor2index(query, chrS)
target <- coor2index(target, chrS)

## search on other strand 
if ( antisense ) 
    query <- switchStrand(query, chrS)



## TODO: reverse upstream-downstream relative positions?
## TODO: allow upstream/downstream ranges
target <- annotateTarget(query=query, target=target, col=qcol,
                         prefix=prefix, details=details)


## convert back to chromosome coordinates
## TODO: handle strands better here! Expecting +/- factors
if ( details ) {
    relcol <- ifelse(prefix=="", "qpos",
                     paste(paste(prefix,"qpos",sep="_")))
    target <- index2coor(target, chrS, strands=c("+","-"), relCol=relcol) 
} else {
    target <- index2coor(target, chrS, strands=c("+","-"))
}

if ( verb>0 )
    cat(paste("Writing target:", outfile, "\n"))
write.table(target, file=outfile, sep="\t",quote=FALSE,row.names=FALSE)
