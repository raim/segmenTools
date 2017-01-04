#!/usr/bin/env Rscript

## ANNOTATING SEGMENTS BY GENOMIC FEATURES THEY COVER

library(segmenTools)
##segtools <- "~/programs/segmenTools/"
##source(file.path(segtools,"R/segmenTools.R")) # for segment analysis
##source(file.path(segtools,"R/coor2index.R")) # coor2index

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")

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
#  make_option(c("--dcol"), type="character",
#              default=c("query,intersect,qlen,qpos"), 
#              help="overlap statistics to copy to best matching target if details is set to TRUE"),
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search matches on reverse strand"),
  make_option(c("--only.best"), action="store_true", default=FALSE,
              help="include only the top-ranking query hit (highest jaccard=intersect/union); if FALSE all matching queries will be collapsed into ;-separated lists; multiple best hits for one target will always be collapsed into ;-separated lists"),
  make_option(c("--each.hit"), action="store_true", default=FALSE,
              help="multiple query hit will be exported as separate rows instead of collapsing them into ;-separated lists"),
  ## TARGET OPTIONS
  make_option(c("-t", "--target"), type="character", default="", 
              help="target set of chromosomal segments, stdin is used if missing, allowing for command line pipes"),    
  make_option(c("--ttypes"), type="character", default="", 
              help="sub-sets of testset to annotate"),
  make_option(c("--ttypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
  make_option(c("--tcol"), type="character", default="", 
              help="columns in targets to write to result table"),
  make_option(c("--include.empty"), action="store_true", default=FALSE,
              help="include targets without query hits from result table; if TRUE and each.hit is FALSE, the result table will have matching rows with the target table"),
  ## OUTPUT
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
#              dcol="character", # TODO: implement?
              tcol="character")
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

## messages
msg <- function(x)
    cat(x, file=msgfile)

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
query <- read.table(query,sep="\t",header=TRUE, stringsAsFactors=FALSE)

if ( verb>0 ) msg(paste("Loading target:", target, "\t\n"))

if ( target=="" ) {
    target <- file("stdin")
}

target <- read.table(target,sep="\t",header=TRUE, stringsAsFactors=FALSE)

## filter by type
if ( length(qtypes)>0 )
  query <- query[query[,qtypcol]%in%qtypes,]
if ( length(ttypes)>0 )
  target <- target[target[,ttypcol]%in%ttypes,]

if ( verb>0 )
    msg(paste("Loaded", nrow(target), "TARGETS &\t\n",
              "      ", nrow(query), "QUERIES\t\n"))

if ( nrow(query)==0 | nrow(target)==0 )
    stop("Empty query (",nrow(query),") or target (", nrow(target), ")")

## converting both to continuous index
query <- coor2index(query, chrS)
target <- coor2index(target, chrS)

## search on other strand 
if ( antisense ) 
    query <- switchStrand(query, chrS)



if ( verb>0 )
    msg(paste("CALCULATE OVERLAPS\t",time(),"\n",sep=""))

## TODO: allow upstream/downstream ranges
## TODO: allow collapse as argument / requires to add row number
## of target and use merge
result <- annotateTarget(query=query, target=target,
                         collapse=!each.hit, 
                         details=details, only.best=only.best,
                         qcol=qcol, prefix=prefix, msgfile=msgfile)

## TODO: QUALITY FILTERS FOR RESULT?

if ( verb>0 )
    msg(paste("TRANSLATE COORDINATES\t",time(),"\n",sep=""))

## TRANSLATE LEFT/RIGHT TO UPSTREAM/DOWNSTREAM
## convert back to chromosome coordinates

## bind coordinates and results to map coordinates and overlaps
resCol <- colnames(result) # store requested query/result columns
trgCol <-  ifelse(prefix=="", "target", paste(paste(prefix,"target",sep="_")))
resCol <- resCol[resCol!=trgCol]
tidx <- as.numeric(result[,trgCol])

## TODO: only required if details
## TODO: reduce tmp to coordinate columns and add these to options
tmp <- cbind(target[tidx,], result[,-which(colnames(result)==trgCol)]) 

## TODO: handle strands better here! Expecting +/- factors
if ( details ) {
    relCol <- ifelse(prefix=="", "qpos",
                     paste(paste(prefix,"qpos",sep="_")))
    tmp <- index2coor(tmp, chrS, strands=c("+","-"), relCol=relCol) 
    ## TRANSLATE RELATIVE POSITION TO TARGET POSITION
    orig <- strsplit(tmp[,relCol],";")
    new <- unlist(lapply(orig, function(x) {
        new <- x
        new[x=="inside"] <- "covers"
        new[x=="covers"] <- "inside"
        new[x=="upstream"] <- "downstream"
        new[x=="downstream"] <- "upstream"
        paste(new,collapse=";")}))
    new[new=="NA"] <- NA
    tmp[,relCol] <- new
} else {
    tmp <- index2coor(tmp, chrS, strands=c("+","-"))
}


## FINAL RESULT TABLE
## select target columns
if ( length(tcol)==0 )
    tcol <- colnames(target)
## and bind selected target and selected query/result columns
result <- cbind(target[tidx,tcol,drop=FALSE], tmp[,resCol,drop=FALSE])

## include empty?
if ( !include.empty ) {
    qCol <- ifelse(prefix=="", "qlen",
                   paste(paste(prefix,"qlen",sep="_")))
    result <- result[result[,qCol]!=0,]
}

## final coordinate mapping
## NOTE: without relCol !!
## TODO: allow coorCols as option!
coorCols <- c("start", "end", "coor", chrCol = "chr", strandCol = "strand")
if ( any(colnames(result)%in%coorCols) )
    result <- index2coor(result, chrS, strands=c("+","-"))

if ( verb>0 )
    msg(paste("DONE. WRITING RESULTS\t",time(),"\n",sep=""))

if ( verb>0 )
    msg(paste("Writing result:", outfile, "\t\n"))
write.table(result, file=outfile, sep="\t",quote=FALSE,row.names=FALSE, na="")
