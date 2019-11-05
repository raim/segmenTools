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
  make_option(c("--chrmap"), action="store_true", default=FALSE,
              help="chromosome are names, and will be mapped to index using the first column of the chromosome index file (--chrfile)"),
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
  make_option(c("--add.qstrand"), type="character", default="", 
              help="use this character as strand information in queries"),
  make_option(c("--add.tstrand"), type="character", default="", 
              help="use this character as strand information in targets"),
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search matches on reverse strand"),
  make_option(c("--qrange"), type="character", default="",  # TODO: implement!
              help="instead of start:end coordinates, use a range around either start (start;-1000:100) or end (end;-1000:100) coordinates"),
  ## divergent can be done with a range on target
#  make_option(c("--divergent"), type="integer", default=0,
#              help="if divergent > 0: search matches on reverse strand transcribed in divergent direction, with maximal distance between pairs set by the argument (>0); can be combined with argument shift to include overlapping divergent pairs [default: %default]"),
  make_option(c("--only.best"), action="store_true", default=FALSE,
              help="include only the top-ranking query hit (highest jaccard=intersect/union); if FALSE all matching queries will be collapsed into ;-separated lists; multiple best hits for one target will always be collapsed into ;-separated lists"),
  make_option(c("--each.query"), action="store_true", default=FALSE,
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
  make_option(c("--each.target"), action="store_true", default=FALSE,
              help="include targets without query hits from result table; if TRUE and each.query is FALSE, the result table will have matching rows with the target table"),
  ## diverse parsing/handling options
  make_option(c("--cchar"), type="character", default="", 
              help="comment character in input files"),  
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


## READ SEGMENTS TO BE TESTED 
if ( verb>0 ) msg(paste("Loading query:", query, "\t\n"))
query <- read.table(query, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                    comment.char=cchar, quote="")

if ( verb>0 ) msg(paste("Loading target:", target, "\t\n"))

if ( target=="" ) {
    target <- file("stdin")
}

target <- read.table(target, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                     comment.char=cchar, quote="")

## filter by type
if ( length(qtypes)>0 )
  query <- query[query[,qtypcol]%in%qtypes,]
if ( length(ttypes)>0 )
  target <- target[target[,ttypcol]%in%ttypes,]

if ( verb>0 )
    msg(paste("TARGETS\t", nrow(target), "\n",
              "QUERIES\t", nrow(query), "\n",sep=""))

if ( nrow(query)==0 | nrow(target)==0 )
    stop("Empty query (",nrow(query),") or target (", nrow(target), ")")

## adding strand if missing
if ( add.qstrand!="" )
    if ( !"strand"%in%colnames(query) )
        query$strand <- add.qstrand
if ( add.tstrand!="" )
    if ( !"strand"%in%colnames(query) )
        target$strand <- add.tstrand

## load chromosome index - DOESNT WORK WITHOUT
## NOTE: file with sorted chromosomes and their lengths
## in the last column; segment query and target files
## refer to these chromosomes in column "chr"
if ( verb>0 )
    msg(paste("Loading chromosome index file:", chrfile, "\t\n"))
cf <- read.table(chrfile,sep="\t",header=FALSE, stringsAsFactors=FALSE, quote="")
chrS <- c(0,cumsum(as.numeric(cf[,ncol(cf)]))) ## index of chr/pos = chrS[chr] + pos

## name information in first column
chrMap <- cf[,1]

## TODO: segmenTools wrapper starting from un-indexed chromosome
## searching for specific rules (tandem; convergent/antisense/divergent)
## with self or target
## TODO: convert to function in R/segmenTools from here:

## converting both to continuous index
## TODO: use chrMap if chromosome columns are names and don't use
## the index
if ( chrmap ) {
    query <- coor2index(query, chrS, chrMap=cf[,1])
    target <- coor2index(target, chrS, chrMap=cf[,1])
} else {
    query <- coor2index(query, chrS)
    target <- coor2index(target, chrS)
}

## search on other strand 
if ( antisense ) 
    query <- switchStrand(query, chrS)

if ( verb>0 )
    msg(paste("CALCULATE OVERLAPS\t",time(),"\n",sep=""))

## check if requested columns are present
## TODO: check this early in annotateTarget as well
if ( any(!qcol%in%colnames(query)) ) {
    warning(paste(qcol[!qcol%in%colnames(query)],collapse=";"),
            " not in query columns (argument --qcol)! Skipped!\n")
    qcol <- qcol[qcol%in%colnames(query)]
    if ( length(qcol)==0 )
        stop("no query columns (argument --qcol) found.")
}

## TODO: allow upstream/downstream ranges
## TODO: allow collapse as argument / requires to add row number
## of target and use merge
result <- annotateTarget(query=query, target=target,
                         collapse=!each.query, 
                         details=details, only.best=only.best,
                         qcol=qcol, prefix=prefix, msgfile=msgfile)
#save.image("tmp.RData")
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
tmp <- cbind(target[tidx,],
             result[,-which(colnames(result)==trgCol),drop=FALSE]) 

## TODO: handle strands better here! Expecting +/- factors
## TODO - 20180320: smarter handling of antisense option?
##                  currently assumes same strand
##                  also: handle convergent/divergent, tandem neighbors
if ( details ) {
    relCol <- ifelse(prefix=="", "qpos",
                     paste(paste(prefix,"qpos",sep="_")))
    if ( chrmap ) 
        tmp <- index2coor(tmp, chrS, strands=c("+","-"), relCol=relCol, chrMap=chrMap)
    else
        tmp <- index2coor(tmp, chrS, strands=c("+","-"), relCol=relCol) 
    ## TRANSLATE RELATIVE POSITION TO TARGET POSITION
    orig <- strsplit(as.character(tmp[,relCol]),";")
    new <- unlist(lapply(orig, function(x) {
        new <- x
        new[x=="inside"] <- "covers"
        new[x=="covers"] <- "inside"
        new[x=="upstream"] <- "downstream"
        new[x=="downstream"] <- "upstream"
        paste(new,collapse=";")}))
    new[new=="NA"] <- NA
    tmp[,relCol] <- new
} 

## FINAL RESULT TABLE
## select target columns
if ( length(tcol)==0 )
    tcol <- colnames(target)
## and bind selected target and selected query/result columns
result <- cbind(target[tidx,tcol,drop=FALSE], tmp[,resCol,drop=FALSE])

## remove empty targets (no hit)
## TODO: only works well if details==TRUE, since query length field can be used
##       implement this smarter for details==FALSE !
if ( !each.target )
    if (  details ) {
        qCol <- ifelse(prefix=="", "qlen",
                   paste(paste(prefix,"qlen",sep="_")))
        result <- result[as.character(result[,qCol])!="0",] # char allows collapse
    } else{
        #qCol <- ifelse(prefix=="", qcol[1],
        #           paste(paste(prefix,qcol[1],sep="_")))
        #result <- result[as.character(result[,qCol])!="",] # char allows collapse
    }


## final coordinate mapping
## NOTE: without relCol !!
## TODO: allow coorCols as option!
coorCols <- c("start", "end", "coor", chrCol = "chr", strandCol = "strand")
if ( any(colnames(result)%in%coorCols) )
    if ( chrmap ) {
        result <- index2coor(result, chrS,
                             strands=c("+","-"), chrMap=chrMap)
    } else {
        result <- index2coor(result, chrS, strands=c("+","-"))
    }

if ( verb>0 )
    msg(paste("DONE. WRITING RESULTS\t",time(),"\n",sep=""))
if ( verb>0 )
    msg(paste("Writing result:", outfile, "\t\n"))

write.table(result, file=outfile, sep="\t",quote=FALSE,row.names=FALSE, na="")

if ( verb>0 )
    msg(paste("RESULTS\t", nrow(result), "\n",sep=""))
