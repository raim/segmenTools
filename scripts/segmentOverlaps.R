#!/usr/bin/env Rscript

## SEGMENT CLASSES OVERLAP STATISTIC
## by permutation test analysis

library(segmenTools)

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=msgfile)


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
              help="name of column target classes to test"),
  make_option(c("--ttypes"), type="character", default="", 
              help="sub-sets of testset to use for testing"),
  make_option(c("--ttypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search target matches on reverse strand (if target is empty; search will be done for sense query vs. antisense query!"),
  make_option(c("--upstream"), type="integer", default=0,
              help="search range upstream of target (in nt.)"),
  make_option(c("--perm"), type="integer", default=100, 
              help="number of permutations"),
  ## OUTPUT
  make_option(c("--fig.type"), type="character", default="png",
              help="figure type (png, pdf, eps) [default %default]"),
  make_option(c("-o", "--outfile"), type="character", default="", 
              help="file name to write annotated target list"),
  make_option(c("-v", "--verb"), type="integer", default=1, 
              help="0: silent, 1: main messages, 2: warning messages"))

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## process comma-separated list arguments
lst.args <- c(ttypes="character")
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
cf <- read.table(chrfile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
chrMap <- cf[,2]
chrL <- cf[,3]
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos
total <- 2*sum(chrL) # both strands!

## READ SEGMENTS TO BE TESTED 
if ( verb>0 ) msg(paste("Loading query:", query, "\t\n"))
query <- read.delim(query, stringsAsFactors=FALSE)

if ( verb>0 ) msg(paste("Loading target:", target, "\t\n"))

## TODO: align strand columns to - -> -1 and + -> 1 

## plot axis labels
qlab <- paste0("query: ", qclass)
tlab <- paste0("target: ", tclass)

## target = antisense of query?
frw.str <- c("1","+")
rev.str <- c("-1","-")
## comparison with self on reverse strand!
## compare forward with reverse strand!
self <- FALSE
if ( target=="" & (antisense|upstream!=0) ) {
    self <- TRUE
    target <- query
    if ( tclass=="" )
      tclass <- qclass
    if ( antisense ) {
        query <- query[as.character(query[,"strand"])%in%frw.str,]
        target <- target[as.character(target[,"strand"])%in%rev.str,]
        ## total length: only one strand
        total <- sum(chrL)
    }    
} else {
    target <- read.delim(target, stringsAsFactors=FALSE)
    ## TODO: map chromosome name to index
    ## perhaps allow to pass chromosome map to coor2index
    
    ## FILTER targets
    if ( length(ttypes)>0 )
        if( ttypes!="" )
            target <- target[target[,ttypcol]%in%ttypes,]
}

## scan for range around targets
if ( upstream!=0 ) {

    if ( antisense )
        stop("Options `--antisense' and `--upstream' are not compatible!")

    ## make sure start < end for both strands
    ## (only valid for non-circular DNA!!)
    utarget <- target[,c("start","end","strand")]
    str <- as.character(utarget[,"strand"])
    start <- utarget[,"start"]
    end <- utarget[,"end"]
    ## order start<end 
    if ( any(end<start & str%in%frw.str) ) 
        stop("`end' must be larger then segment `start' on forward strand")
    start[str%in%rev.str] <-
        apply(utarget[str%in%rev.str,c("start","end")],1,min)
    end[str%in%rev.str] <- apply(utarget[str%in%rev.str,c("start","end")],1,max)
    
    utarget[str%in%frw.str, c("start","end")] <-
        cbind(start-upstream,start-1)[str%in%frw.str,]
    utarget[str%in%rev.str, c("start","end")] <-
        cbind(end+1,end+upstream)[str%in%rev.str,]
    target[,c("start","end","strand")] <- utarget
}

## only compare forward and reverse strands for auto-target antisense
if ( antisense & self ) {
    ## plot axis labels
    qlab <- paste0("query: ", qclass, ", strand ", frw.str[2])
    tlab <- paste0("target: ", tclass, ", strand ", rev.str[2])
}
if ( antisense & !self ) {
    tlab <- paste(tlab, "- antisense")
}
if ( upstream!=0 ) {
    ## plot axis labels
    qlab <- paste("query:", qclass)
    tlab <- paste("target:", tclass, "- upstream", upstream)
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

## calculate Jaccard Index and permutation test
ovl <- segmentJaccard(query=query, target=target,
                      qclass=qclass, tclass=tclass, perm=perm, total=total,
                      verb=1)

if ( verb>0 )
  msg(paste0("DONE\t",time(),"\n"))
if ( verb>0 )
  msg(paste0("writing results\n"))

## write out results
if ( tclass=="" ) tclass <- "all"
if ( qclass=="" ) qclass <- "all"

file.name <- paste0(outfile,"_",qclass,"_",tclass,
                    ifelse(antisense,"_antisense",""),
                    ifelse(upstream!=0, paste0("_upstream",upstream),""))
## store data
if ( !interactive() ) 
  save(ovl, file=paste0(file.name,".RData"))

## plot
if ( perm>0 ) {
    plotdev(paste0(file.name),type=fig.type)
    plotOverlaps(ovl,p.min=.001,main="Jaccard Index (*1000) & permutation test",ylab=qlab,xlab=tlab,scale=1000,round=0)
    dev.off()
}
