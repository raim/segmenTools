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
  make_option(c("--setup"), action="store_true", default=FALSE,
              help="set up queries for the given settings; this will generate randomized queries for re-use in real runs, and use a minimal pseudo target to test the whole pipeline"),
  make_option(c("--merge"), action="store_true", default=FALSE,
              help="don't merge query and target segments (merge is required until bedtools --noOverlapping works)"),
  ## TARGET OPTIONS
  make_option(c("-t", "--target"), type="character", default="", 
              help="target set of chromosomal segments, stdin is used if missing, allowing for command line pipes"),    
  make_option(c("--tclass"), type="character", default="", 
              help="name of column target classes to test"),
  make_option(c("--ttypes"), type="character", default="", 
              help="sub-sets of testset to use for testing"),
  make_option(c("--ttypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
  make_option(c("--nostrand"), action="store_true", default=FALSE,
              help="ignore strand information in query and target"),
  ## SEARCH RANGE OPTIONS
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search target matches on reverse strand (if target is empty; search will be done for sense query vs. antisense query!"),
  make_option(c("--upstream"), type="integer", default=0,
              help="search upstream (<0) or downstream (>0) of target (in bp)"),
  make_option(c("--range"), type="integer", default=0,
              help="search range around target (in bp)"),
  make_option(c("--convergent"), type="integer", default=0,
              help="search convergent (>0) or divergent (<0) overlaps between query and target (in bp)"),
  ## ALGORITHM OPTIONS
  make_option(c("--perm"), type="integer", default=100, 
              help="number of permutations"),
  make_option(c("--count"),  action="store_true", default=FALSE,
              help="count individual overlaps (only required for --bedtools FALSE"),
  make_option(c("--bedtools"),  action="store_true", default=FALSE,
              help="use UCSC bedtools instead of segmenTools overlap functions"),
  make_option(c("--random"),  type="character", default="",
              help="directory where (re-usable) bedtools randomized queries are stored"),
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

## catch incompatible options
if ( antisense & convergent!=0 )
    stop("options --antisense and --convergent are incompatible!")
if ( antisense & upstream!=0 )
    stop("options --antisense and --upstream are incompatible!")
if ( nostrand & antisense )
    stop("options --nostrand and --antisense are incompatible")
if ( nostrand & upstream!=0 )
    stop("options --nostrand and --upstream are incompatible")
if ( nostrand & convergent!=0 )
    stop("options --nostrand and --convergent are incompatible")
if ( upstream!=0 & convergent!=0 )
    stop("options --upstream and --convergent are incompatible")


if ( verb>0 )
    msg(paste("LOADING DATA FILES\t",time(),"\n",sep=""))

## load chromosome index - DOESNT WORK WITHOUT
if ( verb>0 )
    msg(paste("Loading chromosome index file:", chrfile, "\t\n"))
cf <- read.table(chrfile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
#chrMap <- cf[,2]
chrL <- cf[,ncol(cf)]
chrS <- c(0,cumsum(cf[,ncol(cf)])) ## index of chr/pos = chrS[chr] + pos
total <- 2*sum(chrL) # both strands!

## READ SEGMENTS TO BE TESTED 
if ( verb>0 ) msg(paste("Loading query:", query, "\t\n"))
query <- read.delim(query, stringsAsFactors=FALSE)

if ( verb>0 ) msg(paste("Loading target:", target, "\t\n"))

## TODO: align strand columns to - -> -1 and + -> 1 


## target = antisense of query?
frw.str <- c("1","+")
rev.str <- c("-1","-")

## on the fly flags for further processing
mergeq <- FALSE  # should query be merged?
merget <- FALSE  # should target be merged?
self <- FALSE    # symmetric comparison between forward and reverse strand?


## COMPARISON WITH SELF on reverse strand!
## compare forward with reverse strand!
if ( target=="" ) {
    self <- TRUE
    target <- query
    if ( tclass=="" )
      tclass <- qclass
    if ( antisense | convergent!=0 ) {
        ## NOTE: query on forward strand is permutated
        query <- query[as.character(query[,"strand"])%in%frw.str,]
        target <- target[as.character(target[,"strand"])%in%rev.str,]
        ## total length: only one strand
        total <- sum(chrL)
    }  
} else {
    target <- read.delim(target, stringsAsFactors=FALSE)
    
    ## FILTER targets
    if ( length(ttypes)>0 )
        if( all(ttypes!="") ) {
            target <- target[target[,ttypcol]%in%ttypes,]
        } else {
            stop("--ttypes must be a string") 
        }
}

## add type columns, if not provided: all segments are the same class
if ( tclass=="" ) {
    tclass <- "type"
    if ( tclass%in%colnames(target) ) # remove existing
        target <- target[, colnames(target)!="type"]
    target$type <- "all"
}
if ( qclass=="" ) {
    qclass <- "type"
    if ( qclass%in%colnames(query) ) # remove existing
        query <- query[, colnames(query)!="type"]
    query$type <- "all"
}


## scan downstream (>0) or upstream (<0) of target
if ( upstream!=0 ) {

    ## get new coordinates
    target <- segmentUpstream(x=target, upstream=upstream)
    merget <- TRUE
}

## scan for convergent (>0) or divergent (<0) overlaps
if ( convergent!=0 ) {

    ## same as upstream<0 but for both query and target

    target <- segmentUpstream(x=target, upstream=convergent)
    merget <- TRUE

    query <- segmentUpstream(x=query, upstream=convergent)
    mergeq <- TRUE
}


## switch strand for antisense and convergent tests
sstr <- c("+"=-1, "1"=-1, "+1"=-1, "-"=1, "-1"=1)
if ( (antisense|convergent!=0) ) 
    target$strand <- sstr[as.character(target$strand)]
    

## CONSIDER ONLY ONE STRAND:
## either by command line argument, or by absence of strand
## column

if ( !"strand"%in%colnames(query) | !"strand"%in%colnames(target) ) {
    if ( !nostrand )
        warning(paste0("missing strand information, ",
                       "strand information is ignored, --nostrand activated"))
    nostrand <- TRUE
}
if ( nostrand ) {
    total <- total/2

    ## merge if queries or targets from both strands are present
    if ( "strand"%in%colnames(query) )
        mergeq <- ifelse(length(unique(query$strand))>1, TRUE, FALSE)

    query$strand <- "1"
    target$strand <- "1"
}

## PRUNE, SORT AND MERGE both query and target
## to make sure randomizations work properly (assumed non-overlapping!)
## and overlap statistics are exact.

## TODO: prune by default? pruning may fail to catch errors in the input data.
## perhaps, this is too tolerant?

msg(paste("Pruning target:\n"))
target <- segmentPrune(x=target, chrL=chrL, remove.empty=TRUE, verb=1)
msg(paste("Pruning query:\n"))
query <- segmentPrune(x=query,  chrL=chrL, remove.empty=TRUE, verb=1)

## sorting both
target <- segmentSort(target)
query <- segmentSort(query)

## PRUNING CHROMOSOMES
## only use number of chromosomes in the chromosome index file;
## this is e.g. useful to just skip the mitochondrial genome
nchr <- length(chrL)
rmq <- which(query$chr > nchr)
rmt <- which(target$chr > nchr)
if ( length(rmq)+length(rmt)>0 ) {
    msg(paste0("Removing undefined chromosomes:\n\tquery: ", length(rmq),
               "\n\ttarget: ", length(rmt), "\n"))
    query <- query[-rmq,]
    target <- target[-rmt,]
}

## merging if required by the processing steps
if ( merge | merget )  {
    msg(paste("Merging target:\n"))
    target <- segmentMerge(x=target, type=tclass, verb=1)
}
if ( merge | mergeq ) {
    msg(paste("Merging query:\n"))
    query <- segmentMerge(x=query, type=qclass, verb=1)
}

### TODO: setup run, generate mini-target
if ( setup ) {

    msg("SETUP RUN: only generating randomized queries and running a test\n")

    if ( random=="" ) 
        stop(paste("setup run makes no sense without permanently storing",
                   "randomized queries; use --random option."))
    
    target <- data.frame(chr=1, start=1, end=10, strand=1,
                         type="all", ID="test")
    tclass <- "type"
}


## CHECK FINAL SIZES
if ( verb>0 )
    msg(paste("TARGETS\t", nrow(target), "\n",
              "QUERIES\t", nrow(query), "\n",sep=""))

if ( nrow(query)==0 | nrow(target)==0 )
    stop("Empty query (",nrow(query),") or target (", nrow(target), ")")



if ( verb>0 )
    msg(paste("CALCULATE OVERLAPS\t",time(),"\n",sep=""))

## symmetric: only for antisense of self!
symmetric <- (antisense|convergent!=0) & self
if ( symmetric )
    msg(paste("\n\tNOTE: symmetric test of",
              "antisense|convergent with self!\n"))

## calculate Jaccard Index and permutation test
delete.data.message <- FALSE

if ( bedtools ) { ## OVERLAP via bedtools

    ## temporary file directory
    if ( random=="" ) 
        random <-  tempdir()
    else   delete.data.message <- TRUE
    if ( !dir.exists(random) )
        dir.create(random)
    
    ovl <- segmentOverlaps_bed(query=query, target=target, chrL=chrL,
                               prefix="socl_", symmetric=symmetric,
                               qclass=qclass, tclass=tclass, perm=perm, 
                               verb=1, tmpdir=random, save.permutations=TRUE,
                               runid=basename(outfile))
    
} else { ## OVERLAP via internal functions
    
    ## converting both to continuous index
    query  <- coor2index(query, chrS)
    target <- coor2index(target, chrS)
    
    ovl <- segmentOverlaps(query=query, target=target,
                           qclass=qclass, tclass=tclass, perm=perm, total=total,
                           symmetric=symmetric, verb=1)
    
    
    ## ADD COUNTS
    if ( count ) {
        tcol <- tclass
        qcol <- qclass
        if ( !"ID"%in%colnames(query) )
            query <- cbind(query, ID=1:nrow(query))
        qcol <- c("ID",qcol)
        qcol <- qcol[qcol!=""]
        
        if ( !"ID"%in%colnames(target) ) 
            target <- cbind(target, ID=1:nrow(target))
        tcol <- c("ID",tcol)
        tcol <- tcol[tcol!=""]
        
        ann <- annotateTarget(query=query, target=target,
                              collapse=FALSE,  details=FALSE, only.best=FALSE,
                              qcol=qcol, tcol=tcol, prefix="query")
        
       
        ## FILTER EMPTY
        ann <- ann[!is.na(ann[,paste0("query_ID")]),]
        
        ## count table
        ovl$count <- ovl$jaccard
        ovl$count[] <- 0
        tab <- as.matrix(table(ann[,paste0("query_",qclass)], ann[,tclass]))
        
        if ( symmetric ) {
            tab[upper.tri(tab)] <-
                tab[upper.tri(tab)] + t(tab)[upper.tri(tab)]
            tab[lower.tri(tab)] <- 0
        }
        
        ## empty dim names
        rstr <- paste(sample(c(letters, LETTERS),12, replace=TRUE),
                      collapse="")
        colnames(tab)[colnames(tab)==""] <- rstr
        colnames(ovl$count)[colnames(ovl$count)==""] <- rstr
        rownames(tab)[rownames(tab)==""] <- rstr
        rownames(ovl$count)[rownames(ovl$count)==""] <- rstr
        
        ## copy overlap count
        ovl$count[rownames(tab),colnames(tab)] <- tab
        
        ## reset empty dim names
        colnames(ovl$count)[colnames(ovl$count)==rstr] <- ""
        rownames(ovl$count)[rownames(ovl$count)==rstr] <- ""
        
        ## ADD TO RESULT STRUCTURE
        ovl$annotation <- ann
    }
}


if ( verb>0 )
    msg(paste0("DONE\t",time(),"\n"))

if ( setup ) {
    msg(paste0("SETUP DONE; check log files and re-use randomized",
               " queries in directory '",random,"'\n"))
    quit(save="no")
}

if ( verb>0 )
  msg(paste0("writing results\n"))

## store data

## store settings
parameters <- list()
parameters$permutations <- perm
parameters$bedtools <- bedtools
parameters$genomelength <- total
parameters$range <- range
parameters$upstream <- upstream
parameters$convergent <- convergent
parameters$antisense <- antisense
parameters$symmetric <- symmetric
ovl$parameters <- parameters

##OUTFILE NAME
if ( !interactive() ) {
    save(ovl, file=paste0(outfile,".rda"))
    if ( "annotation" %in% names(ovl) )
        write.table(ann, file=paste0(file.name,"_annotation.tsv"),
                    sep="\t", row.names=FALSE, quote=FALSE)
}
  
## plot
if ( perm>0 ) {

    ## only compare forward and reverse strands for auto-target antisense
    ## plot axis labels
    ## plot axis labels
    qlab <- paste0("query: ", qclass)
    tlab <- paste0("target: ", tclass)
    if ( (antisense|convergent!=0) & self ) {
        qlab <- paste0("query: ", qclass, ", strand ", frw.str[2])
        tlab <- paste0("target: ", tclass, ", strand ", rev.str[2])
    }
    if ( antisense & !self ) {
        tlab <- paste(tlab, "- antisense")
    }
    if ( upstream!=0 ) {
        qlab <- paste("query:", qclass)
        tlab <- paste("target:", tclass, "- up/downstream", upstream)
    }
    
    hbase <- 0.25
    wbase <- 1.5*hbase
    wd <- (ncol(ovl$p.value))*wbase + (.75+.6)
    ht <- (nrow(ovl$p.value))*hbase + (.75+.6)
    
    plotdev(paste0(outfile),type=fig.type,width=wd,  height=ht)
    par(mai=c(.75,.75,.6,.6),mgp=c(1.75,.3,0),tcl=-0.1)
    plotOverlaps(ovl,p.min=.001, values="count",
                 show.total=TRUE, short=TRUE, ylab=qlab, xlab=tlab)
    dev.off()
}
if ( delete.data.message )
    warning(paste0("don't forget to delete randomized sequence files ",
                   "in directory '", random, "'"))
