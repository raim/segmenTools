#!/usr/bin/env Rscript

## SEGMENT CLASSES OVERLAP STATISTIC
## by permutation test analysis


## interactive debug: command-line call with parameters
## R --args  -q /home/raim/work/yeastSeq2016/data/segmentation/segmentTest/20170307_merge/analysis4/D.dft1-7.dcash.snr_T.raw_K.12_S.icor_E.3_M.200_nui.3/D.dft1-7.dcash.snr_T.raw_K.12_S.icor_E.3_M.200_nui.3_annotated.csv --qclass mCL --tclass mCL --antisense --chrfile /home/raim/programs/tataProject/yeast/chromosomes/sequenceIndex_R64-1-1_20110208.csv -o test.Rdata
## R --args -q /home/raim/work/yeastSeq2016/data/segmentation/segmentTest/20170307_merge/analysis4/D.dft1-7.dcash.snr_T.raw_K.12_S.icor_E.3_M.200_nui.3/D.dft1-7.dcash.snr_T.raw_K.12_S.icor_E.3_M.200_nui.3_annotated.csv --qclass mCL -t /home/raim/programs/tataProject/yeast/feature_R64-1-1_20110208_withclusters.csv --tclass CL_rdx --antisense --chrfile /home/raim/programs/tataProject/yeast/chromosomes/sequenceIndex_R64-1-1_20110208.csv -o test.RData

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
  make_option(c("--perm"), type="integer", default=100, 
              help="number of permutations"),
  ## OUTPUT
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
cf <- read.table(chrfile,sep="\t",header=FALSE)
chrL <- cf[,3]
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos

## READ SEGMENTS TO BE TESTED 
if ( verb>0 ) msg(paste("Loading query:", query, "\t\n"))
query <- read.delim(query, stringsAsFactors=FALSE)

if ( verb>0 ) msg(paste("Loading target:", target, "\t\n"))

## target = antisense of query?
samesame <- FALSE # only do forward or reverse strand if testing against self?
frw.str <- 1
rev.str <- -1
if ( target=="" & antisense ) {
    target <- query
    samesame <- TRUE
    ## only compare forward and reverse strands for auto-target
    query <- query[query[,"strand"]==frw.str,]
    target <- target[target[,"strand"]==rev.str,]
    tclass <- qclass
    ## plot axis labels
    qlab <- paste0("query: ", qclass, ", strand", frw.str)
    tlab <- paste0("target: ", tclass, ", strand", rev.str)
} else {
    target <- read.delim(target, stringsAsFactors=FALSE)
    ## FILTER targets
    if ( ttypes!="" )
      target <- target[target[,ttypcol]%in%ttypes,]

    ## plot axis labels
    qlab <- paste0("query: ", qclass)
    tlab <- paste0("target: ", tclass)
   
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

getJaccard <- function(query, target, qclass, tclass) {

    ## query classes
    qcls <- as.factor(query[,qclass])
    qcls.srt <- sort(unique(qcls))
    qN <- length(qcls.srt)
    
    ## target classes
    tcls <- as.factor(target[,tclass])
    tcls.srt <- sort(unique(tcls))
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
    colnames(jc) <- tcls.srt
    rownames(jc) <- qcls.srt
    for ( i in 1:qN ) {
        for ( j in 1:tN ) { ## asymmetric only for target = antisense of query?
            is <- length(intersect(qcls.rng[[i]],tcls.rng[[j]]))
            un <- length(union(qcls.rng[[i]],tcls.rng[[j]]))
            jc[i,j] <- is/un
        }
    }
    jc
}

J.real <- getJaccard(query, target, qclass, tclass)

## randomize queries

J.pval <- J.real
J.pval[] <- 0
for ( i in 1:perm ) {

    cat(paste(i/perm," "))
    
    ## sort
    query <- query[order(query$start),]
    
    ## segment lengths
    sglen <- apply(query[,c("start","end")], 1, diff)
    sglen <- sglen+1
    sgcls <- as.factor(query[,qclass])

    ## inter-segment lengths
    qnum <- nrow(query)
    islen <- apply(cbind(query[1:(qnum-1),"end"],query[2:qnum,"start"]), 1, diff)
    islen <- c(query[1,"start"], islen, 2*sum(chrL)-query[qnum,"end"]+1)
    islen <- islen-1
    
    ## check (if start and end were added
    if ( sum(islen)+sum(sglen) != 2*sum(chrL) )
        stop()
    
    ## sample segment lengths
    ridx <- sample(1:qnum)
    rcls <- sgcls[ridx] # store cluster to keep cluster length dist!
    rsglen <- sglen[ridx]
    ## sample intersegment lengths
    rislen <- sample(islen)
    
    ## construct randomized segmentation
    tot <- length(rsglen)+length(rislen)
    cumlen <- rep(NA,tot)
    cumlen[seq(1,tot,2)] <- rislen
    cumlen[seq(2,tot,2)] <- rsglen
    cumlen <- cumsum(cumlen)
    
    rquery <- data.frame(start=cumlen[seq(1,tot-1,2)],
                         end=cumlen[seq(2,tot-1,2)],
                         type=rcls)
    
    ## TODO: test cluster length distribution?

    J.rnd <- getJaccard(rquery, target, qclass="type", tclass)
    J.pval <- J.pval + as.numeric(J.rnd >= J.real)
}
cat(paste("\n"))

## p-value
J.pval <- J.pval/perm


## plot
ovl <- list()
ovl$overlap <- round(1000*J.real)
ovl$p.value <- J.pval

pdf(paste0(sub(".RData","",basename(outfile)),"_",qclass,"_",tclass,".pdf"))
plotOverlaps(ovl,p.min=.001,main="Jaccard Index (*1000) & permutation test",ylab=qlab,xlab=tlab)
dev.off()

## store
save.image(file=outfile)
