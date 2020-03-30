#!/usr/bin/env Rscript

## SEGMENT CLASSES OVERLAP STATISTIC
## by permutation test analysis

## TODO: implement downstream scan

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
  make_option(c("--intersegment"), type="character", default="", 
              help="use this value as additional 'qclass' for inter-segment stretches"),
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
  make_option(c("--antisense"), action="store_true", default=FALSE,
              help="search target matches on reverse strand (if target is empty; search will be done for sense query vs. antisense query!"),
  make_option(c("--upstream"), type="integer", default=0,
              help="search range upstream of target (in nt.)"),
##  make_option(c("--downstream"), type="integer", default=0,
##              help="search range downstream of target (in nt.)"),
  make_option(c("--perm"), type="integer", default=100, 
              help="number of permutations"),
  make_option(c("--count"),  action="store_true", default=FALSE,
              help="count individual overlaps"),
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
if ( antisense & upstream!=0 )
    stop("options --antisense and --upstream are incompatible!")
if ( nostrand & antisense )
    stop("options --nostrand and --antisense are incompatible")
if ( nostrand & upstream!=0 )
    stop("options --nostrand and --upstream are incompatible")


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

## consider only one strand
if ( nostrand ) {
    total <- total/2
    query$strand <- "+"
    target$strand <- "+"
}

## converting both to continuous index
query <- coor2index(query, chrS)
target <- coor2index(target, chrS)

## search on other strand 
if ( antisense ) 
    target <- switchStrand(target, chrS)

## upstream bug yeast - missing 2 micron plasmid - catch here
if ( any(is.na(query[,"start"])) ) {
    rm <- which(is.na(query[,"start"]))
    warning("removing ", length(rm), " queries due to missing coordinates",
            "(in yeast: upstream error, missing chr column for 2-micron genes)")
    query <- query[-rm,]
}
if ( any(is.na(target[,"start"])) ) {
    rm <- which(is.na(target[,"start"]))
    warning("removing ", length(rm), " queries due to missing coordinates",
            "(in yeast: upstream error, missing chr column for 2-micron genes)")
    target <- target[-rm,]
}
    

if ( verb>0 )
    msg(paste("CALCULATE OVERLAPS\t",time(),"\n",sep=""))

## add inter-segments here
if ( intersegment!="" ) {

    ## convert to numeric, if qclass is numeric
    
    if ( is.numeric(query[,qclass]) )
        intersegment <- as.numeric(intersegment)

    maxL <- ifelse(nostrand, max(chrS), 2*max(chrS))
 
    emptyseg <- data.frame(ID=paste("is",sprintf("%04d",1:(nrow(query)+1)),
                                    sep=""),
                           chr=rep(1,nrow(query)+1),
                           start=c(1,query[1:nrow(query),"end"]+1),
                           end=c(query[1:nrow(query),"start"]-1,maxL),
                           stringsAsFactors=FALSE)
    emlen <- emptyseg[,"end"] - emptyseg[,"start"] +1 
    emptyseg <- emptyseg[emlen>0,] # rm 0-length and negative overlaps
    ## split chromosome end segments
    emptyseg <- splitsegs(emptyseg,chrS,idcol="ID",verb=verb)

    emq <- as.data.frame(matrix(NA, ncol=ncol(query), nrow=nrow(emptyseg)))
    colnames(emq) <- colnames(query)
    emptyseg <- emptyseg[,colnames(emptyseg)%in%colnames(emq)]
    emq[,colnames(emptyseg)] <- emptyseg
    if ( qclass%in%colnames(emq) )
        emq[,qclass] <- intersegment ## add intersegment cluster
    if ("ID"%in%colnames(query) )
        query <- rbind(query, emq)
}
## calculate Jaccard Index and permutation test
ovl <- segmentJaccard(query=query, target=target,
                      qclass=qclass, tclass=tclass, perm=perm, total=total,
                      verb=1)


## ADD COUNTS
if ( count ) {
    tcol <- tclass
    qcol <- qclass
    if ( "ID"%in%colnames(query) )
        qcol <- c("ID",qcol)
    if ( "ID"%in%colnames(target) ) {
        tcol <- c("ID",tcol)
        tcol <- tcol[tcol!=""]
    }
    ann <- annotateTarget(query=query, target=target,
                          collapse=FALSE,  details=FALSE, only.best=FALSE,
                          qcol=qcol, tcol=tcol, prefix="query")
    
    if ( tclass=="" ) {
        ann <- cbind(ann, tclass="target")
        tcol <- tclass <- "tclass"
    }
    if ( qclass=="" ) {
        ann <- cbind(ann, qclass="query")
        qcol <- qclass <- "qclass"
    }

    ## FILTER EMPTY
    ann <- ann[!is.na(ann[,paste0("query_",qclass)]),]
    
    ## count table
    ovl$count <- ovl$p.value
    ovl$count[] <- 0
    tab <- as.matrix(table(ann[,paste0("query_",qclass)], ann[,tclass]))

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


if ( verb>0 )
  msg(paste0("DONE\t",time(),"\n"))
if ( verb>0 )
  msg(paste0("writing results\n"))

## write out results
if ( tclass=="" ) tclass <- "all"
if ( qclass=="" ) qclass <- "all"

## store data
##OUTFILE NAME
file.name <- paste0(outfile,"_",qclass,"_",tclass,
                    ifelse(antisense,"_antisense",""),
                    ifelse(upstream!=0, paste0("_upstream",upstream),""))
if ( !interactive() ) {
    save(ovl, file=paste0(file.name,".RData"))
    if ( "annotation" %in% names(ovl) )
        write.table(ann, file=paste0(file.name,"_annotation.tsv"),
                    sep="\t", row.names=FALSE, quote=FALSE)
}
  
## plot
if ( perm>0 ) {
    plotdev(paste0(file.name),type=fig.type)
    plotOverlaps(ovl,p.min=.001,main="Jaccard Index (*1000) & permutation test",ylab=qlab,xlab=tlab,scale=1000,round=0)
    dev.off()
}
