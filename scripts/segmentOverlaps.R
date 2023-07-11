#!/usr/bin/env Rscript

## SEGMENT CLASSES OVERLAP STATISTIC
## by permutation test analysis

## TODO: test downstream implementation, just using negative upstream

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
              help="search upstream (>0) or downstream (<0) of target (in bp)"),
  make_option(c("--range"), type="integer", default=0,
              help="search range around target (in bp)"),
  make_option(c("--convergent"), type="integer", default=0,
              help="search range downstream of query and target (in bp)"),
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

## plot axis labels
qlab <- paste0("query: ", qclass)
tlab <- paste0("target: ", tclass)

## target = antisense of query?
frw.str <- c("1","+")
rev.str <- c("-1","-")
## comparison with self on reverse strand!
## compare forward with reverse strand!
self <- FALSE
if ( target=="" & (antisense|convergent!=0) ) {
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
    ## TODO: map chromosome name to index
    ## perhaps allow to pass chromosome map to coor2index
    
    ## FILTER targets
    if ( length(ttypes)>0 )
        if( ttypes!="" )
            target <- target[target[,ttypcol]%in%ttypes,]
}

## scan for range around targets
if ( range!=0 ) {

    ## make sure start < end for both strands
    ## (only valid for non-circular DNA!!)
    utarget <- target[,c("start","end","strand")]
    str <- as.character(utarget[,"strand"])
    start <- utarget[,"start"]
    end <- utarget[,"end"]
    
    ## get start and end: independent of start/end definition; using strand info
    start[str%in%rev.str] <-
        apply(utarget[str%in%rev.str,c("start","end")],1,min)
    end[str%in%rev.str] <- apply(utarget[str%in%rev.str,c("start","end")],1,max)

    utarget[str%in%frw.str, c("start","end")] <-
        cbind(start-range, end+range)[str%in%frw.str,]
    
    utarget[str%in%rev.str, c("start","end")] <-
        cbind(start+range, end-range)[str%in%rev.str,]
        
    target[,c("start","end","strand")] <- utarget
}
## scan for upstream (>0) or downstream (<0) of target
if ( upstream!=0 ) {

    ## make sure start < end for both strands
    ## (only valid for non-circular DNA!!)
    utarget <- target[,c("start","end","strand")]
    str <- as.character(utarget[,"strand"])
    start <- utarget[,"start"]
    end <- utarget[,"end"]
    
    ## get start and end: independent of start/end definition; using strand info
    start[str%in%rev.str] <-
        apply(utarget[str%in%rev.str,c("start","end")],1,min)
    end[str%in%rev.str] <- apply(utarget[str%in%rev.str,c("start","end")],1,max)

    if ( upstream > 0 ) {
        
        utarget[str%in%frw.str, c("start","end")] <-
            cbind(start-upstream, start-1)[str%in%frw.str,]

        utarget[str%in%rev.str, c("start","end")] <-
            cbind(end+1, end+upstream)[str%in%rev.str,]

    } else { # use negative upstream as downstream from end!

        utarget[str%in%frw.str, c("start","end")] <-
            cbind(end+1, end-upstream)[str%in%frw.str,]

        utarget[str%in%rev.str, c("start","end")] <-
            cbind(start+upstream,  start-1)[str%in%rev.str,]
    
    }
    target[,c("start","end","strand")] <- utarget
}

if ( convergent!=0 ) {
    ## same as upstream<0 but for both query and target
    ## forward strand: max(start,end) + range
    ## reverse strand: min(start,end) - range
    for ( tmp in c("query","target") ) {
        
        ## make sure start < end for both strands
        ## (only valid for non-circular DNA!!)
        udat <- get(tmp)[,c("start","end","strand")]
        str <- as.character(udat[,"strand"])
        start <- udat[,"start"]
        end <- udat[,"end"]
        
        ## get start and end: independent of start/end definition; using strand info
        start[str%in%rev.str] <-
            apply(udat[str%in%rev.str,c("start","end")],1,min)
        end[str%in%rev.str] <- apply(udat[str%in%rev.str,c("start","end")],1,max)
        
        udat[str%in%frw.str, c("start","end")] <-
            cbind(end+1, end+convergent)[str%in%frw.str,]

        udat[str%in%rev.str, c("start","end")] <-
            cbind(start-convergent,  start-1)[str%in%rev.str,]
        assign(x=paste0("u",tmp), value=udat)
        
    }
    target[,c("start","end","strand")] <- utarget
    query[,c("start","end","strand")] <- uquery
}

## only compare forward and reverse strands for auto-target antisense
## plot axis labels
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

## cut chromosome ends:
## upstream and convergent may have produced non-existent coordinates!
target <- pruneSegments(target, chrL=chrL, remove.empty=TRUE, verb=1)
query  <- pruneSegments(query,  chrL=chrL, remove.empty=TRUE, verb=1)

## converting both to continuous index
query  <- coor2index(query, chrS)
target <- coor2index(target, chrS)

## search antisense and convergent targets on other strand
## for self: targets are reverse strand features, see above
if ( (antisense|convergent!=0) ) 
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

## symmetric: only for antisense of self!
symmetric <- (antisense|convergent!=0) & self
if ( symmetric )
    cat(paste("\n\tNOTE: symmetric test of antisense|convergent with self!\n"))

## TODO: replace this (optionally) with segmentoverlaps_bed.sh
## write out bed files here (before coor2index), using coor2bed,
## and construct script call, and parse results.
ovl <- segmentJaccard(query=query, target=target,
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
    
    if ( tclass=="" ) {
        tcol <- tclass <- "tclass"
        ann <- cbind(ann, tclass="target")
    }
    if ( qclass=="" ) {
        qcol <- qclass <- "qclass"
        ann <- cbind(ann, query_qclass="query")
    }

    ## FILTER EMPTY
    ann <- ann[!is.na(ann[,paste0("query_ID")]),]
    
    ## count table
    ovl$count <- ovl$jaccard
    ovl$count[] <- 0
    tab <- as.matrix(table(ann[,paste0("query_",qclass)], ann[,tclass]))

    if ( symmetric ) {
        tab[upper.tri(tab)] <-
            tab[upper.tri(tab)] + tab[lower.tri(tab)]
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


if ( verb>0 )
  msg(paste0("DONE\t",time(),"\n"))
if ( verb>0 )
  msg(paste0("writing results\n"))

## write out results
if ( tclass%in%c("","tclass") ) tclass <- "all"
if ( qclass%in%c("","qclass") ) qclass <- "all"

## store data
##OUTFILE NAME

streamid <- ifelse(upstream>0, "_upstream","_downstream")
file.name <- paste0(outfile,"_",qclass,"_",tclass,
                    ifelse(antisense,"_antisense",""),
                    ifelse(upstream!=0, paste0(streamid,upstream),""))

## store settings
parameters <- list()
parameters$permutations <- perm
parameters$genomelength <- total
parameters$range <- range
parameters$upstream <- upstream
parameters$convergent <- convergent
parameters$antisense <- antisense
parameters$symmetric <- symmetric
ovl$parameters <- parameters

if ( !interactive() ) {
    save(ovl, file=paste0(file.name,".rda"))
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
