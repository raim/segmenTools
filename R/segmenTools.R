#' segmenTools : analysis of genome segmentations by segmenTier
#'@author Rainer Machne \email{raim@tbi.univie.ac.at}
#'@docType package
#'@name segmenTools
#'@section Dependencies: basic (\code{stats}, \code{graphics}, \code{grDevices}), clustering, \code{flowClust}, \code{flowMerge}
#'@importFrom graphics image axis par plot matplot points lines legend arrows strheight strwidth text abline hist
#'@importFrom grDevices png dev.off rainbow gray xy.coords
#'@importFrom stats ecdf loess predict qt quantile runmed sd var
NULL # this just ends the global package documentation


### GENOME SEGMENT UTILS & ANALYSIS

### DATA STAT & TRANSFORMATION UTILS

#' perform Discrete Fourier Transformation using \code{\link[stats]{mvfft}},
#' and returning the non-redundant (for real numbers) first half of
#' the transform, i.e., from the DC (direct current) component to the
#' Nyquist frequency
#' @param x data to be transformed
#' @export
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n]
    colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}

#' asinh trafo, an alternative to log transformation that has less
#' (compressing) effects on the extreme values (low and high values),
#' and naturally handles negative numbers and 0
#' @param x data to be transformed
#' @export
ash <- function(x) log(x+sqrt(x^2+1))
#' log trafo handling zeros by adding 1
#' @param x data to be transformed
#' @return log(x+1)
#' @export
log_1 <- function(x) log(x+1)

#' moving average using \code{\link[stats]{filter}}
#' @param x data vector along which a moving average will be calculated
#' @param n moving average window size
#' @param circular logical see help of function \code{\link[stats]{filter}}
#' @export
ma <- function(x, n=5, circular=FALSE) {
    stats::filter(x,rep(1/n,n), sides=2, circular=circular)
}

#' calculate 95% confidence intervals for the given
#' data vector using a t-distribution
#' @param x data vector for which the 95% confidence interval
#' will be calculated
#' @param na.rm remove NA values before calculation
#' @export
ci95 <- function(x,na.rm=FALSE) {
    if ( na.rm ) x <- x[!is.na(x)]
    n <- length(x)
    if ( n<2 ) return(NA)
    error <- qt(0.975, df=n-1) * sd(x)/sqrt(n)
    return(error)
}



### PLOT UTILS
#' Switch between plot devices
#' @param file.name file name without suffix (.png, etc)
#' @param type plot type: pdf, png or eps
#' @param width figure width in inches
#' @param height figure height in inches
#' @param res resolution in ppi (pixels per inch), only for 'png'
#' @export
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100) {
  file.name <- paste(file.name, type, sep=".")
  if ( type == "png" )
    grDevices::png(file.name, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    grDevices::postscript(file.name, width=width, height=height,paper="special")
  if ( type == "pdf" )
    grDevices::pdf(file.name, width=width, height=height)
}

#' Some plotting routines for overlap statistics from
#' \code{\link{getOverlapStats}}
#' @export
plot.cdfLst <- function(x=seq(0,2,.05), CDF, type="rcdf", col, lty, h=c(.2,.8), v=c(0.8,1.2), ylab="cum.dist.func.", ...) {
    
    plot(x, CDF[[1]]$rcdf(x),type="l",col=NA,main=NA, ylim=c(0,1), ylab=ylab,
         ...)
    abline(v=v,lty=2)
    abline(h=h,lty=2)
    abline(h=0:1, lty=2, col="gray",lwd=.75)
    for ( i in 1:length(CDF) ) 
        if ( type %in% names(CDF[[i]]) ) 
            lines(x,CDF[[i]][[type]](x),col=col[i],lty=lty[i]) 
}

#' Wrapper around \code{\link[graphics]{image}} to plot a matrix as
#' it is displayed in R console, i.e. the field \code{dat[1,1]}
#' is at the top left corner. It further allows to plot text
#' into individual fields and have colored axis tick labels.
#' @param dat the numeric data matrix to be plotted
#' @param text a matrix of characteres corresponding to \code{dat}
#' which will be plotted on the image
#' @param axis integer vector, sets whether bottom (1) and/or left
#' (2) axis are draw; the column and row names of \code{dat} will
#' be used as tick labels
#' @param axis1.col invididual colors for x-axis tick labels, length must
#' equal the number of columns of \code{dat}
#' @param axis2.col invididual colors for y-axis tick labels, length must
#' equal the number of rows of \code{dat}
#' @param ... further arguments to \code{\link[graphics]{image}}, e.g., col
#' to select colors
#' @export
image_matrix <- function(dat, text, axis, axis1.col, axis2.col, ...) {

    ## reverse columns and transpose
    imgdat <- t(apply(dat, 2, rev))
    image(x=1:ncol(dat), y=1:nrow(dat), z=imgdat , ...)

    ## add text
    if ( !missing(text) )
        text(x=rep(1:ncol(dat),nrow(dat)), y=rep(nrow(dat):1,each=ncol(dat)),
             paste(t(text)))

    ## add axes
    if ( !missing(axis) ) {
        if ( 1 %in% axis ) 
            if ( !missing(axis1.col) ) # colored ticks
                for ( i in 1:nrow(dat) )
                    axis(1, at=i,colnames(dat)[i],
                         col.axis=axis1.col[i], col=axis1.col[i],
                         las=2,cex.axis=1.5, lwd=4)
            else
                axis(1, at=1:ncol(dat), colnames(dat),las=2)
        if ( 2 %in% axis )
            if ( !missing(axis2.col) ) # colored ticks
                for ( i in 1:nrow(dat) )
                        axis(2, at=nrow(dat)-i+1,rownames(dat)[i],
                             col.axis=axis2.col[i], col=axis2.col[i],
                             las=2,cex.axis=1.5, lwd=4)
            else
                axis(2, at=nrow(dat):1, rownames(dat),las=2)        
    }
} 


### SEGMENT UTILS

#' which segment covers a position
#' @param pos chromosome position in continuous index
#' @param seg segment (genomic interval)
#' @export
whichSegment <- function(pos, seg) 
    which(seg[,"start"]<= pos & seg[,"end"]>=pos)

## calculate overlap between two sets of genome segments, e.g.
## test a segmentation of RNA-seq data vs. known features or transcripts
## NOTE: all chromosome coordinates must be mapped to a continuous index
## via coor2index
## NOTE: mod from $TATADIR/yeast/scripts/analyzeSeq2013_utils.R

#' specifically tailored to segment strings in segmenTier;
#' splits a list of strings by a separator and constructs
#' a table from the entries.
#' @param sgtypes a list of strings, that is converted to a table
#' of classes based on a string separator (\code{sep})
#' @param sep the separator for classes in the string
#' @param gsep a separator within classes
#' to be used as classification ID (column header of the produces
#' class table)
#' @export
getSegmentClassTable <- function(sgtypes, sep="_", gsep=":") {
    ## get classes
    cllst <- strsplit(sgtypes,"_")
    if ( length(unique(unlist(lapply(cllst, length))))>1 )
        cat(paste("WARNING: inconsistent segment classes\n"))
    ## get class ids
    clid <- unique(unlist(lapply(cllst, function(x) sub(":.*","",x))))

    ## fill data.frame
    cltab <- data.frame(matrix(NA, nrow=length(sgtypes), ncol=length(clid)))
    colnames(cltab) <- clid
    rownames(cltab) <- sgtypes
    for ( i in 1:length(cllst) ) {
        tmp <- strsplit(cllst[[i]],":")
        class <- unlist(lapply(tmp, function(x) x[2]))
        names(class) <- unlist(lapply(tmp, function(x) x[1]))
        cltab[i, names(class)] <- class
    }
    ## convert numeric to numeric
    for ( i in 1:ncol(cltab) ) {
        col <- cltab[,i]
        if ( suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
            cltab[,i] <- as.numeric(col)
        } else {
            cltab[,i] <- as.factor(col)
        }
    }
    
    cltab
}

#' specifically tailored to segment strings in segmenTier
#' @param sgtypes a list of strings, that is converted to a table
#' of classes based on a string separator (\code{sep})
#' @param sep the separator for classes in the string
#' @param gsep currently not used, a separator within classes
#' to be used as classification ID (column header of the produces
#' class table)
#' @export
getSegmentClasses <- function(sgtypes, sep="_", gsep=":") {
    sgclasses <- sort(unique(unlist(strsplit(sgtypes,"_"))))
    ## filter classes that are not varied
    rm <- NULL
    for ( class in sgclasses ) 
        if ( length(grep(class,sgtypes))==length(sgtypes) )
            rm <- c(rm, which(sgclasses==class))
    if ( !is.null(rm) )
        sgclasses <- sgclasses[-rm]
    sgclasses
}

### WRAPPERS of segmentOverlap for specific purposes

## wrapper around \code{\link{segmentOverlap}}, used to 
## annotate the query set by a column in the target set
## @param query the query set of segments (genomic intervals)
## @param target the target set of segments (genomic intervals)
## @param col column names to copy from target to it's best match in query
annotateQuery <- function(query, target, col) {
    ## TODO: fix this, make real annotate query
}

#' wrapper around \code{\link{segmentOverlap}}, used to 
#' annotate the target set by a column in the query set
#' @param query the query set of segments (genomic intervals)
#' @param target the target set of segments (genomic intervals)
#' @param col column names to copy from target to it's best match in query
#' @param prefix column name prefix for the copied columns
#' @param details set to \code{TRUE} to add details of the used
#' match, i.e., union and intersect, and relative position of the query
#' to the matching targets
#' @param duplicates how to handle multiple matches; if "collapse" the
#' multiple hits are collapsed into ;-separated strings
#' @export
annotateTarget <- function(query, target, col, prefix, details=FALSE,
                           duplicates="collapse") {

    ## TODO: use details flag to also bind details of overlap (left/right)
    #cltr <- annotateQuery(query, target, col)
    cltr <- segmentOverlap(query=query,target=target,
                           add.na=TRUE,details=details,sort=FALSE)

    ## bind query column to overlap table
    cltr <- cbind(cltr, query[cltr[,"query"],col,drop=FALSE])

    ## reduce to top-ranking query hit
    best <- cltr[cltr[,"qrank"]==1,,drop=FALSE]
    
    ## collapse or remove duplicated with qrank==1
    if ( sum(duplicated(best[,"target"])) )
        cat(paste("handling", sum(duplicated(best[,"target"])), "duplicated:",
                  duplicates, "\n"))
    dups <- which(duplicated(best[,"target"]))
    if ( duplicates=="collapse" ) {
        for ( d in dups ) {
            id <- best[d,"target"]
            idx <- which(best[,"target"]==id) ## find all 
            ovl <- best[idx,]
            covl <- rep(NA, ncol(ovl))
            names(covl) <- colnames(ovl)
            for ( j in 1:ncol(ovl) )
                covl[j] <- paste(ovl[,j],collapse=";")
            covl["target"] <- ovl[1,"target"]
            covl["tlen"] <- ovl[1,"tlen"]
            covl["qrank"] <- ovl[1,"qrank"]
            best[idx[1],] <- covl
        }
    }
    if ( length(dups)>0 )
        best <- best[-dups,]

    ## e.g., details = c("qpos", "intersect", "union")
    if ( details )
        col <- c(col, "qpos", "qlen", "intersect", "union")
    
    ## get and optionally rename requested columns
    idx <- match(best[,"target"],1:nrow(target))
    addcol <- best[idx,col, drop=FALSE]

    ## add prefix to column ids
    if ( !missing(prefix) )
        if ( prefix!="" )
            colnames(addcol) <- paste(prefix,"_",colnames(addcol),sep="")

    ## bind to target
    target <- cbind(target, addcol)

    ## add index to column ID if already present
    ## TODO: this is from clustering in segmenTier, adapt
    ## and make smarter: count existing indexes
    ##if ( any(duplicated(colnames(target))) ) {
    ##    sel <- paste(colnames(target),".1",sep="")
    ##    cnt <- 2
    ##    while( sum(duplicated(sel)) ) {
    ##        sel[duplicated(sel)] <- sub("\\..*",paste(".",cnt,sep=""),
    ##                                    sel[duplicated(sel)])
    ##       cnt <- cnt+1
    ##    }
    ##    colnames(target) <- sub("\\.1$","",sel)
    ##}

    target
}

#' Collect stats from a nested list of overlap statistics in lists
#' and vectors for a given segment \code{type}
#' @param ovlStatLst list of overlap statistics, where individual
#' entries come from \code{\link{segmentOverlap}}
#' @param type name of the overlap statistics list first level
#' @export
collectOvlStats <- function(ovlStatLst, type) {
    CDF <- lapply(ovlStatLst[[type]], function(x) x$CDF)
    class(CDF) <- "cdfLst"
    DIST <- lapply(ovlStatLst[[type]], function(x) x$DIST)
    height <- matrix(unlist(lapply(ovlStatLst[[type]],
                                   function(x) x$height)),ncol=2,byrow=TRUE)
    jaccard <- unlist(lapply(ovlStatLst[[type]],
                             function(x) x$jaccard))
    hitnum <- unlist(lapply(ovlStatLst[[type]],
                                function(x) x$hitnum))
    numhit <- unlist(lapply(ovlStatLst[[type]],
                            function(x) x$numhit))
    qnum <- unlist(lapply(ovlStatLst[[type]],
                          function(x) x$qnum))
    tnum <- unique(unlist(lapply(ovlStatLst[[type]],
                                 function(x) x$tnum)))
    nms <- names(CDF)
    names(jaccard) <- names(hitnum) <- names(numhit) <-
        names(qnum) <- names(CDF) <- names(DIST) <- NULL
    if ( length(tnum)>1 ) # temporary until sure that it works
        stop("ERROR,",type,"tnum should be unique for types")
    res <- list(CDF=CDF, DIST=DIST,
                jaccard=jaccard, height=height, hitnum=hitnum, numhit=numhit,
                qnum=qnum, tnum=tnum, nms=nms)
    class(res) <- "overlapStatLst"
    res
}

## TODO
plotOverlap <- function(ovlstats,type="rcdf",file.name) {
    if ( !missing(file.name) )
        png(file.name,width=400,height=400)
    ## TODO
   if ( !missing(file.name) ) dev.off()   
}

### MAIN ALGORITHM - FINDS OVERLAPS BETWEEN TWO SETS
### CHROMOSOMAL SEGMENTS

#' Simple sweeping algorithm to find overlapping intervals. Both
#' intervals must be in continuous index (coor2index) and start<end,
#' with strand information implied in the continuous index position.
#' see \code{bedtools} \code{intersects}, the Binary Interval
#' Search (BITS) algorithm, and the \code{NClist} algorithm (implemented
#' in packages IRanges and GenomicRanges) for similar tools.
#' It loops over targets and collects all queries that match;
#' optionally adds NA lines for non-matched targets
#' @param query the query set of segments (genomic intervals)
#' @param target the target set of segments (genomic intervals)
#' @param details add details on the relative positions of the query
#' wrt the target (covers,left,inside,right); note, that  keywords 'left'
#' and 'right' will be converted to 'upstream' or 'downstream' by
#' index2coor, accounting for strand info.
#' @param add.na add an empty overlap line for targets without overlaps
#' @param untie if several queries have equal top rank of 1, the rank
#' can be  replaced by simple order, such that only the first query
#' will have rank 1
#' @param collapse if \code{TRUE} multiple query hits are collapsed
#' into a single row, with ;-separated strings in the respective fields.
#' If both \code{add.na} and \code{collapse} are set to TRUE the resulting
#' overlap matrix will be of the same dimension as the input target, i.e.,
#' it will contain one row for each target.
#' @param sort sort query hits by their \code{rank}
#' @export
segmentOverlap <- function(query, target, details=FALSE, add.na=FALSE, untie=FALSE, collapse=FALSE, sort=FALSE) {

    ## get target and query ID - only required for messages
    if ( "ID" %in% colnames(target) ) {
        tids <- as.character(target[,"ID"])
    } else if ( !is.null(rownames(target)) ){
        tids <- rownames(target)
    } else {
        tids <- 1:nrow(target)
    }
    if ( "ID" %in% colnames(query) ) {
        qids <- as.character(query[,"ID"])
    } else if ( !is.null(rownames(query)) ){
        qids <- rownames(query)
    } else {
        qids <- 1:nrow(query)
    }

    overlaps <- NULL
    for ( k in 1:nrow(target) ) {

        ## overlapping 
        idx <- which(query[,"end"] >= target[k,"start"])
        idx <- idx[which(query[idx,"start"] <= target[k,"end"])]
      
        ## detailed overlap structure
        ## TODO: rm redundant code above
        ## TODO: allow upstream and downstream ranges!
        if ( details ) {
            ## 1: query ends right of target start (ALL overlaps)
            RhLt <- query[,"end"]   >= target[k,"start"]
            ## 2: query starts left of target end (ALL overlaps)
            LlRt <- query[,"start"] <= target[k,"end"]
            ## 3: query starts right of target start (INSIDE or RIGHT)
            LhLt <- query[,"start"] >= target[k,"start"]
            ## 4: query ends left of target end (INSIDE or LEFT)
            RlRt <- query[,"end"]   <= target[k,"end"]
            ## 5: query ends right of target end (INSIDE, COVERS or RIGHT)
            ##    ~ !RlRt, but incl. end
            RhRt <- query[,"end"]   >= target[k,"end"]  
            ## 6: query starts left of target start (INSIDE, COVERS or LEFT)
            ##    ~ !LhLt, but incl. start
            LlLt <- query[,"start"] <= target[k,"start"]
            
            ## combinations:
            overlp <- LlRt & RhLt # general overlap: should contain all below
            inside <- LhLt & RlRt # is inside target
            covers <- LlLt & RhRt # covers target
            inside[inside & covers] <- FALSE ## EXACT HIT
            left   <- !LhLt & !RhRt & RhLt # left overlapping
            right  <- !RlRt & !LlLt & LlRt # right overlapping

            ## debugging: these should NOT happen
            if ( sum(overlp,na.rm=T)!=length(idx) )
                stop(k,"ERROR, overlap vs. idx")
            if ( sum(overlp,na.rm=T) !=
                sum((covers|inside)|(left|right),na.rm=T) )
                stop(k,"ERROR, overlap vs. detail")
            if ( any(rowSums(cbind(inside,covers,right,left),na.rm=T)>1) )
                stop(k,"ERROR, too many hits")
        }

        ## add empty target if requested
        if ( length(idx)==0 ) {
            if ( add.na ) {
                result<-c(target=k,query=NA,intersect=0, union=0, qlen=0,tlen=0)
                if ( details ) result <- c(result,qrank=1,qpos=NA)
                else result <- c(result,qrank=1)
                overlaps <- rbind(overlaps,result)
            }
            next
        }
        ## calculate overlap %
        tgrng <- target[k,"start"]:target[k,"end"]
        tlen <- length(tgrng)
        ovl <- NULL
        tps <- NULL
        for ( id in idx ) {
            qurng <- query[id,"start"]:query[id,"end"]
            qlen <- length(qurng)
            intersect <- length(intersect(tgrng,qurng)) 
            union <- diff(range(c(qurng,tgrng))) +1 
            ## store overlaps of target k with query
            ## values to calculate ratios in \code{getOverlapStats(ovl)}:
            ## Jaccard measure: J(query,target) = intersect/union
            ## ratio:           ratio = qlen/tlen
            ## target coverage: tgcov=intersect/tlen
            ## query coverage:  qucov=intersect/qlen
            result <- c(target=k, query=id, 
                        intersect=intersect, union=union, 
                        qlen=qlen, tlen=tlen)
            ## relative query position wrt target
            if ( details ) {
                qpos <- NULL #factor(levels=c("left","covers","inside","right"))
                ## TODO: number positions from left to right
                if ( id %in% which(left) )   qpos <- c(qpos,"left")
                if ( id %in% which(covers) ) qpos <- c(qpos,"covers")
                if ( id %in% which(inside) ) qpos <- c(qpos,"inside")
                if ( id %in% which(right) )  qpos <- c(qpos,"right")

                ## debugging: these should NOT happen
                if ( is.null(qpos) )
                    stop(k, "ERROR: must have one qpos")
                if ( length(qpos)>1 )
                    stop(k, "ERROR: can only have one qpos")
                tps <- c(tps,qpos)
            }
            ovl <- rbind(ovl, result)
        }
        ## RANK: 1 for best match!
        ## rank by Jaccard=intersect/union
        jc <- ovl[,"intersect"]/ovl[,"union"]
        qrank <- rank( -unlist(jc), ties.method="min") 
        ## attempt to untangle ties by second criterium: query coverage
        ##if ( sum(duplicates(qrank))>0 ) {
        if ( sum(qrank==1)>1 ) {
            cat(paste("WARNING:",sum(qrank==1),
                      "with qrank 1 (highest target coverage) for target",
                      tids[k],":", paste(qids[idx[qrank==1]],collapse=";")))
            if ( untie ) {
                cat(paste(" - un-tieing by order"))
                ## NOTE: assignment of first hit to rank 1
                qrank[order(qrank)] <- 1:length(qrank)
                cat(paste(":",paste(qids[idx[qrank==1]],collapse=";")))
            }
            cat(paste(".\n"))
        }
        
        rownames(ovl) <- NULL
        if ( details )  ovl <- data.frame(ovl,qrank=qrank,qpos=tps)
        else ovl <- data.frame(ovl,qrank=qrank)

        ## sort by qrank, top-hit (qrank==1) first
        ## TODO: test this
        if ( sort ) 
            ovl <- ovl[order(ovl[,"qrank"]),]
        
        ## collapse all query hits into ;-separated lists
        if ( collapse ) {
            covl <- rep(NA, ncol(ovl))
            names(covl) <- colnames(ovl)
            for ( j in 1:ncol(ovl) )
                covl[j] <- paste(ovl[,j],collapse=";")
            covl["target"] <- k
            covl["tlen"] <- tlen
            ovl <- covl
        }
        
        overlaps <- rbind(overlaps, ovl)
    }
    rownames(overlaps) <- NULL

    ## TODO - add target rank
    
    data.frame(overlaps, stringsAsFactors=FALSE)
}

#' Statistics of overlaps between two segment sets. 
#' @param ovl overlap table from \code{\link{segmentOverlap}}
#' @param ovlth threshold fraction of overlap of target and query to be
#' counted as 'good' hit
#' @param hrng lower and upper thresholds of the ratio CDF (query length/target
#' length); the fraction of 'best' hits within this range will be reported
#' as 'height'
#' @param tnum number of targets in the original call, required for total CDF
#' ('acdf')
#' @param qnum number queries in the original call, just passed to results
#' @param tid ID for the target set, just passed on to results and used in plot
#' @param qid ID for the query set, just passed on to results and used in plot
#' @export
getOverlapStats <- function(ovl, ovlth=.8, hrng=c(.8,1.2), tnum=NA, qnum=NA, qid=NA, tid=NA) {

    ## returning NULL if no query had been found
    if ( !any(!is.na(ovl[,"query"])) )
      return(NULL)
    
    ## input ovl is from \code{segmentOverlap(ovl)}
    ## Calculate ratios:
    ## Jaccard measure: J(query,target) = intersect/union
    ## ratio:           ratio = qlen/tlen
    ## target coverage: tgcov=intersect/tlen
    ## query coverage:  qucov=intersect/qlen
    jaccard <- ovl[,"intersect"]/ovl[,"union"] 
    ratio <- ovl[,"qlen"]/ovl[,"tlen"]
    tgcov <- ovl[,"intersect"]/ovl[,"tlen"] # handle NA correctly?
    qucov <- ovl[,"intersect"]/ovl[,"qlen"] # -"-


    ## ALL HITS - total target coverage
    ## CDF of all hits, i.e. summing target coverages from multiple segments
    ## ?? todo: get total length and calculate ratio
    ## TODO: required?
    alltrg <- NULL
    if ( !is.na(tnum) ) {
        ## calculate total coverage of each target
        alltrg <- unlist(lapply(1:tnum,function(x)
            sum(tgcov[ovl[,"target"]==x])))
        alltrg[is.na(alltrg)] <- 0
        acdf <- ecdf(alltrg)
    }

    ## BEST HITS: rank==1 - hit with maximal coverage of target 
    ## TODO: rm duplicates properly, i.e. also multiple qrank==1!
    rnk <- ovl[,"qrank"]

    ## BEST HIT - total Jaccard
    J <- sum(ovl[rnk==1,"intersect"])/sum(ovl[rnk==1,"union"])

    ## total Jaccard
    ## todo: total Jaccard requires to account for unmachted qlen/tlen
    ##Jtot <- sum(ovl[,"intersect"])/sum(ovl[,"union"])
    
    ## 'GOOD' HIT: minimal mutual coverages > threshold
    ##             only used to calculate sum below
    sgt <- tgcov[rnk==1] >= ovlth & qucov[rnk==1]>= ovlth


    ## best hit target coverage distributions
    ## TODO: required?
    tgcovdist <- hist(tgcov[rnk==1],breaks=0:20/20,plot=FALSE)
    
    ## cumulative dist function of target coverage

    ## CDF of best hit target coverage
    ## TODO: required?
    tcdf <- ecdf(tgcov[rnk==1]) 

    ## relative best hit segment/test length ratio
    ## CDF of relative best hit ratio, i.e. w/o non-covered targets
    rat <- ratio[rnk==1]
    rrcdf <- ecdf(rat)   ## NA not counted!


   
    ## absolute best hit segment/test length ratio
    ## CDF of absolute best hit ratio, i.e., including non-covered targets
    ## TODO: align this with tnum
    if ( !is.na(tnum) ) {
        rat[is.na(rat)] <-0  # NOTE: difference only for add.na=TRUE
        if ( length(rat)<tnum ) # todo: why?
            rat <- c(rat,rep(0,tnum-length(rat) ))
        rcdf <- ecdf(rat)   
    }

    ## HEIGHTS: target recovery fraction within threshold - get from CDF
    height <- c(rcdf(hrng[1]),rcdf(hrng[2]))
    
    ## relative jaccard (intersect/union)
    ## CDF of relative jaccard ratio, i.e. w/o non-covered targets
    jac <- jaccard[rnk==1]
    rjcdf <- ecdf(jac)   ## NA not counted!

    ## absolute jaccard segment/test length ratio
    ## CDF of absolute jaccard ratio, i.e., including non-covered targets
    ## TODO: align this with tnum
    if ( !is.na(tnum) ) {
        jac[is.na(jac)] <-0  # NOTE: difference only for add.na=TRUE
        if ( length(jac)<tnum ) # todo: why?
            jac <- c(jac,rep(0,tnum-length(jac) ))
        jcdf <- ecdf(jac)   
    }
   
    ## SINGLE OPTIMIZATION MEASURES
    ## distrubtion of queries per target
    tab <- table(ovl[,"target"])
    brks <- max(100, max(tab))
    tgnumdist <- hist(table(ovl[,"target"]), breaks=brks, plot=FALSE)
    ## TODO is mean a good measure?
    numhit <- mean(tab)            ## MINIMIZE: ~ too much fragmentation 
    hitnum <- sum(sgt, na.rm=TRUE) ## MAXIMIZE: ~ # of 'good' hits!

    ## LIST OF CDFs
    CDF <- list(rcdf=rcdf,rrcdf=rrcdf,tcdf=tcdf,acdf=acdf,
                  jcdf=jcdf,rjcdf=rjcdf,
                  qid=qid, tid=tid) # required when results are resorted
    ## LIST OF DISTRIBUTIONS
    DIST <- list(tgcovdist=tgcovdist, # best hit target coverage distribution
                 tgnumdist=tgnumdist) # num. of hits per target distribution
                
    ## TODO: other values? 
    ## TODO: which of these is actually useful except rcdf and rrcdf?
    res <- list(CDF=CDF, # list of CDFs of different measures
                DIST=DIST, # list of distributions
                height=height, # fraction of best hit ratios within hrng
                jaccard=J,     # MAX: jaccard measure of best hits
                hitnum=hitnum, # MAX: num. of hits with minimal mutual coverage
                numhit=numhit, # MIN: avg. num. of hits per target
                qnum=qnum, tnum=tnum,# number of queries/targets 
                qid=qid, tid=tid)
    class(res) <- "overlapStats"
    return(res)
}



### SEGMENT READ STATISTICS

#' Calculates phase distributions in segments, i.e., the circular mean and
#' standard deviation and the R statistics, copied from from package
#' \code{circstats} and after
#' http://en.wikipedia.org/wiki/Directional_statistics#Measures_of_location_and_spread and
#' https://en.wikipedia.org/wiki/Mean_of_circular_quantities
#' TODO: recognize bimodal?
#' @param phs vector of phases in degrees
#' @param w optional weights for weighted mean
#' @export
phaseDist <- function(phs, w) {
  avg <- c(mean=NA, r=NA, na=sum(is.na(phs))/length(phs))
  if ( !missing(w) ) w <- w[!is.na(phs)]
  phs <- phs[!is.na(phs)] * pi/180
    
  if ( length(phs)>0 ) {
    if ( !missing(w) ) {
      c <- sum(cos(phs)*w)
      s <- sum(sin(phs)*w)
      n <- sum(w)
    }else {
      c <- sum(cos(phs))
      s <- sum(sin(phs))
      n <- length(phs)
    }
    r <- sqrt(c^2 + s^2)/n
    mean <- atan2(s, c) *180/pi
    mean <- ifelse(mean<=  0,mean + 360, mean)
    avg[1:2] <- c(mean=mean, r=r)
  }
  if ( !missing(w) ) names(avg) <- paste("w",names(avg),sep=".")
  avg
}
#' calculates the mean of p-values and number of significant
#' p-values (<\code{threshold}) and number of NA values
#' @param pvs vector of p-values
#' @param threshold p-value threshold to be counted as signficant
#' @export
pvalDist <- function(pvs,threshold=0.01) {
  pvs[pvs==1] <- NA # convert back to NA
  avg <- c(p.mean=mean(pvs,na.rm=TRUE),
           p.signif=sum(pvs<threshold,na.rm=TRUE)/length(pvs),
           p.na=sum(is.na(pvs))/length(pvs))
  avg
}

#' distribution of read-counts, i.e., the mean, var, min&max and the number
#' of NA values
#' @param rds vector of read-counts, e.g., within a genomic interval (segment)
#' @export
readDist <- function(rds) {
  avg <- c(r.mean=mean(rds,na.rm=TRUE),
           r.var=var(rds,na.rm=TRUE),
           r.min=min(rds,na.rm=TRUE),
           r.max=max(rds,na.rm=TRUE),
           r.na=sum(is.na(rds)))
  avg
}

#' Average read-count  time-series of a genomic interval (segment). Note
#' that the diverse attempts to smooth or filter read-counts before
#' taking the averages did not give good results, and we remained using
#' avg="mean" and no other transformations.
#' @param rds matrix of read-count time-series
#' @param avg function to use as average statistics, \code{mean} or
#' \code{median}
#' @param mean.ratio calculate the mean ratio over time of each read-count
#' time-series before taking the total average
#' @param rm.extreme remove extreme value reads,
#' i.e., where the temporal mean is \code{<0.1} or \code{>0.9} of all means
#' @param endcut fraction of the ends the genomic interval to remove before
#' taking the total average
#' @param k parameter \code{k} for the running median
#' function \code{\link[stats]{runmed}}; the running median will be calculated
#' only if \code{k>1}
#' @param endrule parameter  \code{endrule} for the running median
#' function \code{\link[stats]{runmed}}
#' @export
segmentAverage <- function(rds, avg="mean", mean.ratio=FALSE,
                           rm.extreme=FALSE,
                           endcut=0, k=1, endrule="median") {

    ## not used - doesnt work well
    if ( mean.ratio ) {
        ## store mean
        mn <- mean(rds)
        ## mean-ratio: scale each read time-series by its mean
        rds <- rds/apply(rds,1,mean)
    }

    ## not used: cut ends
    if ( endcut>0 ) {
        N <- nrow(rds)
        end <- round(N*endcut)
        idx <- (1+end):(N-end)
        rds <- rds[idx,]
    }
    ## not used: smooth each position
    if ( k>1 )
      rds <- t(apply(rds,1,runmed,k=k,endrule=endrule))

    ## not used: remove extrema
    if ( rm.extreme ) {
        rmn <- apply(rds,1,mean)
        qu <- quantile(rmn,c(.1,.9))
        rm <- rmn<qu[1] | rmn>qu[2]
        rds[rm,] <- NA
    }

    ## not used: loess - not really working!
    if ( avg=="loess" ) { 
      y <- c(rds)
      X <- rep(1:ncol(rds),each=nrow(rds))
      data <- data.frame(x=X,y=y)
      fit <- loess(y ~ X, data=data, degree = 2) # slow!
      #plot(y ~ x, data=data,pch=19,cex=0.1)
      #j <- order(data$x)
      #lines(data$x[j],fit$fitted[j],col="red",lwd=3)
      return(predict(fit,1:ncol(rds)))
    }
    ##matplot(t(rds),type="l",lty=1,lwd=.5,col="black",log="")
    ##lines(apply(rds,2,mean,na.rm=T),col=2,lwd=3)

    ## not used (see above): ... and scale by real mean
    if ( mean.ratio ) avg <- avg*mn
    
    ## return average (mean or median)
    ## median not used
    apply(rds,2,get(avg,mode="function"), na.rm=TRUE)
}



### TODO - NOT WORKING CODE

## more general approach (faster/better?) of segmentOverlap
## loops over segments instead!
## Here, we join query and target segments and get the order
## of all start and end sides. Then, we go through queries. If a query
## has the same order of start and end sites, then it does NOT overlap
## with any target. [problem: order of equal positions?]
## TODO - NOT WORKING: atm this does NOT find all overlaps!
## requires additional sort of query/target start vs. target/query ends
## to find overlapping - not sure whether this would be any time improvement
## see e.g BIT algo http://europepmc.org/articles/PMC3530906
segmentOverlap.v2 <- function(query, target, details=FALSE, add.na=FALSE) {
    cols <- c("start","end")
    all<-rbind(query[,cols],target[cols]);
    n <- nrow(query)
    m <- nrow(all)
    stord <- order(all[,"start"])
    enord <- order(all[,"end"])
    ## search and classify all insert between queries 1:n
    overlaps <- data.frame()
    ## go through queries and get overlapping targets
    all <- NULL
    for ( j in 1:n ) {
        i <- which(stord == j)
        k <- which(enord == j)
        if ( i<k ) {
            ovl <- stord[(i+1):k] - n
            #ovl2 <- enord[(i-1):k] - n
            ## TODO: instead search end of query in target starts
            ## and start of query in target ends
        
            ## details
            ## scan each target
            qu <- query[j,"start"]:query[j,"end"]
            res <- NULL
            for ( t in ovl ) {
                tg <- target[t,"start"]:target[t,"end"]
                is <- intersect(tg,qu) # intersection
                len <- length(is)   
                totlen <- diff(range(c(qu,tg))) +1 # union
                ## store overlaps of target k with query
                result<-c(query=j, target=t, intersect=len,
                          ratio=round(length(qu)/length(tg),3),#query/target
                          jaccard=round(len/totlen,3), #intersect/union (Jaccard)
                          tgcov=round(len/length(tg),3), #intersect/target
                          qucov=round(len/length(qu),3)) #intersect/query
                ## is isect left/right/inside/over target?
                if ( details ) {}
                res <- rbind(res,result)
            }
            if ( !is.null(res) ) 
                all <- rbind(all,res)
            ## TODO: this doesn't cover all!! diff. stord and enord
            ovl <- cbind(query=rep(j,length(ovl)),target=ovl)
            overlaps <- rbind(overlaps,ovl)
        }        
    }
    #overlaps[order(overlaps[,"target"]),]
    rownames(all) <- NULL
    all[order(all[,"target"]),]
}
