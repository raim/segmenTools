#!/usr/bin/env Rscript

## TESTING SEGMENT LENGTH DIST AND AGAINST ANNOTATED GENES AND TRANSCRIPTS

library(segmenTools)
##segtools <- "~/programs/segmenTools/"
##source(file.path(segtools,"R/segmenTools.R")) # for segment analysis
##source(file.path(segtools,"R/coor2index.R")) # coor2index

library(cluster) # for pam clustering
library(optparse) # command-line options

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")

### OPTIONS
option_list <- list(
    make_option(c("-i", "--infile"), type="character", default="", 
                help="chromosome coordinates of primary segments as produced by clusterSegments.R but without header ('allsegs.csv')"),    
    make_option(c("--chrfile"), type="character", default="",
                help="chromosome index file, providing a sorted list of chromosomes and their lengths in column 3 [default %default]"),
    make_option(c("--qtypes"), type="character", default="", 
                help="sub-set testset in column 'type'"),
    make_option(c("--qtypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
    ##chrfile = $YEASTDAT/chromosomes/sequenceIndex_R64-1-1_20110208.csv
    ## SEGMENT TEST SETTINGS
    make_option(c("--fuse.segs"), action="store_true", default=FALSE,
                help="use FUSE tag from clusterSegments to fuse adjacent segments"),
    make_option(c("--target"), type="character", default="", 
              help="target set of chromosomal segments"),    
    make_option(c("--ttypes"), type="character", default="", 
                help="sub-set testset in column 'type'"),
    make_option(c("--ttypcol"), type="character", default="type", 
              help="name of column with sub-set annotation"),
    make_option(c("--tcolcol"), type="character", default="", 
              help="name of column with sub-set colors for plots"),
    make_option("--ovlth", default=0.8,
                help="overlap threshold (mutual coverage) to be counted as a direct hit; at least ovlth*length must be reached for both query (segments) and targets (test set) [default %default]"),
    ## OUTPUT OPTIONS
    make_option(c("--out.path"), type="character", default=".", 
                help="directory path for output data (figures, csv files)"),
    make_option(c("--testid"), type="character", default="testset", 
                help="ID of the testset, used for folder and file names"),
    make_option(c("-v", "--verb"), type="integer", default=1, 
                help="0: silent, 1: main messages, 2: warning messages"),
    make_option(c("--fig.type"), type="character", default="png",
                help="figure type, png or pdf [default %default]"),
    make_option(c("--save"), action="store_true", default=FALSE,
                help="save overlap data as RData file (big!)"))

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## process comma-separated list arguments
lst.args <- c(ttypes="character",qtypes="character")
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


### START

## load chromosome index
if ( verb>0 )
    cat(paste("Loading chromosome index file:", chrfile, "\n"))
cf <- read.table(chrfile,sep="\t",header=FALSE)
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos

## READ SEGMENTS TO BE TESTED 
if ( verb>0 )
    cat(paste("LOADING SEGMENTS:", infile, "\n"))
segs <- read.table(infile,sep="\t",header=TRUE)



#' uses a logical column to fuse adjacent segments, i.e. the
#' lower segment i with segment i-1, where segments
#' is an ordered list of segments
fuseSegments <- function(segments, col="fuse", val=1) {
    fuse <- segments[2:nrow(segments),col] == val
    cat(paste("NOTE: FUSING", sum(fuse), "SEGMENTS, from segment types:\n",
              paste(unique(segments[fuse,"type"]),collapse="; "),"\n"))
    fsegs <- segments[c(TRUE,!fuse),]
    
    fsegs[,"end"] <- segments[c(!fuse,TRUE),"end"]
    fsegs
}
## use "fuse" column from clusterSegments.R fuse filter!
if ( fuse.segs ) {
    if ( verb>0 )
        cat(paste("fusing segments by fuse tag\n"))
    segs <- fuseSegments(segs,col="fuse",val=1)
}

## replace genome coordinates by continuous index
segs <- coor2index(segs,chrS)

## filter types
qtypes <- qtypes[qtypes!=""]
if ( length(qtypes)==0 ) {
    qtypes <- as.character(segs[,qtypcol])
    qtypes[qtypes==""] <- "na"
    qtypes[is.na(qtypes)] <- "na"
} else {
    ## filter to types given by cmdline arg qtypes
    segs <- segs[as.character(segs[,tpcol])%in%qtypes,]
    qtypes <- as.character(segs[,tpcol]) # get remaining
}

## split by type
lst <- split(segs,segs[,qtypcol])
sgtypes <- names(lst)



debug <- FALSE
if ( debug ) { ## test only large number
    sgnum <- unlist(lapply(lst,nrow))
    sgtypes <- sgtypes[sgnum>1000] ### TODO: remove this and fuse types
}

## define colors and pch for segment types
sgcols <- rainbow(length(sgtypes))
sgpchs <- rep(1:17, len=length(sgtypes)) #only 17 of 18 to avoid multiple of lty
sgltys <- rep(1:6, len=length(sgtypes))
names(sgcols) <- names(sgpchs) <- names(sgltys) <- sgtypes

## get segment classes 
sgclasses <- getSegmentClasses(sgtypes, sep="_")

## col, pch and lty for classes
sgclcols <- rainbow(length(sgclasses))
sgclpchs <- rep(1:17, len=length(sgclasses)) #only 17 to avoid lty multiple
sgclltys <- rep(1:6, len=length(sgclasses))
names(sgclcols) <- names(sgclpchs) <-names(sgclltys) <- sgclasses

## new: segment class table
sgcltab <- getSegmentClassTable(sgtypes, sep="_")

### START ANALYSIS

## output dirs - generate locally in cwd
dir.create(out.path) # length distributions
## output dir
dir.create(file.path(out.path,testid),recursive=TRUE) # test set


### TEST AGAINST OTHER DATA SETS

if ( verb>0 )
    cat(paste("LOADING TEST SETS:", testid, "-", target, "\n"))

trgs <- read.table(target,sep="\t",header=TRUE,
                   stringsAsFactors=FALSE, comment.char = "")
trgs <- coor2index(trgs,chrS) # map to continuous index!

tpcol <- ttypcol
cpcol <- tcolcol

## test types
ttypes <- ttypes[ttypes!=""]
if ( length(ttypes)==0 ) {
    ttypes <- as.character(trgs[,tpcol])
    ttypes[ttypes==""] <- "na"
    ttypes[is.na(ttypes)] <- "na"
} else {
    ## filter to types given by cmdline arg ttypes
    trgs <- trgs[as.character(trgs[,tpcol])%in%ttypes,]
    ttypes <- as.character(trgs[,tpcol]) # get remaining
}
test.types <- sort(names(table(ttypes))) # sort by number

if ( cpcol=="" ) {
    # TODO: rgb range1:length(test.types)
    tcols <- sub("00$","",rainbow(length(test.types),alpha=0))
    names(tcols) <- test.types
} else {
    ## just take first color for each type,
    ## user must ensure that these are unique
    tcols <- sapply(test.types,
                    function(x)
                        as.character(trgs[which(trgs[,tpcol]==x)[1],cpcol]))
}
tcols["na"] <- "#939393"


## results
ovlstats <- rep(list(NA),length(test.types))
names(ovlstats) <- test.types

## loop over targets (transcript and ORF data sets from SGD)
for ( test.type in test.types ) {

    if ( verb> 0 )
        cat(paste(test.type, "\t",time(),"\n"))
    
    target <- trgs[ttypes==test.type,]

    ## result list: overlap statistics
    ostat <- rep(list(NA),length(sgtypes))
    names(ostat) <- sgtypes
    ovlstats[[test.type]] <- ostat

    ## loop over queries (segment types)
    for ( type in sgtypes ) {
        if ( verb>0 )
            cat(paste("#",which(sgtypes==type), type,
                      "overlap with",test.type,"\t",time(),"\n"))
        sgs <- lst[[type]]
        
        ## for testing
        devel <- FALSE
        if ( devel )
            sgs <- sgs[1:1000,]

        
        ovl <- segmentOverlap(query=sgs, target=target,
                              add.na=TRUE, details=TRUE,
                              untie=FALSE, collapse=FALSE, sort=FALSE,
                              msgfile=file("stdout"))
        if ( !any(!is.na(ovl[,"query"])) ) next
        sts <- getOverlapStats(ovl,ovlth=ovlth, hrng=c(.8,1.2),
                               tnum=nrow(target),qnum=nrow(sgs),
                               qid=type, tid=test.type)

        ovlstats[[test.type]][[type]] <- sts

        ## TODO: save results as files 
    }
}

if ( save ) 
  save.image(file.path(out.path,testid,"overlaps.RData"))

## COLLECT DATA AND WRITE OUT RESULTS
## export hitnum, jaccard, numhit and ratio thresholds 
for ( test.type in test.types ) {

    #if ( is.null(ovlstats[[test.type]]) ) next
    
    covlStats <- collectOvlStats(ovlstats, type=test.type)

    ids <- covlStats$nms
    jaccard <- covlStats$jaccard # jaccard of best hits
    height <- covlStats$height # target recovery fraction within threshold
    hitnum <- covlStats$hitnum # total number of 'good' hits 
    numhit <- covlStats$numhit # average number of hits per target
    tnum <- covlStats$tnum # number of tested targets
    tnum <- covlStats$tnum # number of tested targets

    result <- data.frame(ID=ids, tnum=tnum, hits=hitnum, Jaccard=jaccard,
                         hits.per.target=numhit,
                         ratio.low=height[,1], ratio.high=height[,2])
    
    ## write out table of segmentation characteristics
    file.name <- file.path(out.path,testid,paste("segmentRecovery_",test.type,
                                          ".csv",sep=""))
    write.table(result,file=file.name, sep="\t",
                col.names=TRUE,row.names=FALSE,quote=FALSE)
}

#### PLOT OF OVERLAP STATISTICS

## remove empty results!
rm <- unlist(lapply(ovlstats, function (x)
                    !any(unlist(lapply(x, function(y) y$hitnum))>0)))
ovlstats <- ovlstats[!rm]
test.types <- test.types[!rm]

### PLOT BY TEST TYPES
for ( test.type in test.types ) {

    covlStats <- collectOvlStats(ovlstats, type=test.type)

    CDF <- covlStats$CDF
    jaccard <- covlStats$jaccard
    height <- covlStats$height
    hitnum <- covlStats$hitnum
    numhit <- covlStats$numhit
    qnum <- covlStats$qnum
    tnum <- covlStats$tnum
    nms <- covlStats$nms

    ## TODO: mv this to getOverlapStats
    
    ## MAX vs. MIN CLUSTERING  cluster jaccard vs. numhits vs segment classes
    ## TODO: also include height in clustering?
    K <- 6
    dat <- cbind((jaccard-min(jaccard))/(max(jaccard)-min(jaccard)),
                 (numhit-min(numhit))/(max(numhit)-min(numhit)))
    pm <- pam(dat, K)
    cllst <- apply(sgcltab,2,unique)
    allcl <- unlist(sapply(1:length(cllst),
                           function(x) paste(names(cllst)[x],
                                             cllst[[x]],sep=".")))
    enum <- matrix(NA, nrow=K, ncol=length(allcl))
    colnames(enum) <- allcl
    rownames(enum) <- 1:K
    pval <- enum
    ## get cumulative hypergeometric distribution of clustering vs. segments
    for ( i in 1:K ) 
      for ( j in 1:ncol(sgcltab) ) {
          clcl <-  clusterCluster(pm$clustering==i,sgcltab[,j],
                                  plot=FALSE,verbose=FALSE)
          cln <- colnames(sgcltab)[j]
          cln <- paste(cln,colnames(clcl$overlap),sep=".")
          if ( "TRUE" %in% rownames(clcl$overlap) ) {
              enum[i,cln] <- clcl$overlap["TRUE",]
              pval[i,cln] <- clcl$p.value["TRUE",]
          }
      }

    file.name <- file.path(out.path,testid,
                           paste(test.type,"_segmentationClusters",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(1,.7,.1,.1),mgp=c(1.75,.5,0))
    image_matrix(-log2(pval) ,text=enum, axis=1:2,
                 col=c("#FFFFFF",rev(grey.colors(20))),
                 axis2.col=1:nrow(pval),
                 xlab=NA,ylab=NA)
    axis(1, at=cumsum(unlist(lapply(cllst, length)))+.5, tck=-1,labels=NA) 
    dev.off()
   
    ## TODO: plot by segment classes as length dist

    ## new lower and upper threshold of ratio
    ## OPT: uppler left, 
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_ratioTotal_lh_clustered",sep=""))
    plotdev(file.name,width=10,height=5,type=fig.type)
    par(mfcol=c(1,2),mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0))
    plot(height,xlab="fraction: ratio < 0.8",ylab="fraction: ratio < 1.2",
         col=NA)
    legend("topleft","good",bty="n")
    points(height,col=pm$clustering,pch=sgpchs[nms])
    abline(v=.2,lty=2)
    abline(h=.8,lty=2)
    par(mai=c(1,.7,.1,.1))
    image_matrix(-log2(pval) ,text=enum, axis=1:2,
                 col=c("#FFFFFF",rev(grey.colors(20))),
                 axis2.col=1:nrow(pval),
                 xlab=NA,ylab=NA)
    axis(1, at=cumsum(unlist(lapply(cllst, length)))+.5, tck=-1,labels=NA)
    dev.off()

    file.name <- file.path(out.path,testid,
                           paste(test.type,"_ratioTotal_lh",sep=""))
    plotdev(file.name,width=10,height=5,type=fig.type)
    par(mfcol=c(1,2),mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0))
    plot(height,xlab="fraction: ratio < 0.8",ylab="fraction: ratio < 1.2",
         col=sgcols[nms],pch=sgpchs[nms])
    abline(v=.2,lty=2)
    abline(h=.8,lty=2)
    par(mai=c(.1,.1,.1,.1))
    plot(0, col=NA,axes=FALSE,ylab=NA,xlab=NA)
    legend("topleft", legend=nms, col=sgcols[nms],pch=sgpchs[nms],
           cex=.5,bty="n",ncol=2)
    dev.off()
    
    ## CDF of absolute best hit CDF (rcdf)
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_ratioTotal",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot_cdfLst(x=seq(0,2,.05), CDF=CDF, type="rcdf", col=sgcols, lty=sgltys,
                h=c(.2,.8), v=c(ovlth,2-ovlth), #c(0.8,1.2),
                xlab="ratio: query length/target length")
    legend("topleft",paste(test.type,"-",tnum))
    dev.off()

    ## CDF of absolute best hit CDF (rcdf)  - cluster colors
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_ratioTotal_clustered",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot_cdfLst(x=seq(0,2,.05), CDF=CDF, type="rcdf",
                col=pm$clustering, lty=sgltys,
                h=c(.2,.8), v=c(ovlth,2-ovlth), #c(0.8,1.2),
                xlab="ratio: query length/target length")
    legend("topleft",paste(test.type,"-",tnum))
    dev.off()
      
    ## MAX vs. MIN: jaccard of best hits vs. num hits per target 
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_jaccard_fragmentation",sep=""))
    plotdev(file.name,width=10,height=4,type=fig.type)
    par(mfcol=c(1,2),mai=c(1,1,.1,.1))
    plot(jaccard,numhit,
         col=sgcols[nms],pch=sgpchs[nms],
         ylab="average hits per target sequence",
         xlab="jaccard: intersect/union")#,
    par(mai=c(.1,.1,.1,.1))
    plot(0, col=NA,axes=FALSE,ylab=NA,xlab=NA)
    legend("topleft", legend=nms, col=sgcols[nms],pch=sgpchs[nms],
           cex=.5,bty="n",ncol=2)
    dev.off()

    ## MAX vs. MIN CLUSTERS - JACCARD
    file.name <- file.path(out.path,testid,
                     paste(test.type,"_jaccard_fragmentation_clustered",sep=""))
    plotdev(file.name,width=10,height=4,type=fig.type)
    par(mfcol=c(1,2),mai=c(1,1,.1,.1))
    plot(jaccard,numhit,col=NA,         
         ylab="average hits per target sequence",
         xlab="jaccard: intersect/union")#,
    legend("bottomright","good",bty="n")
    points(jaccard,numhit,col=pm$clustering,pch=sgpchs[nms])
    par(mai=c(1,.7,.1,.1))
    image_matrix(-log2(pval) ,text=enum, axis=1:2,
                 col=c("#FFFFFF",rev(grey.colors(20))),
                 axis2.col=1:nrow(pval),
                 xlab=NA,ylab=NA)
    axis(1, at=cumsum(unlist(lapply(cllst, length)))+.5, tck=-1,labels=NA) 
    dev.off()

    ## MAX vs. MIN: ratio heights vs. num hits per target
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_ratio_fragmentation",sep=""))
    plotdev(file.name,width=10,height=4,type=fig.type)
    par(mfcol=c(1,2),mai=c(1,1,.1,.1))
    plot(apply(height,1,diff), numhit,
         col=sgcols[nms],pch=sgpchs[nms],
         ylab="average hits per target sequence",
         xlab="fraction: 0.8 < ratio < 1.2")#,
    par(mai=c(.1,.1,.1,.1))
    plot(0, col=NA,axes=FALSE,ylab=NA,xlab=NA)
    legend("topleft", legend=nms, col=sgcols[nms],pch=sgpchs[nms],
           cex=.5,bty="n",ncol=2)
    dev.off()
    
    ## MAX vs. MIN CLUSTERS - good hits per target vs. num hits per target 
    file.name <-file.path(out.path,testid,
                          paste(test.type,"_ratio_fragmentation_clustered",
                                sep=""))
    plotdev(file.name,width=10,height=4,type=fig.type)
    par(mfcol=c(1,2),mai=c(1,1,.1,.1))
    plot(apply(height,1,diff), numhit,
         col=pm$clustering,pch=sgpchs[nms],
         ylab="average hits per target sequence",
         xlab="fraction: 0.8 < ratio < 1.2")#,
    imgdat <- t(apply(-log2(pval), 2, rev))
    par(mai=c(1,.7,.1,.1))
    image_matrix(-log2(pval) ,text=enum, axis=1:2,
                 col=c("#FFFFFF",rev(grey.colors(20))),
                 axis2.col=1:nrow(pval),
                 xlab=NA,ylab=NA)
    axis(1, at=cumsum(unlist(lapply(cllst, length)))+.5, tck=-1,labels=NA)
    dev.off()

### BELOW PERHAPS NOT REQUIRED

    ## CDF PLOTS

    ## CDF of jaccard (jcdf)
    file.name <- file.path(out.path,testid,paste(test.type,"_jaccard_cdf_clustered",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot_cdfLst(x=seq(0,1.1,.05), CDF=CDF, type="rjcdf",
                col=pm$clustering, lty=sgltys,
                h=c(.2,.8), v=c(ovlth,2-ovlth), #c(0.8,1.2),
                xlab="cumulative jaccard: intersect/union")
    legend("topleft",paste(test.type,"-",tnum))
    dev.off()

    ## CDF of relative best hit CDF (rrcdf)
    file.name <- file.path(out.path,testid, paste(test.type,"_ratio",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot_cdfLst(x=seq(0,2,.05), CDF=CDF, type="rrcdf", col=sgcols, lty=sgltys,
                h=c(.2,.8), v=c(ovlth,2-ovlth), #c(0.8,1.2),
                xlab="relative ratio: query length/target length")
    legend("topleft",paste(test.type,"-",tnum))
    dev.off()


    ## CDF of best hit target coverage
    file.name <-file.path(out.path,testid,
                          paste(test.type,"_totalCoverage",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot_cdfLst(x=seq(0,1.1,.05), CDF=CDF, type="tcdf", col=sgcols, lty=sgltys,
                h=c(.2,.8), v=c(ovlth,2-ovlth), #c(0.8,1.2),
                xlab="coverage of test set")
    legend("topleft",paste(test.type,"-",tnum))
    dev.off()
   
     ## CDF of best hit target coverage - cluster colors
    file.name <-file.path(out.path,testid,
                          paste(test.type,"_totalCoverage_clustered",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.1,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot_cdfLst(x=seq(0,1.1,.05), CDF=CDF, type="tcdf",
                col=pm$clustering, lty=sgltys,
                h=c(.2,.8), v=c(ovlth,2-ovlth), #c(0.8,1.2),
                xlab="coverage of test set")
    legend("topleft",paste(test.type,"-",tnum))
    dev.off()
    
    ## SUMMARY OF CDF of absolute best hit ratio
    ## fraction of mutual coverage between 0.8 and 1.2
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_ratioTotal_rng_clustered",sep=""))
    plotdev(file.name,width=2+.2*length(CDF),height=6,type=fig.type)
    par(mai=c(3.3,.75,.1,.1))
    plot(0,col=NA,ylim=c(0,1.2),xlim=c(0,length(CDF)+1),
         axes=FALSE,xlab=NA,ylab=NA)
    axis(2)
    abline(h=c(.2,.8),lty=2)
    leg <- NULL
    for ( i in 1:length(CDF) ) {
        if ( "rcdf" %in% names(CDF[[i]]) ) {
            nm <- CDF[[i]]$qid
            if ( diff(height[i,]) > 0 )
                arrows(x0=i,x1=i,y0=height[i,1],y1=height[i,2],
                       col=pm$clustering[i],lty=sgltys[nm],
                       length=0.1, angle=90, code=3)
            text(i, 1.1, round(diff(height[i,]),3),srt=90)
        }
    }
    axis(1,at=1:length(CDF),labels=nms,las=2)
    dev.off()
   
    ## MAX vs. MIN: good hits per target vs. num hits per target
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_coverage_fragmentation",sep=""))
    plotdev(file.name,width=10,height=4,type=fig.type)
    par(mfcol=c(1,2),mai=c(1,1,.1,.1))
    plot(hitnum/tnum, numhit,
         col=sgcols[nms],pch=sgpchs[nms],
         ylab="average hits per target sequence",
         xlab="total number of 'good' hits")#,
    par(mai=c(.1,.1,.1,.1))
    plot(0, col=NA,axes=FALSE,ylab=NA,xlab=NA)
    legend("topleft", legend=nms, col=sgcols[nms],pch=sgpchs[nms],
           cex=.5,bty="n",ncol=2)
    dev.off()

    ## MAXIMIZE COVERAGE: good hits per target
    cvg <- hitnum/tnum
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_coverage",sep=""))
    plotdev(file.name,,width=2+.2*length(CDF),height=6,type=fig.type)
    par(mai=c(3.3,.75,.1,.1),mgp=c(1.75,.5,0))
    plot(1:length(cvg),cvg,type="b",col=sgcols[nms],pch=sgpchs[nms],axes=FALSE,
         xlab=NA,ylab=paste("number of 'good' hits",ovlth))
    axis(2)
    axis(1,at=1:length(cvg),labels=nms,las=2)
    dev.off()
    
    ## MAXIMIZE JACCARD: good hits per target
    file.name <- file.path(out.path,testid,
                           paste(test.type,"_jaccard",sep=""))
    plotdev(file.name,,width=2+.2*length(CDF),height=6,type=fig.type)
    par(mai=c(3.3,.75,.1,.1),mgp=c(1.75,.5,0))
    plot(1:length(jaccard),jaccard,type="b",col=sgcols[nms],pch=sgpchs[nms],
         axes=FALSE,xlab=NA,ylab="jaccard: intersect/union")
    axis(2)
    axis(1,at=1:length(jaccard),labels=nms,las=2)
    dev.off()

    ## MINIMIZE FRAGMENTATION: num hits per target 
    nmh <- numhit
    file.name <-file.path(out.path,testid,
                          paste(test.type,"_fragmentation",sep=""))
    plotdev(file.name,,width=2+.2*length(CDF),height=6,type=fig.type)
    par(mai=c(3.3,.75,.1,.1),mgp=c(1.75,.5,0))
    plot(1:length(nmh),nmh,type="b",col=sgcols[nms],pch=sgpchs[nms],axes=FALSE,
         xlab=NA,ylab=paste("avg segments per target"))
    axis(2)
    axis(1,at=1:length(nmh),labels=nms,las=2)
    dev.off()
}


dir.create(file.path(out.path,testid,"segtypes"))
### PLOT BY SEGMENT TYPES
for ( type in sgtypes ) {

    ## re-order result list
    hitnum <- rep(NA,length(test.types))
    CDF <- rep(list(NA),length(test.types))
    height <- matrix(NA, ncol=2, nrow=length(test.types))
    rownames(height) <- names(hitnum) <- names(CDF) <- test.types
    numhit <- tnum <- hitnum
    for ( test.type in test.types ) {
        CDF[[test.type]] <- ovlstats[[test.type]][[type]]$CDF
        CDF[[test.type]]$name <- test.type
        numhit[test.type] <-  ovlstats[[test.type]][[type]]$numhit
        hitnum[test.type] <-  ovlstats[[test.type]][[type]]$hitnum
        tnum[test.type] <-  ovlstats[[test.type]][[type]]$tnum
        height[test.type,] <- ovlstats[[test.type]][[type]]$height
    }
    
    ## plot
    file.name <- file.path(out.path,testid,"segtypes",
                           paste(type,"_ratioTotal",sep=""))
    plotdev(file.name,width=5,height=5,type=fig.type)
    par(mai=c(.75,.75,.5,.1),mgp=c(1.75,.5,0),xaxs="i")
    leg <- NULL
    xmax <- 3
    x <- seq(0,xmax,.05)
    plot(x, CDF[[1]]$rcdf(x),type="l",main=NA,col=NA,lty=1,
         xlab="ratio: query length/target length",xlim=c(0,xmax),
         ylab="cum.dist.func.",ylim=c(0,1.1))
    abline(v=c(ovlth,2-ovlth),lty=2)
    abline(h=c(.2,.8),lty=2)
    abline(h=0:1, lty=2, col="gray",lwd=.75)
    for ( i in 1:length(CDF) ) {
      col <- tcols[CDF[[i]]$name] # todo: tid
      if ( !is.null(CDF[[i]]$rcdf) ) {
        lines(x,CDF[[i]]$rcdf(x),col=paste(col,"AA",sep=""),lty=1,lwd=2) 
        leg <- c(leg, CDF[[i]]$name) # todo: tid
      }
    }
    mtext(paste(testid,"recovery by segments"),3,0,cex=1.5)
    legend("bottomright",leg,col=tcols[leg],lty=1,lwd=2,bty="n")
    legend("topleft",paste(type))
    dev.off()

    file.name <- file.path(out.path,testid,"segtypes",
                           paste(type,"_ratioTotal_rng",sep=""))
    plotdev(file.name,width=2+.2*length(CDF),height=6,type=fig.type)
    par(mai=c(3.3,.75,.1,.1))
    plot(0,col=NA,ylim=c(0,1.2),xlim=c(0,length(CDF)+1),
         axes=FALSE,xlab=NA,ylab=NA)
    axis(2)
    abline(h=c(.2,.8),lty=2)
    leg <- NULL
    for ( i in 1:length(CDF) ) {
        if ( "rcdf" %in% names(CDF[[i]]) ) {
            if ( diff(height[i,]) > 0 )
                arrows(x0=i,x1=i,y0=height[i,1],y1=height[i,2],
                       col=tcols[names(CDF)[i]],lty=1,length=0.1,
                       angle=90,code=3)
            text(i, 1.1, round(diff(height[i,]),3),srt=90)
        }
    }
    axis(1,at=1:length(CDF),labels=names(CDF),las=2)
    dev.off()
}

quit(save="no")


### TEST-SET ANALYSES
#### TODO; allow segmentRevery to work with transcripts vs features

## ORF recovery by transcripts!
trcl <- segmentOverlap(query=clorf,target=genes,
                       add.na=TRUE,details=TRUE,sort=TRUE)
## plot cluster gene recovery by published transcripts
trst <- getOverlapStats(trcl,tnum=nrow(genes),
                        qnum=nrow(clorf),qid="clorf", ovlth=ovlth)
file.name <- file.path(paste("orfXtranscripts_ratioTotal",sep=""))
plotdev(file.name,width=5,height=5,type=fig.type)
par(mai=c(.75,.75,.5,.1),mgp=c(1.75,.5,0),xaxs="i")
xmax <- 2
x <- seq(0,xmax,.05)
plot(x, trst$CDF$rcdf(x),type="l",col=1,lty=1,xlim=c(0,xmax),xlab="transcript length/ORF length",ylab="cum.dist.func.",main=NA,ylim=c(0,1))
abline(v=c(ovlth,2-ovlth),lty=2)
abline(h=c(.2,.8),lty=2)
abline(h=0:1, lty=2, col="gray",lwd=.75)
mtext(paste("ORF recovery by ORF transcripts"),3,0,cex=1.5)
legend("bottomright","ORF recovery by ORF transcripts",bty="n")
legend("topleft","")
dev.off()

## cluster ORF recovery by transcripts
## CLUSTER ORFS
test <- "clgenes"
trgs <- genes
#trgs[trgs[,"CL_rdx"]=="","CL_rdx"] <- "na"
tpcol <- "CL_rdx"
cpcol <- "CL_rdx_col"
## output dir
dir.create(file.path(out.path,test),recursive=TRUE)

## test types
ttypes <- as.character(trgs[,tpcol])
ttypes[ttypes==""] <- "na"
ttypes[is.na(ttypes)] <- "na"
test.types <- sort(names(table(ttypes)))
tcols <- sapply(test.types,
                function(x)
                    as.character(trgs[which(trgs[,tpcol]==x)[1],cpcol]))

tcols["na"] <- "#939393"
                
## results
ovlstats <- rep(list(NA),length(test.types))
names(ovlstats) <- test.types

query <- transcripts[transcripts[,"type"]=="ORF",]
for ( test.type in test.types ) {

    target <- trgs[ttypes==test.type,,drop=FALSE]
     
    cat(paste("transcript overlap with",test.type,"\t",time(),"\n"))
    
    ovl <- segmentOverlap(query=query,target=target, add.na=TRUE,details=TRUE)
    sts <- list(CDF=NA,hitnum=0,qid=test.type,tid="ORF")

    if ( any(!is.na(ovl[,"query"])))
        sts <- getOverlapStats(ovl,ovlth=.8,qid=test.type,tid="ORF",
                               tnum=nrow(target),qnum=nrow(query))
    
    ovlstats[[test.type]] <- sts
}
rm <- unlist(lapply(ovlstats, function (x)
                    !any(x$hitnum>0)))
ovlstats <- ovlstats[!rm]
test.types <- test.types[!rm]

## fuse results to primary result list - TODO: not fused?
CDF <- rep(list(NA),length(test.types))
names(CDF) <- test.types
for ( test.type in test.types ) {
    CDF[[test.type]] <- ovlstats[[test.type]]$CDF
    CDF[[test.type]]$name <- test.type
}
 
## plot
file.name <- file.path(paste("clusterXtranscripts_ratioTotal",sep=""))
plotdev(file.name,width=5,height=5,type=fig.type)
par(mai=c(.75,.75,.5,.1),mgp=c(1.75,.5,0),xaxs="i")
leg <- NULL
xmax <- 3
x <- seq(0,xmax,.05)
plot(x, CDF[[1]]$rcdf(x),type="l",col=NA,lty=1,xlim=c(0,xmax),xlab="ratio: query length/target length",ylab="cum.dist.func.",main=NA,ylim=c(0,1))
abline(v=c(ovlth,2-ovlth),lty=2)
abline(h=c(.2,.8),lty=2)
abline(h=0:1, lty=2, col="gray",lwd=.75)
for ( i in 1:length(CDF) ) {
    col <- tcols[CDF[[i]]$name]
    if ( "rcdf" %in% names(CDF[[i]]) )
        lines(x,CDF[[i]]$rcdf(x),col=paste(col,"AA",sep=""),lty=1,lwd=2) 
    leg <- c(leg, CDF[[i]]$name)
}
mtext(paste(test,"recovery by ORF transcripts"),3,0,cex=1.5)
legend("topleft","",bty="n")
legend("bottomright",leg,col=tcols[leg],lty=1,lwd=2,bty="n")
dev.off()


quit(save="no")

