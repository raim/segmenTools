#!/usr/bin/env Rscript

## TESTING SEGMENT LENGTH DIST (from testSegments.R 20161221)

## genomeBrowser utils

library(segmenTools)
#segtools <- "~/programs/segmenTools/"
#source(file.path(segtools,"R/segmenTools.R")) # for segment analysis
#source(file.path(segtools,"R/coor2index.R")) # coor2index

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")

## shape : alpha
## rate: beta 
## scale : 1/beta
## mean = shape * scale = alpha/beta
#' calculates gamma distribution for given parameters
#' @param x vector of x values for which the gamma distribution will
#' be calculated
#' @param start list containing gamma distribution parameters:
#' \code{a} (shape) and \code{b} (rate=1/scale)
#' @export
get_gamma <-  function(x, start) {
    a <- start$a # shape
    b <- start$b # rate = 1/scale
    b^a*x^(a-1)*exp(-x*b)/gamma(a)
}

### OPTIONS
suppressPackageStartupMessages(library(optparse))
option_list <- list(
    make_option(c("-i", "--infile"), type="character", default="", 
                help="chromosome coordinates of primary segments as produced by clusterSegments.R but without header ('allsegs.csv')"),    
    make_option(c("--chrfile"), type="character", default="",
                help="chromosome index file, providing a sorted list of chromosomes and their lengths in column 3 [default %default]"),
    ##chrfile = $YEASTDAT/chromosomes/sequenceIndex_R64-1-1_20110208.csv
    ## SEGMENT TEST SETTINGS
    make_option(c("--fuse.segs"), action="store_true", default=FALSE,
                help="use FUSE tag from clusterSegments to fuse adjacent segments"),
    ## OUTPUT OPTIONS
    make_option(c("--out.path"), type="character", default=".", 
                help="directory path for output data (figures, csv files)"),
    make_option(c("-v", "--verb"), type="integer", default=1, 
                help="0: silent, 1: main messages, 2: warning messages"),
    make_option(c("--fig.type"), type="character", default="png",
                help="figure type, png or pdf [default %default]")
)

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

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
cf <- read.table(chrfile,sep="\t",header=TRUE)
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos

## READ SEGMENTS TO BE TESTED 
if ( verb>0 )
    cat(paste("LOADING SEGMENTS:", infile, "\n"))
segs <- read.table(infile,sep="\t",header=TRUE, comment.char="")

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

## split by type
lst <- split(segs,segs$type)
sgtypes <- names(lst)

## define colors and pch for segment types
sgcols <- rainbow(length(sgtypes))
sgpchs <- rep(1:17, len=length(sgtypes)) #only 17 of 18 to avoid multiple of lty
sgltys <- rep(1:6, len=length(sgtypes))
names(sgcols) <- names(sgpchs) <- names(sgltys) <- sgtypes

## define segment classes by their type
sgcltab <- getSegmentClassTable(sgtypes,sep="_") ## TODO: use
sgclasses <- getSegmentClasses(sgtypes,sep="_")

## col, pch and lty for classes
sgclcols <- rainbow(length(sgclasses))
sgclpchs <- rep(1:17, len=length(sgclasses)) #only 17 to avoid lty multiple
sgclltys <- rep(1:6, len=length(sgclasses))
names(sgclcols) <- names(sgclpchs) <-names(sgclltys) <- sgclasses



### START ANALYSIS

## output dirs - generate locally in cwd
dir.create(out.path, showWarnings = FALSE) # length distributions

## total segment number
segnum <- unlist(lapply(lst, nrow))


### TEST SEGMENT LENGTH DISTRIBUTIONS

if ( verb>0 )
    cat(paste("Calculating segment length distribution.\n"))

## plot segment length distributions for each segment type
## store sglen summary plots!
sglen <- list() 
sgcdf <- list()
sggam <- matrix(NA, ncol=2, nrow=length(sgtypes)) # gamma distribution
rownames(sggam) <- sgtypes
colnames(sggam) <- c("a","mu") ## gamma distribution parameters

xmax <- 4.5e3
ymax <- 6e3
brks <- seq(0,35e4,100)
for ( type in sgtypes ) {
    if ( verb>0 )
        cat(paste("segment", type, "distribution\t",time(),"\n"))
    sgs <- lst[[type]]
    len <- sgs[,"end"] - sgs[,"start"] +1
    ## cumulative dist. function
    sgcdf[[type]]  <- ecdf(len)

    ## gamma dist. nls fit
    tmp <- hist(len,breaks=brks,plot=FALSE)
    mn <- mean(len)
    ## mean = scale * shape
    xy <- data.frame(x=tmp$mids, y=tmp$density)
    a <- 1.05
    start <- list(a=a, b=a/mn) #  start <- list(a=20, b=1/10) looks good, mn=200
    fit <- nls(y ~ b^a*x^(a-1)*exp(-x*b)/gamma(a), data=xy, start=start)
    a <- coefficients(fit)["a"]
    mu <- a/coefficients(fit)["b"]
    sggam[type,] <- c(a,mu)

    ## cut max. for distribution plots
    len[len>xmax] <- xmax
    ## calculate length distribution histogram
    tmp <- hist(len,breaks=brks,plot=FALSE)
    sglen[[type]]  <- tmp

    ## plot length distribution and gamma dist, and length cumulative dist.
    file.name <- file.path(out.path,paste("length_segment_",type,sep=""))
    plotdev(file.name,width=4.5,height=4.5, type=fig.type)
    par(mfcol=c(1,1), mai=c(.75,.75,.75,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot(tmp,border=sgcols[type],freq=FALSE,#ylim=c(0,ymax),
         xlim=c(0,xmax*1.05), xlab="length, bp",main=type)
    legend("topright",legend=c(nrow(sgs),paste("fuse",sum(sgs[,"fuse"]))),
           col=c(sgcols[type],NA),pch=c(sgpchs[type],NA),lty=c(sgpchs[type],NA))
    #legend("topright",legend=type)
    lines(xy$x, get_gamma(xy$x, as.list(coefficients(fit))),type="l",col=4)# GAMMA
    legend("right", paste(c("a","mu"),":", signif(sggam[type,],3)),
           lty=c(1,NA),col=4)
    dev.off()
    #high <- which(tmp$counts>ymax)
    #axis(3, at= tmp$mids[high], label=tmp$counts[high],las=2, cex.axis=.7)
    file.name <- file.path(out.path,paste("length_segment_",type,"_cum",sep=""))
    plotdev(file.name,width=4.5,height=4.5, type=fig.type)
    par(mfcol=c(1,1), mai=c(.75,.75,.75,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot(sgcdf[[type]],xlim=c(0,xmax),col=sgcols[type],main=NA,
         xlab="length, bp",ylab="cum.dist.fun.")
    dev.off()
}

## calculate average length dist for segment classes
x <- sglen[[1]]$mids

for ( i in 1:ncol(sgcltab) ) {

    scl <- colnames(sgcltab)[i]
    if ( verb>0 )
        cat(paste("class", scl, "distribution\t",time(),"\n"))

    classes <- as.character(unique(sgcltab[,i]))

    ## convert stored histograms and cumulative dist. functions (CDF)
    ## to matrices of mean and sd values for segment classes
    tot <- rep(NA, length(classes))
    means <- matrix(NA,nrow=length(x),ncol=length(classes))
    names(tot) <- colnames(means) <- classes
    cdfmeans <- cdfci <- ci <- means
    totsd <- tot
    for ( class in classes ) {
        types <- rownames(sgcltab)[as.character(sgcltab[,i])==class]
        lens <- matrix(NA,ncol=length(types),nrow=length(x))
        colnames(lens) <- types
        cdfs <- lens
        for ( type in types )  {
            lens[,type] <- sglen[[type]]$counts
            cdfs[,type] <- sgcdf[[type]](x)
        }
        ## mean distributions of this class
        means[,class] <- apply(lens,1,mean) # MEAN
        ci[,class] <- apply(lens,1,sd)      # STANDARD DEVIATION
        cdfmeans[,class] <- apply(cdfs,1,mean) # MEAN CDF
        cdfci[,class] <- apply(cdfs,1,sd)      # SD of CDF
        ## total segment number
        tot[class] <- mean(segnum[types])
        totsd[class] <- sd(segnum[types])
    }

    ## plot mean/sd of histograms and of CDF
    file.name <- file.path(out.path,paste("length_class_",scl,sep=""))
    plotdev(file.name,width=9,height=4.5,type=fig.type)
    par(mfcol=c(1,1),mai=c(.75,.75,.5,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot(0,col=NA, xlim=c(0,xmax), ylim=c(0,ymax),
         xlab="length, bp",ylab="count",main=paste("class:", scl))
    for ( class in classes ) {
        sgtype <- paste(scl,":",class,sep="")
        lines(x, means[,class],col=sgclcols[sgtype],lty=sgclltys[sgtype],lwd=2)
        px <- c(x,rev(x))
        py <- c(means[,class]-ci[,class],rev(means[,class]+ci[,class]))
        polygon(x=c(px,px[1]),y=c(py,py[1]),border=NA,
                col=sub("FF$","33",sgclcols[sgtype]))
    }
    tmp <- paste(scl,":",classes,sep="")
    leg <- paste(classes, ": ",
                 round(tot/1000,1), "+/-",
                 round(totsd/1000,1), " k", sep="")
    legend("topright",legend=leg, col=sgclcols[tmp],lty=sgclltys[tmp])
    dev.off()
    
    file.name <- file.path(out.path,paste("length_class_",scl,"_CDF",sep=""))
    plotdev(file.name,width=9,height=4.5,type=fig.type)
    par(mfcol=c(1,1),mai=c(.75,.75,.5,.1),mgp=c(1.75,.5,0),xaxs="i")
    plot(0,col=NA, xlim=c(0,xmax), ylim=c(0,1),
         xlab="length, bp",ylab="cum. dist. fun.",main=NA)
    abline(h=0:1,lty=2,col="gray")
    for ( class in classes ) {
        sgtype <- paste(scl,":",class,sep="")
        lines(x, cdfmeans[,class],col=sgclcols[sgtype],
              lty=sgclltys[sgtype],lwd=2)
        px <- c(x,rev(x))
        py <- c(cdfmeans[,class]-cdfci[,class],
                rev(cdfmeans[,class]+cdfci[,class]))
        polygon(x=c(px,px[1]),y=c(py,py[1]),border=NA,
                col=sub("FF$","33",sgclcols[sgtype]))
    }
    dev.off()
}

sgtab <- table(segs[,"type"]) 
file.name <- file.path(out.path,paste("segmentNumber",sep=""))
plotdev(file.name,width=10,height=5,type=fig.type,res=300)
par(mai=c(1.8,.75,.15,.1),mgp=c(1.75,.5,0))
plot(sgtab,type="h",col=sgcols,xlab=NA,ylab="# of segments",axes=FALSE)
axis(2)
axis(1,at=1:length(sgtab),labels=names(sgtab),las=2,cex.axis=.5)
dev.off()

## write out table with number of segments and gamma distribution fit parameters
result <- data.frame(ID=names(segnum),segnum,sggam)
file.name <- file.path(out.path,paste("segmentLengths.csv",sep=""))
write.table(result,file=file.name, sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

## TODO: more analyses of segment number vs. classes
## TODO: segment length vs. dynamics - load result
## from analyzeSegments.R
