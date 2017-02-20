#!/usr/bin/env Rscript

## cut time-series into primary segments based on level
## of expression

library("segmenTools") # coor2index and presegment function
library("stringr") # for 0-padded filenames

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")

### OPTIONS
suppressPackageStartupMessages(library(optparse))
option_list <- list(
    ## INPUT
    make_option(c("-i", "--infile"), type="character", default="", 
                help="RData file containing matrix 'ts' and chromosome coordinates 'coor'"),    
    ## PRIMARY SEGMENT SETTINGS
    make_option(c("--avg"), type="integer", default=1000,  
                help="window length for moving average of read-count presence [default: %default]"),
    make_option(c("--favg"), type="integer", default=100,  
                help="window length for primary segment extension [default: %default]"),
    make_option(c("--minrd"), default=8,
                help="minimal count of reads in main prim.seg. cutoff [default %default]"),
    make_option(c("--minds"), type="integer", default=5000,  
                help="minimal distance to be treated as separate segments [default: %default]"),
    make_option(c("--minsg"), type="integer", default=5000,  
                help="minimal segment length for fusion [default: %default]"),
    make_option(c("--rmlen"), type="integer", default=1,  
                help="minimal segment length for removal [default: %default]"),
    ## OUTPUT OPTIONS
    make_option(c("-o", "--outdir"), type="character", default=".", 
                help="directory path for output data (figures, csv files"),
    make_option(c("-v", "--verb"), type="integer", default=1, 
                help="0: silent, 1: main messages, 2: warning messages"),
    make_option(c("--fig.type"), type="character", default="png",
                help="figure type, png or pdf [default %default]"),
    make_option(c("--plot.borders"), action="store_true", default=FALSE,
                help="plot borders between primary segments"),
    make_option(c("--plot.summary"), action="store_true", default=FALSE,
                help="plot length distribution of primary segments "),
    make_option(c("--write.segments"),  action="store_true", default=FALSE,
                help="write out data for each segment [default %default]")
)

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))


## promote options to main environment and print all arguments
if ( opt$verb>0 )
    cat(paste("SETTINGS:\n"))
for ( i in 1:length(opt) ) {
    if ( opt$verb>0 )
        cat(paste(names(opt)[i], "\t", #typeof(opt[[i]]),
                  paste(opt[[i]],collapse=", "), "\n",sep=""))
    arg <- names(opt)[i]
    assign(arg, opt[[arg]])
}
if ( verb>0 )
    cat(paste("\n"))


## load time-series and oscillation data: coor, ts
if ( verb>0 )
    cat(paste("Loading data:", infile, "\n"))
load(infile)

## generate chromosome mapping index
## required to map between data row and chromosome position
## for faster indexing
chr <- unique(coor[,"chr"])
chrS <- rep(0,length(chr)+1)
for ( i in 2:length(chrS) ) {
    idx <- which(coor[,"chr"]==chr[i-1])
    chrS[i] <- idx[length(idx)] - nrow(coor)/2
}


## PRE-SEGMENT TIME-SERIES

## pass fig.path to induce border plotting (in dir fig.path)
## and seg.path to write files with data for all segments
if ( verb>0 )
    cat(paste("Calculating pre-segmentation\t",time(),"\n",sep=""))
if ( plot.borders ) {
    if ( verb>0 ) cat(paste("... and plotting border scans.\n"))
    primseg <- presegment(ts=ts, chrS=chrS, map2chrom=FALSE,
                          avg=avg, favg=favg,
                          minrd=minrd, minds=minds, minsg=minsg, rmlen=rmlen,
                          fig.path=outdir, verb=verb)
} else {
    primseg <- presegment(ts=ts, chrS=chrS, map2chrom=FALSE,
                          avg=avg, favg=favg,
                          minrd=minrd, minds=minds, minsg=minsg, rmlen=rmlen,
                          verb=verb)
}

## column chr will be filled by index2coor
primseg <- data.frame(ID=paste("sg",str_pad(1:nrow(primseg),4,pad="0"),sep=""),
                      chr=rep(NA,nrow(primseg)),primseg)

## write out segments!
if ( write.segments ) {
    if ( verb>0 )
        cat(paste("Writing segment data files\t",time(),"\n",sep=""))
    writeSegments(data=ts, segments=primseg, name="primseg", path=outdir)
}

## inter-segments
emptyseg <- data.frame(ID=paste("is",str_pad(1:(nrow(primseg)+1),4,pad="0"),
                         sep=""),
                       chr=rep(NA,nrow(primseg)+1),
                       start=c(1,primseg[1:nrow(primseg),"end"]+1),
                       end=c(primseg[1:nrow(primseg),"start"]-1,nrow(ts)))
emlen <- emptyseg[,"end"] - emptyseg[,"start"] +1 
emptyseg <- emptyseg[emlen>0,] # rm 0-length
## TODO: remove >max(chrS)

## write out inter-segments!
if ( write.segments ) {
    if ( verb>0 )
        cat(paste("Writing inter-segment data files\t",time(),"\n",sep=""))
    writeSegments(data=ts, segments=emptyseg, name="primseg", path=outdir)
}

## TODO: write out presegmentation stats
## total coverage?

## (5): map back to chromosome coordinates
## and write to file
if ( verb>0 )
    cat(paste("Writing segments and intersegments\t",time(),"\n",sep=""))

sgcoors <- index2coor(primseg,chrS)
file.name <- file.path(outdir,"primseg.csv")
write.table(sgcoors,file.name,row.names=FALSE,sep="\t",quote=FALSE)

emcoors <- index2coor(emptyseg,chrS)
file.name <- file.path(outdir,"primseg_interseg.csv")
write.table(emcoors,file.name,row.names=FALSE,sep="\t",quote=FALSE)


### PLOTTING
if ( !plot.summary )
    quit(save="no")
if ( verb>0 )
    cat(paste("Summary plot\t",time(),"\n",sep=""))

## calculate numts as in function presegment
numts <- rowSums(ts > 0) ## timepoints with reads


## (5) analyze segment length and read-count presence distributions
## primary segment lengths
sglen <- primseg[,"end"] - primseg[,"start"] +1
emlen <- emptyseg[,"end"] - emptyseg[,"start"] +1 

## genome coverage
total <- max(chrS)*2
sgcvg <- sum(sglen)/total
emcvg <- sum(emlen)/total



## average number of expressed time-points in emptyseg
emexpr <- apply(emptyseg,1,function(x)
    sum(numts[x["start"]:x["end"]])/length(x["start"]:x["end"]))
sgexpr <- apply(primseg,1,function(x)
    sum(numts[x["start"]:x["end"]])/length(x["start"]:x["end"]))


max.sg <- max(c(sglen))
max.em <- max(c(emlen))
mx <- max(c(max.sg,max.em)) *1.01
xl<- 2.5*mean(sglen)

## plot segments and inter-segments
## average time-point presence and segment lengths
file.name <- file.path(outdir,"primseg_lengths")
plotdev(file.name,width=5,height=4,type=fig.type)
par(mfcol=c(2,1),mai=c(.5,.75,.15,.1),mgp=c(1.5,.5,0))
hist(emexpr,breaks=seq(0,24,.5),border=2,xlab="mean number of expressed time points",
     main=NA)#"mean number of expressed time points")
hist(sgexpr,breaks=seq(0,24,.5),add=TRUE)
legend("right",legend=c("pre-segments","inter-segments"),
       col=1:2, pch=15,bty="n")
hist(emlen,breaks=seq(0,mx,500),border=2,xlim=c(0,xl),xlab="length, bp",
     main=NA)#"segment length distribution")
hist(sglen,breaks=seq(0,mx,500),add=TRUE)
legend("right",legend=c(paste(length(sglen), "segments"),
                        paste(100*c(round(sgcvg,2)),"% coverage",sep=""),
                        paste("max:", round(c(max.sg)/1e3), "kb"),
                        paste("tail: ", sum(sglen>xl),sep="")),
       col=c(1,NA,NA,NA), pch=15,bty="n")       
dev.off()
cat(paste("maximum segment, bp\t", max.sg, "\n"))
cat(paste("maximum inter-segment, bp\t", max.em, "\n"))

## time is N^2
sg.eta <- ecdf(sglen^2)
em.eta <- ecdf(emlen^2)
x <- seq(0,sum(sglen^2),1e4)
mxx <- x[which(sg.eta(x)==1)][1]
x <- seq(0,mxx,mxx/1000) 
file.name <- file.path(outdir,"primseg_times")
plotdev(file.name,width=5,height=2.5,type=fig.type)
par(mfcol=c(1,1),mai=c(.5,.75,.05,.05),mgp=c(1.5,.5,0))
plot(x, sg.eta(x),type="l",log="x",
     xlab="time ~ N^2, bp^2", ylab="E.C.D.F.", ylim=c(0,1))
lines(x, em.eta(x), col=2)
abline(v=mxx)
legend("right",paste("max", signif(mxx,3)))
dev.off()
