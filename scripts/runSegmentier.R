#!/usr/bin/env Rscript


library("segmenTools") # coor2index
library("segmenTier")  # main segmentation algorithm
suppressPackageStartupMessages(library(optparse))
#library("Rcpp")
#source("~/programs/segmenTier/R/cluster.R")
#source("~/programs/segmenTier/R/segment.R")
#sourceCpp("~/programs/segmenTier/src/cluster.cpp")
#sourceCpp("~/programs/segmenTier/src/segment.cpp")

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=stdout()) # until piping is implemented

### PARSE OPTIONS

option_list <- list(
    ## INPUT OPTIONS
    make_option(c("-i", "--infile"), type="character", default="", 
                help="chromosome coordinates of primary segments"),    
    make_option(c("--chrfile"), type="character", default="",
                help="chromosome index file, providing a sorted list of chromosomes and their lengths in column 3 [default %default]"),
    make_option(c("--datafile"), type="character", default="",
                help="full data set, RData that contains the time series in matrix 'ts' and its genomic coordinates in matrix 'coor' [default %default]"),
    make_option(c("--primdir"), type="character", default="primarysegments",
                help="individual data sets for primary segments [default %default]"),
    make_option(c("--primfiles"), type="character", default="primseg_",
                help="individual data sets for primary segments [default %default]"),    
    make_option(c("--use.data"), action="store_true", default=FALSE,
                help="use the full data set (helper, will be obsolete)"),
    ## OUTPUT OPTIONS
    make_option(c("-v", "--verb"), type="integer", default=1, 
                help="0: silent, 1: main messages, 2: warning messages"),
    make_option(c("--fig.type"), type="character", default="png",
                help="figure type, png or pdf [default %default]"),
    make_option(c("--do.test"), action="store_true", default=FALSE,
                help="only do testsets"),
    make_option(c("--only.plot"), action="store_true", default=FALSE,
                help="only plot existing segmentations; useful to check ongoing results during longer runs"),
    make_option(c("--redo"), action="store_true", default=FALSE,
                help="overwrite existing segmentations; leave FALSE if you merely want to plot existing segmentations"),
    make_option(c("--save.matrix"), action="store_true", default=FALSE,
                help="NOT IMPLEMENTED: store the total score matrix S(i,c) and the backtracing matrix K(i,c) in the RData file"),
    make_option(c("--plot"), action="store_true", default=FALSE,
                help="plot segment figures"),
    make_option(c("--genbro"), type="character", default="",
                help="path for the R genome browser sources"),
    make_option(c("--gendat"), type="character", default="",
                help="path for the R genome data sources"),
    make_option(c("--idsuffix"), type="character", default="",
                help="suffix for segment IDs [default %default]"),
    make_option(c("--type.name"),  type="character", default=c("T,D"),
                help="parameters to use for IDs [default %default]"),
    make_option(c("--short.name"),  action="store_true", default=FALSE,
                help="strip down segment types to varied values [default %default]"),
    make_option(c("--genome"), type="character", default="yeast_R64-1-1",
                help="in plots, additionally show annoted features and reported transcripts; this requires installation of the `R Genome Browser' and data sets and is currently only available for `yeast_R64-1-1"),
    make_option(c("-o", "--outdir"), type="character", default=".", 
                help="path to output directory"),    
    ## PRIMARY SEGMENT SELECTION
    make_option(c("--segs"), type="character", default="",  
                help="primary segment numbers, all or test if empty, comma-separated list [default: %default]"),
    make_option(c("--use.empty"), action="store_true", default=FALSE,
                help="cluster inter-segment regions instead"),
    make_option(c("--min.sz"), type="integer", default=2,  
                help="minimal size of segments (non cluster-0 values) [default: %default]"),
    make_option(c("--max.sz"), type="integer", default=3e6,  
                help="maximal size of segments (non cluster-0 values) [default: %default]"),
    ## TIMESERIES PROCESSING
    make_option(c("--trafo"), type="character", default="raw",
                help="time-series transformation function, R base functions like 'log', and 'ash' for asinh is available [default %default]"),
    make_option("--low.thresh", default=-Inf,
                help="cut-off to set low values to 0 [default %default]"),
    ## DISCRETE FOURIER TRANSFORM
    make_option(c("--use.fft"), action="store_false", default=TRUE, 
                help="do DFT of time-series [default %default]"),
    make_option(c("--dft.range"), type="character", default="1,2,3,4,5,6,7", 
                help="DFT components to use, comma-separated list of integers [default %default]"),
    make_option(c("--use.snr"), action="store_false", default=TRUE,
                help="do SNR of the DFT [default %default]"),
    make_option(c("--dc.trafo"), type="character", default="raw", 
                help="DC component transformation function, see --trafo [default %default]"),
    ## CLUSTERING K-MEANS
    make_option(c("--K"), type="character", default="12,16,20", 
                help="number of clusters to use in k-means, comma-separated [default %default]"),
    make_option(c("--iter.max"), type="integer", default=100000,
                help="max. iterations in k-means [default %default]"),
    make_option(c("--nstart"), type="integer", default=100,
                help="initial k-means configurations [default %default]"),
    make_option("--nui.thresh", default=0.6,
                help="cluster correlation threshold to re-assign points to nuissance [default %default]"),
    ## SEGMENTATION PARAMETERS
    make_option("--nui.cr", type="character", default="2",
                help="correlation of nuissance cluster with others (-) and itself (+), comma-separated characters [default %default]"),
    make_option(c("--scores"), type="character", default="ccor,icor", 
                help="scoring functions to use, comma-separated characters [default %default]"),
    make_option(c("--scales"), type="character", default="1,3,5", 
                help="scale exponent for similarity matrices, comma-separated doubles [default %default]"),
    make_option(c("--M"), type="character", default="175", 
                help="minimal length penalty, comma-separated integers [default %default]"),
    make_option(c("--Mn"), type="character", default="15", 
                help="minimal length penalty [default %default]"),
    make_option(c("--multi"), type="character", default="max", 
                help="handling of multiple max. score k in scoring, min or max [default %default]"),
    make_option(c("--multib"), type="character", default="max", 
                help="handling of multiple max. score k in back-tracing, min or max [default %default]"),
    make_option(c("--nextmax"), type="character", default="T", 
                help="in back-tracing, search for the next non-decreasing S(i,c) before initializing a new segment [default %default]"),
    ## POST-PROCESSING
    make_option("--fuse.thresh", default=0.2,
                help = "correlation threshold between clusters to fuse adjacent segments [default \"%default\"]")
)

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## process comma-separated list arguments
lst.args <- c(segs="integer",
              type.name="character",
              dft.range="integer",
              K="integer",
              scores="character",
              scales="numeric",
              M="integer", Mn="integer", nui.cr="integer",
              nextmax="logical",
              multi="character",multib="character")
for ( i in 1:length(lst.args) ) {
    idx <- which(names(opt)==names(lst.args)[i])
    ## get individual values
    tmp <- as.list(unlist(strsplit(opt[[idx]], ",")))
    ## expand ranges
    if ( lst.args[i]=="numeric" |  lst.args[i]=="integer" )
        for ( j in 1:length(tmp) ) { # only for numeric modes
            tmp2 <- unlist(strsplit(tmp[[j]], ":"))
            if ( length(tmp2)>1 ) {
                tmp2 <- as.numeric(tmp2)
                tmp[[j]] <- tmp2[1]:tmp2[2]
            }
        }
    if ( length(tmp)>0 )
        opt[[idx]] <- unlist(tmp)
    mode(opt[[idx]]) <- lst.args[i]
}


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


### START
cat(paste("LOADING SEQUENCING DATA\t", time(), "\n",sep=""))

## TODO: move this out and pre-write primary segment data
## via a separate function writing out primarySegments
## or use the infile approach?

## load time-series and oscillation data: coor, ts, osc
## only "ts", the original time-series and "coor", the corresponding
## chromosome coordinates are used here!
if ( chrfile != "" ) {
    cat(paste("Using chromosome index file:",chrfile, "\n"))
    cf <- read.table(chrfile, sep="\t",header=FALSE)
    chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos
} else {
    ##datafile <- file.path("/data/yeast/RNAseq/results/genomeData/Sacchromyces.cerevisiae.clean.S288C_oscillation.RData")
    cat(paste("Using full data file:", datafile, "\n"))
    load(datafile) 
    
    chr <- unique(coor[,"chr"])
    chrS <- rep(0,length(chr)+1)
    for ( i in 2:length(chrS) ) {
        idx <- which(coor[,"chr"]==chr[i-1])
        chrS[i] <- idx[length(idx)] - nrow(coor)/2
    }
} 

## NOT USED:
#rawfile <- file.path(data.path,"genomeData",
#                    "Sacchromyces.cerevisiae.clean.S288C_timeseries.csv")


## LOAD SELECTED PRIMARY SEGMENTS
cat(paste("Loading primary segments:",infile,"\n"))

primseg <- read.table(infile, sep="\t",header=TRUE)
primseg <- coor2index(primseg, chrS)[,c("ID","start","end")]
prdf <- primseg[,"end"] - primseg[,"start"] + 1
## set segment name (used in segment IDs and file names)
outname <- file.path(outdir,"primseg")


## TESTSETS

## primseg3 tests
tests <- list(known=c(srg1=436,gdh3=4),
              nonsplit=c(3133,16,39,100,2360,2427),
              long5utr=c(34),
              short3utr=c(38),
              splitgene=c(37),
              fiveprimedeg=c(4,14,21,25))
alltests <- unique(unlist(tests)) # new primseg3 20161115

### SELECT SEGMENTS
if ( length(segs)==0 ) {
    sets <- order(prdf) # rev(order(prdf)) ) #
    if ( do.test ) {
        cat(paste("DEVEL OPTION: CALCULATING TESTSETS\n"))
        sets <- alltests
    }
} else { ## command-line
    sets <- segs
} 

if ( !only.plot )
  cat(paste("CALCULATING SEGMENTATIONS\t", time(), "\n",sep=""))
cat(paste("TESTSETS\t", length(sets), "\n",sep=""))

do.sets <- sets
if ( only.plot ) { # skip segmentation; continue at plots
    do.sets <- c()
    plot <- TRUE
}

### RUN SEGMENTATION
for ( i in do.sets ) { 

    ## generate segment id
    segdat <- i
    segid <- as.character(primseg[i,"ID"])
    ##segid <- str_pad(i,5,pad="0")
    ##if ( idsuffix!="" )
    ##    segid <- paste(segid, idsuffix, sep="_")
    
    cat(paste("PRIMARY SEGMENT\t", segid, "\t", which(sets==i), "of",
              length(sets),"\n",sep=""))

    rng <- primseg[i,"start"]:primseg[i,"end"]
    strand <- idx2str(primseg[i,"start"],chrS)
    
    ## minimial size of segment in terms of expressed points!
    if ( length(rng) < 2 ) {
        cat(paste("\tinvalid primary segment\n"))
        next
    }
    ## minimial size of segment in terms of expressed points!
    if ( length(rng) < min.sz ) {
        cat(paste("\ttoo short:",length(rng),"\n"))
        next
    }
    ## maximal size of segment in terms of expressed points!
    if ( length(rng) > max.sz ) {
        cat(paste("\ttoo long:",length(rng),"\n"))
        next
    }
    ## already done?
    file.name <- file.path(paste(outname,"_",segid,"_clusters.csv",sep=""))
    if ( file.exists(file.name) & !redo ) {
        cat(paste("\talready done\n"))
        next
    }
    
    ## load time-series
    #ts <- read.table(rawfile,sep="\t",header=FALSE,colClasses="numeric",
    #                 skip=primseg[i,"start"],
    #                 nrows=primseg[i,"end"]-primseg[i,"start"]+1)
    if ( !use.data ) {
        segdata <- file.path(primdir, paste(primfiles,segid,".csv",sep=""))
        tsd <- read.table(segdata, sep="\t", header=TRUE)
    } else 
        tsd <- ts[rng,]

    ## reverse data from reverse strand!
    ## both strands are handled in direction of transcription!
    if ( strand==-1 )
        tsd <- tsd[nrow(tsd):1,]

    cat(paste("CLUSTERING\t",time(),"\n",sep=""))

    ## process time series, get DFT etc.
    ## TODO: add processing info to tset
    ##       allow multiple processing, and take over IDs in cluster
    tset <- processTimeseries(ts=tsd,
                              trafo=trafo, dc.trafo=dc.trafo,
                              use.fft=use.fft, dft.range=dft.range,
                              use.snr=use.snr, low.thresh=low.thresh)
    ## cluster time series
    cset <- clusterTimeseries(tset,K=K, iter.max=iter.max, nstart=nstart,
                              nui.thresh=nui.thresh, verb=verb)
    
    ## segment all clusterings for different scoring functions
    allsegs <- NULL
    if ( !is.null(cset) ) {

        ## collect varysettings
        vary <- list(## scoring function
            E=scales, S=scores,
            ## scoring params
            M=M, Mn=Mn, a=2, nui=nui.cr,
            ## backtrace params
            nextmax=nextmax, multi=multi,multib=multib)

        ## note: we add time-series processing to the segment
        ## type names, allowing to externally vary (multiple
        ## distinct calls to this script) and later co-analyse
        ## them.
        sset <- segmentCluster.batch(cset, varySettings=vary, 
                                     verb=verb,
                                     fuse.threshold=fuse.thresh,
                                     id=segid, type.name=type.name,
                                        #short.name=short.name,
                                     save.matrix=save.matrix)
        allsegs <- sset$segments # SEGMENTS!
        ## plotSegmentation(tset,cset,sset)

        if ( is.null(allsegs) )
          cat(paste("no segments\n"))
        else {
            ## CHROMOSOMAL COORDINATES
            ## 1) add/subtract to primseg indexed coordinates
            if ( strand==-1 ) {
                segcoors <- cbind(chr = rep(NA, nrow(allsegs)),
                                  primseg[i,"end"] + 1 -
                                  allsegs[, c("start","end"), drop=FALSE],
                                  strand=rep(NA,nrow(allsegs)))
            } else {
                segcoors <- cbind(chr = rep(NA, nrow(allsegs)),
                                  allsegs[, c("start","end"), drop=FALSE] +
                                  primseg[i,"start"] - 1,
                                  strand=rep(NA,nrow(allsegs)))
            }
            ## 2) map back to chromosomes
            segcoors <- index2coor(segcoors,chrS)
            allsegs <- cbind(allsegs[,c("ID","type")],
                             segcoors,
                             allsegs[,c("CL","color","fuse"),drop=FALSE])
            
            ## store clustering data
            centers <- cset$centers
            clusters <- cset$clusters
            if ( strand == -1 ) 
                clusters <- clusters[nrow(clusters):1,]
            
            ## save data!
            file.name <- file.path(paste(outname,"_",segid,sep=""))

            ## as table
            write.table(allsegs,sep="\t",
                        file=paste(file.name,"_segments.csv",sep=""),
                        quote=FALSE,row.names=FALSE,col.names=TRUE)
            ## as RData, mainly for plots below
            sset$segments <- allsegs
            save(opt, segid, tset, cset, sset, 
                 file=paste(file.name,"_segments.RData",sep=""))
        }
    }
}

if ( !plot ) next

cat(paste("PLOTTING\t",time(),"\n",sep=""))

## TODO: load plotting specific parameters
##       and align with opt from saved RData file

columns <- c(name="name", chr="chr", strand="strand",
             start="start", end="end", type="type", color="color")

## LOAD YEAST GENOME BROWSER
##genome <- "yeast_R64-1-1"
if ( genome=="yeast_R64-1-1" ) {

### REQUIRES ENVIRONMENT VARIABLES SET TO YEAST GENOME BROWSER & DATA!
### TODO: make genomeBrowser a package; introduce official IDs for genomes
### see https://gitlab.com/raim/genomeBrowser and raim@tbi.univie.ac.at
### for the required data files

    cat(paste("genome data\t", genome, "\n",sep=""))

    ## for genome data and plots - todo - skip this!
    if ( missing(genbro) )
      genbro <- Sys.getenv("GENBRO")
    source(file.path(genbro,"src/genomeBrowser.R")) ## for loadData
    source(file.path(genbro,"src/genomeBrowser_utils.R")) ## plotFeature
    if ( missing(gendat) )
      gendat <- Sys.getenv("GENDAT")
    gendat <- file.path(gendat,"yeast")

    
    ## TODO: load data from files to make independent of
    ## genomeBrowser/genomeData
    ## PLOT SETTINGS FOR TESTSETS
    fcolumns <- columns
    fcolumns["color"] <- "CL_rdx_col"
    ftypes <- c( "gene_cassette"       , "gene",           
                "five_prime_UTR_intron", "intron",                   
                "dubious"              , "pseudogene",               
                "rRNA"                 , "tRNA",              
                "tRNA.intron"          , "ncRNA",                    
                "snRNA"                , "snoRNA",                   
                "telomere"             , "centromere",               
                "ARS"                  , "transposable_element_gene",
                "LTR_retrotransposon"  , "repeat_region")
    
    ## LOAD DATA SETS
    cat(paste("Loading annotation data\n"))
    testIDs <- c("transcripts","annotation")
    dataSets <- loadData(testIDs, data.path=gendat)
    
    dataSets[["annotation"]]$settings$names <- FALSE
}


### PLOTTING

for ( i in sets ) { 

    ## generate segment id
    segid <- primseg[i,"ID"]
    ##segid <- str_pad(i,5,pad="0")
    ##if ( idsuffix!="" )
    ##    segid <- paste(segid, idsuffix, sep="_")
    
    cat(paste("Primary segment\t", segid, "\t", which(sets==i), "of",
              length(sets),"\n",sep=""))

    file.name <- file.path(paste(outname,"_",segid,sep=""))
 
    ## already plotted?
    if ( file.exists(paste(file.name,fig.type,sep=".")) & !redo ) {
        cat(paste("\talready plotted.\n"))
        next
    }

    ## load data
    dfile <- paste(file.name,"_segments.RData",sep="")
    if ( !file.exists(dfile) ) {
        cat(paste("\tnot yet calculated.\n"))
        next
    }
    load(dfile)

    ## record settings of plotted data
    ## TODO: write setting file  before segmentation run and
    ## use to load into analysis scripts?
    sink(paste(file.name,"_settings.dat",sep=""))
    if ( opt$verb>0 )
        cat(paste("LOADED SETTINGS:\n"))
    for ( j in 1:length(opt) ) {
        if ( opt$verb>0 )
            cat(paste(names(opt)[j], "\t", #typeof(opt[[i]]),
                      paste(opt[[j]],collapse=", "), "\n",sep=""))
        ##arg <- names(opt)[i]
        ##assign(arg, opt[[arg]])
    }
    sink()

    ## plot
    ## TODO:
    ## consolidate to plot functions, aligned with segment_data.R
    ## plot.tset(tset)
    ## plot.cset(cset, k)
    ## plot.sset(sset, k)
    ## plotAll(tset,cset,sset) # plot by k in fitting order
    coors <- index2coor(t(c(chr=1,unlist(primseg[i,c("start","end")]))),chrS)
    strand <- ifelse(coors[,"strand"]==-1, "-", "+")
    xaxis <- coors[,"start"]:coors[,"end"]
    N <- nrow(tset$ts)

    ## revert back timeseries!
    ## TODO: use function or move to plot
    if ( strand=="-" ) {
        tset$tot <- rev(tset$tot)
        tset$zero.vals <- rev(tset$zero.vals)
        tset$ts <- tset$ts[nrow(tset$ts):1,]
    }
    
    ## adapt width to segment length!
    ## TODO: adapt left mai according to longest type name!
    width <- 2.5 + N/1e3 # 1 kb per inch; plut left margin
    if ( fig.type=="png" )
      width <- min(c(width, 326)) # APPARENTLY THE MAX. WIDTH ALLOWED
    nsg <- length(sset$ids)
    nrows <- 3 + ncol(cset$clusters) + ifelse(genome=="yeast_R64-1-1",2,0)
    if ( save.matrix )
        nrows <- nrows + length(sset$SK)
    heights <- rep(0.7,nrows)
    heights[4] <- max(c(0.7, nsg * 0.7/10))
    height <- sum(heights)

    plotdev(file.name,width=width,height=height,type=fig.type)
    ## TODO: replace mfcol by layout and adjust height
    ## with segment number
    par(mai=c(.01,1.25,.01,.01),mgp=c(1.7,.5,0),xaxs="i")
    layout(matrix(1:nrows),widths=1,heights=heights)
    plot(tset,plot=c("total","timeseries"),xaxis=xaxis)
    for ( k in 1:ncol(cset$clusters) )
      plot(cset, k=k)
    plot(sset,plot="segments",xaxis=xaxis, lwd=2)
    axis(1)
    if ( genome=="yeast_R64-1-1" ) {
        tmp <- plotFeatures(dataSets[["annotation"]]$data,
                            coors=coors, strand=strand,
                            typord=TRUE, cuttypes=TRUE,
                            axis1=TRUE, ylab=NA, names=TRUE,
                            columns=fcolumns, types=ftypes)
                                        #axis(1)
        tpy <- plotFeatures(dataSets[["transcripts"]]$data,
                            coors=coors, strand=strand,
                            typord=TRUE, cuttypes=TRUE, ylab=NA)
    }
    ## PLOT SCORING MATRICES
    if ( save.matrix ) 
      plot(sset,plot="S",xaxis=xaxis)
    dev.off()
 
}

cat(paste("DONE AT  \t",time(),"\n",sep=""))
quit(save="no")

## OLD TODO - keep until checked
segmentData <- function(segment=primseg[1,],data,kmax=30) {
  ### cluster DFT 
  ## get data
  ## flowClust, Kmax=30, find best K via BIC, or use cluster DAG/TREE
  ## do on command-line or convert genomeOscillation.R to functions
  ## see clusterSegments.sh
  clsk <- kmeans(data,centers,1:kmax)
  ## calculate cluster similarity
  ## run SegmentDP - use clustering or directly
  ##                 use similarity to cluster medians!?
  ##                 or BEST: cluster probability matrix
  ## back-trace winning clusters for each sub-segment
  ## calculate mean for winning clusters
  ## return: <sub-segment, cluster, segmentmedian>
  ##         <cluster, clustermedian> 
}


### 3) global clustering

clusterSegments <- function(secseg) {
  ## collect all secseg clustermedians
  ## cluster all clustermedians
  ## re-assign each cluster to global cluster
  ##  OR simpler
  ## cluster segmentmedians
  ## return: <sub-segment, cluster, segmentmedian>
  ##         <cluster, datamedian>
}


