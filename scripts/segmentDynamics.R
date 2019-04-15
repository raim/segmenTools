#!/usr/bin/env Rscript

## COLLECTING SEGMENT READ DISTRIBUTIONS AND CLUSTERING OF SEGMENT TIME SERIES

## load segments and time-series, incl. phase/pval
## * analyze phase-dist of adjacent segments, tandem
## select significantly different segments
## * rm all coding segments and calculate overlap, Jaccard measure and
## p-vals for remaining segments vs. CUT/SUT/XUT
## * rm all CUT/SUT/XUT segments and analyze additional
## * compare to ARS
## * compare to introns
## * analyze genes non-consistent with clustering: oscillating n/r, vs.
## wrong or none-oscillating major cluster genes; C vs. D in last cycle

## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")


## segment utils

library("segmenTier") # for processTimeseries - TODO : unify coor2index!
library(segmenTools)
#segtools <- "~/programs/segmenTools/"
#source(file.path(segtools,"R/segmenTools.R")) # for segment analysis
#source(file.path(segtools,"R/coor2index.R")) # coor2index


#suppressPackageStartupMessages(library("stringr")) # for 0-padded filenames
suppressPackageStartupMessages(library(optparse))

### OPTIONS
option_list <- list(
  make_option(c("-i", "--infile"), type="character", default="", 
              help="chromosome coordinates of primary segments as produced by clusterSegments.R but without header ('allsegs.csv')"),    
  make_option(c("--chrfile"), type="character", default="",
              help="chromosome index file, providing a sorted list of chromosomes and their lengths in column 3 [default %default]"),
  ##chrfile = $YEASTDAT/chromosomes/sequenceIndex_R64-1-1_20110208.csv
  make_option(c("--datafile"), type="character", default="",
              help="full data set, RData that contains the time series in matrix 'ts' and its genomic coordinates in matrix 'coor' [default %default]"),
  ## SEGMENT SETTINGS
  make_option(c("--fuse.segs"), action="store_true", default=FALSE,
              help="use FUSE tag from clusterSegments to fuse adjacent segments"),
  make_option(c("--typecol"), type="character", default="type", 
              help="name of the column with segment types"),
  make_option(c("--stypes"), type="character", default="", 
              help="sub-set of segments in column 'type', option --typecol, use --typecol ALL and --stypes ALL to avoid splitting into types (unless `ALL' is an actual colum name)"),
  make_option(c("--idcol"), type="character", default="ID", 
              help="name of the column with unique segment names"),
  ## OSCILLATION SETTINGS
  make_option("--pval.thresh", default=1,
              help="phases above this thresh will be set to NA [default %default]"),
  make_option("--phase.weight", action="store_true", default=FALSE,
              help="use weight for calculating average phases of segments: 1-p.value [default %default]"),
  make_option("--pval.thresh.sig", default=1,
              help="pvals above will be counted as non-significant; optionally used in segment averaging (option --filter.reads), and segment filtering before clustering (options --with.rain or --with.permutation, and --cl.filter) [default %default]"),
  make_option("--read.rng", type="character", default="",
              help="range of time-points for total read-count, Fourier and rain calculations (but not used for clustering!), comma-separated list of integers"),
  make_option("--period", default=24,
              help="period for `rain' oscillation stats [default %default]"),
  make_option("--deltat", default=2,
              help="sampling interval for `rain' oscillation stats [default %default]"),
  ## SEGMENT AVERAGING
  make_option(c("--endcut"), type="integer", default=0, 
              help="fraction at segment ends that will not be considered for average time series"),
  make_option("--avgk", type="integer", default=0,
              help="integer width of running median smoothing window, using stats::runmed; must be odd"),
  make_option(c("--mean.ratio"), action="store_true", default=FALSE,
              help="take mean ratio of each read time-series before calculating segment average"),
  make_option(c("--filter.reads"), action="store_true", default=FALSE,
              help="use only reads with oscillation p.value < pval.thresh.sig for segment average caculation"),
  ## SEGMENT TIME-SERIES PROCESSING 
  make_option(c("--trafo"), type="character", default="raw",
              help="time-series transformation function, R base functions like 'log', and 'ash' for asinh is available [default %default]"),
  make_option(c("--use.snr"), action="store_true", default=FALSE,
              help="do SNR scaling of amplitudes [default %default]"),
  make_option(c("--dc.trafo"), type="character", default="raw", 
              help="DC component transformation function, see --trafo [default %default]"),
  make_option("--perm", type="integer", default=0,
              help="number of permutations used to calculate p-values for all DFT components"),
  make_option("--smooth.time", type="integer", default=1, # best for clustering was 3
              help="span of the moving average for smoothing of individual read time-series"),
  ## SEGMENT CLUSTERING
  make_option(c("--missing"), action="store_true", default=FALSE,
              help="only calculate missing clusterings; useful if SGE jobs were not successful, to only calculate the missing"),
  make_option(c("--dft.range"), type="character", default="2,3,4,5,6,7", 
              help="DFT components to use for clustering, comma-separated [default %default]"),
  make_option(c("--cl.filter"), type="character", default="unsig.p", 
              help="type of segment filter for clustering [default %default]"),
  
  make_option(c("--with.rain"), type="character", default="", 
              help="path for prior calculation of `rain' p-values for segment filtering, p < pval.thresh.sig [default %default]"),
  make_option(c("--with.permutation"), type="character", default="", 
              help="path for prior calculation of `permutation' p-values for segment filtering, p < pval.thresh.sig [default %default]"),
  ##    make_option("--smooth.time.plot", type="integer", default=3, # so far best! 
  ##                help="as smooth.time but only for plotting clusters not used of analysis"),
  ## FLOWCLUST PARAMETERS
  make_option("--ncpu", type="integer", default=1, 
              help="number of available cores for flowClust"),
  ##    make_option("--seed", type="integer", default=1, 
  ##                help="seed for the random number generator before calling flowclust to get stable clustering results (TODO: set.seed, does it work for flowclust?)"),
  make_option(c("--K"), type="character", default="12:20", 
              help="number of clusters to use in flowClust, comma-separated list of integers and colon-separated ranges [default %default]"),
  make_option(c("--fixedK"), type="integer", default=0, 
              help="fixed number of clusters to select in flowClust, flowMerge will start from there [default %default]"),
  make_option("--B", type="integer", default=500, 
              help="maximal number of EM iterations of flowClust"),
  make_option("--tol", default=1e-5, 
              help="tolerance for EM convergence in flowClust"),
  make_option("--lambda", default=1, 
              help="intial Box-Cox transformation parameter in flowClust"),
  make_option("--nu", default=4, 
              help="degrees of freedom used for the t distribution in flowClust; Inf for pure Gaussian"),
  make_option("--nu.est", type="integer", default=0, 
              help="A numeric indicating whether ‘nu’ is to be estimated or not; 0: no, 1: non-specific, 2: cluster-specific estimation of nu"),
  make_option("--trans", type="integer", default=1, 
              help="A numeric indicating whether the Box-Cox transformation
          parameter is estimated from the data; 0: no, 1: non-specific, 2: cluster-specific estim. of lambda"), ## TODO: try 2
  make_option("--randomStart", type="integer", default=1, 
              help="number of kmeans initializations"), 
  make_option(c("--merge"), action="store_true", default=FALSE,
              help="use flowMerge to merge best BIC clustering"),
  make_option(c("--recluster"), action="store_true", default=FALSE,
              help="use k-means to re-cluster best BIC clustering"),
  ## OUTPUT OPTIONS
  make_option(c("--jobs"), type="character", default="distribution,timeseries,fourier,rain,clustering", 
              help=",-separated list of rseults to save as csv: distributions,timeseries,fourier,clustering; default is to save all, clustering only if specified by separate option --cluster, and fourier results will contain p-values only of --perm was specified and >1"),
  make_option(c("--out.name"), type="character", default="", 
              help="file name prefix of summary file"),
  make_option(c("--out.path"), type="character", default=".", 
              help="directory path for output data (figures, csv files)"),
  make_option(c("-v", "--verb"), type="integer", default=1, 
              help="0: silent, 1: main messages, 2: warning messages"),
  make_option(c("--fig.type"), type="character", default="png",
              help="figure type, png or pdf [default %default]"),
  make_option(c("--save.rdata"), action="store_true", default=FALSE,
              help="save complete analysis as RData file (big!)"))

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## process comma-separated list arguments
lst.args <- c(dft.range="numeric",
              read.rng="numeric",
              stypes="character",
              jobs="character",
              K="numeric")
for ( i in 1:length(lst.args) ) {
    idx <- which(names(opt)==names(lst.args)[i])
    ## get individual values
    tmp <- as.list(unlist(strsplit(opt[[idx]], ",")))
    ## expand ranges
    if ( lst.args[i]=="numeric" & length(tmp)>0 )
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
                  paste(opt[[i]],collapse=","), "\n",sep=""))
    arg <- names(opt)[i]
    assign(arg, opt[[arg]])
}
if ( verb>0 )
    cat(paste("\n"))

## LOADING PACKAGES FOR REQUESTED JOBS
if ( "clustering" %in% jobs ) {
    suppressPackageStartupMessages(library("flowClust"))
    suppressPackageStartupMessages(library("flowMerge"))
}
if ( "rain" %in% jobs ) 
    suppressPackageStartupMessages(library("rain"))


### START 

## generate locally in cwd
dir.create(out.path, showWarnings = FALSE)

if ( ncpu>1 ) # load parallel for flowClust
    suppressPackageStartupMessages(library("parallel"))


### LOAD DATA

## load chromosome index
cf <- read.table(chrfile,sep="\t",header=FALSE)
chrS <- c(0,cumsum(cf[,ncol(cf)])) ## index of chr/pos = chrS[chr] + pos

## load time-series and oscillation data: coor, ts, osc
if ( verb>0 )
    cat(paste("Loading data\t",time(),"\n"))
load(datafile)

## set missing read range to all!
if ( length(read.rng)==0 | any(is.na(read.rng)) ) 
  read.rng <- 1:ncol(ts)


## oscillation parameters (for first two cycles only)
phase <- osc[,"phase" ]
pval <- osc[,"rpval"]
phase[pval >= c(min(1,pval.thresh))] <- NA # =1 rm 360° artefact in osci set
logit <- function(p) log(p/(1-p))
lpvl <- logit(pval)
lpvl[is.infinite(lpvl)] <- NA

## PHASE WEIGHTS by p-values
## TODO: use weights? seems to have little effect
wght <- rep(1,length(pval)) # phase.weight default equiv. to no weight
if ( phase.weight ) wght <-  1-pval 
##if ( phase.weight=="log" ) wght <- -log2(pval) # TODO: fix infinite weights!

## LOAD SEGMENTS

if ( verb>0 )
    cat(paste("Loading segments\t",time(),"\n"))
segs <- read.table(infile,sep="\t",header=TRUE, comment.char="",
                   stringsAsFactors=FALSE)

## add type ALL column,
## allows to pass ALL to cmdline option --stypes to avoid typesplitting
if ( stypes[1]=="ALL" & !"ALL"%in%colnames(segs) ) {
    segs <- cbind.data.frame(segs, all="all")
    stypes <- "all"
    typecol <- "all"
}

## reduce to requested segment types
if ( stypes[1]=="" )  
    stypes <- sort(unique(as.character(segs[,typecol])))
segs <- segs[as.character(segs[,typecol])%in%stypes,]



if ( fuse.segs ) {
    fuse <- segs[2:nrow(segs),"fuse"]==1
    cat(paste("NOTE: FUSING", sum(fuse), "SEGMENTS, from segment types:\n",
              paste(unique(segs[fuse,typecol]),collapse="; "),"\n"))
    fsegs <- segs[c(TRUE,!fuse),]
    
    fsegs[,"end"] <- segs[c(!fuse,TRUE),"end"]
    segs <- fsegs
}


## replace genome coordinates by continuous index
segs <- coor2index(segs,chrS)

## split by type
lst <- split(segs,segs[,typecol])
if ( length(stypes)>0 ) 
    lst <- lst[names(lst)%in%stypes]
sgtypes <- names(lst)
segnum <- unlist(lapply(lst,nrow))

## CALCULATE OSCI PARAMETERS FOR SEGMENTS

## cluster number of max BIC clustering
## will be used in segmentation analysis together with results from
## segmentLengths and testSegments
if ( "clustering" %in% jobs ) {
    clnum <- matrix(NA, length(sgtypes),ncol=4) 
    rownames(clnum) <- sgtypes
    colnames(clnum) <- c("K","BIC","NUMCL", "TOT")
}

for ( type in sgtypes ) {

    if ( !exists("type", mode="character") )
      type <- sgtypes[1] ## NOTE: DEVEL HELPER - NOT REQUIRED

    if ( verb>0 )
        cat(paste(type, "\t",time(),"\n"))
    fname <- gsub(":",".",type) # FOR FILENAMES
    
    sgs <- lst[[type]]

    ## READ-COUNT DISTRIBUTIONS OF SEGMENTS
    ## pvs required for filtering before clustering
    phs <- pvs <- rds <- NULL
    if ( any(c("distribution","clustering") %in% jobs) ) {

        if ( verb>0 )
          cat(paste("segment read statistics\t",time(),"\n"))
        ##file.name <- file.path(out.path,paste(fname,"_dynamics",sep=""))
        ##if ( file.exists(file.name) & !redo) {
        ##    cat(paste("\talready done\n")
        ##}
        
        ## add phase distribution, weighted by p-values!
        phs <- t(apply(sgs,1,function(x)
                       phaseDist(phase[x["start"]:x["end"]],
                                 w=wght[x["start"]:x["end"]])))

        ## comparing sg0002_raw_ash_icor_463 with sg0002_raw_ash_icor_464
        ## from D:dft1-7.dcash.snr_T:raw_K:12_S:icor_E:3_M:75_nui:3
        test.circ.test <- FALSE
        if (test.circ.test) {
            idx1 <- grep("sg0002_raw_ash_icor_463",sgs[,idcol],value=F)
            idx2 <- grep("sg0002_raw_ash_icor_464",sgs[,idcol],value=F)
            x1 <- circular(phase[sgs[idx1,"start"]:sgs[idx1,"end"]],
                           type="angles",units="degrees")
            x2 <- circular(phase[sgs[idx2,"start"]:sgs[idx2,"end"]],
                           type="angles",units="degrees")
            plot(x1);points(x2,col=2)
            watson.two.test(x1,x2)
        }
        
        ## raw pvalue distribution
        ## NOTE only p.signif is used, fraction of signif. oscill. reads!
        pvs <- t(apply(sgs,1,function(x)
            pvalDist(pval[x["start"]:x["end"]],pval.thresh.sig)))

        ## total range of expression
        ## NOTE: ONLY TAKING FIRST 19/20 TIMEPOINTS TO SKIP SHIFT IN THE END?
        rds <- t(apply(sgs,1,function(x)
            readDist(c(ts[x["start"]:x["end"],read.rng]))))
        
        ## write out phase, pval and read-count distributions
        ## convert back to chromosome coordinates
        sgdst <- data.frame(ID=sgs[,idcol],rds,phs,pvs)
        file.name <- file.path(out.path,paste(fname,"_dynamics",sep=""))
        write.table(sgdst,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)
    }
    
    if ( !any(c("rain","timeseries","fourier","clustering") %in% jobs) )
      next

    ## AVERAGE SEGMENT TIME-SERIES
    ## required for all following options
    
    ## NOTE: 20161231 - all smoothing attempts failed, and the
    ## simple mean of each nucleotide works best
    ## TODO: rm extrema (take only center 85%)?
    ## runmed (avgk>1) before global mean?
    ## why doesnt median work (clustering fails)?
    ## try to filter for most significant oscillators?
      
    if ( verb>0 )
      cat(paste("segment time series\t",time(),"\n"))
    sgavg <- function(x) {
        rng <- as.numeric(x["start"]):as.numeric(x["end"])
        rds <- ts[rng,]
        if ( filter.reads )
          rds <- rds[pval[rng] < pval.thresh.sig]
        ## get average
        avg <- segmentAverage(rds,
                              avg="mean",endcut=endcut,k=avgk,endrule="median",
                              mean.ratio=mean.ratio)
        
        avg
    }
    avg <- t(apply(sgs,1,sgavg))
    
    
    ## write out average timeseries
    if ( "timeseries" %in% jobs ) {
        sgts <- data.frame(ID=sgs[,idcol], avg)
        file.name <- file.path(out.path,paste(fname,"_timeseries",sep=""))
        write.table(sgts,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)
    }
    

    ## get DFT, use time series processing from segmenTier
    ## this will be written out; re-calculated after filtering
    ## below for clustering
    ## NOTE: using read.rng here as well (above for total read count)
    dft <- NULL
    ##if ( any(c("fourier","clustering") %in% jobs) ) {
    if ( any(c("fourier") %in% jobs) ) {

        if ( verb>0 )
          cat(paste("fourier transform\t",time(),"\n",sep=""))
        if ( perm>0 & verb>0 )
          cat(paste("permutations\t", perm,"\n",sep=""))
        
        tset <- processTimeseries(avg[,read.rng],na2zero=TRUE, use.fft=TRUE,
                                  smooth.time=smooth.time,
                                  trafo=trafo, perm=perm,
                                  dft.range=dft.range, dc.trafo=dc.trafo,
                                  use.snr=use.snr,low.thresh=-Inf, verb=verb)
        
        ## write out phase, pval and DFT from segment averages
        dft <- tset$dft
        if ( perm>0 ) {
            pvl <- tset$pvalues
            colnames(pvl) <- paste(colnames(pvl),"p",sep="_")
            dft <- data.frame(dft,pvl)
        }
        sgdft <- data.frame(ID=sgs[,idcol], dft)
        file.name <- file.path(out.path,paste(fname,"_fourier",sep=""))
        write.table(sgdft,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)

    }

    ## use rain
    ## TODO: establish exact period from DO, then use rain
    ## use only first 20 time-points here as well
    if ( any(c("rain") %in% jobs) ) {
        if ( verb>0 )
          cat(paste("rain osci stastistics\t",time(),"\n"))
        rn <- rain(t(avg[,read.rng]), period=period, deltat=deltat)
        sgrain <- data.frame(ID=sgs[,idcol], rn)
        file.name <- file.path(out.path,paste(fname,"_rain",sep=""))
        write.table(sgrain,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)
    }
    
    if ( !"clustering" %in% jobs ) next
    if ( verb>0 )
        cat(paste("segment clustering\t",time(),"\n"))

    ## CLUSTER DFT OF AVERAGE TIME-SERIES
    ## TODO: instead, collect DFT cluster centers of segments
    ## and use these as centers, or cluster these instead?
    ## OR: trust cluster-segment and take only those reads that
    ## were in major clusters

    ## TODO: use permutation results as filter; optionally load
    ## permutation from previous run!
    ## TODO: move back to unsig! was the most sensible!

    ## FILTERS
    ## NOTE: 20170118
    ## currently only `unsig' works well in filtering segment mass in
    ## exponent>1; however, this is based on prior pvalue of read-counts
    ## while it may be nicer to have a filter based on segment averages

    ##nosig <- dft[,"X2_p"] >= 0.05; nosig[is.na(nosig)] <- TRUE # NO SIG @PERIOD
    #### nonsig - NO SIGNIFICANT PVALUE IN ANY FOURIER COMPONENT!
    ##nonsig <- !apply(tset$pvalues,1,function(x) any(x< .05))
    ##nonsig[is.na(nonsig)] <- TRUE
    ##lowex <- rds[,"r.0"]>.9        # MINIMAL FRACTION OF READ COUNTS>0
    ##len <- sgs[,"end"]-sgs[,"start"]+1
    ##short <- len < 100             # LONGER THEN 150
    
    ## use prior RAIN calculation as filter
    unsigr <- rep(FALSE, nrow(avg))
    if ( with.rain!="" ) {
        sgrain <- read.delim(file.path(with.rain,paste0(fname,"_rain.csv")),
                             row.names=1,stringsAsFactors=FALSE)[as.character(sgs[,idcol]),]

        ## FILTER
        unsigr <- sgrain[,"pVal"] >= pval.thresh.sig

        ## plot
        ## TODO: this is used in paper; make fit for supplementary material
        file.name <- file.path(out.path,paste0(fname,"_filter_rain"))
        plotdev(file.name,width=4,height=4,type=fig.type,res=300)
        rpcdf <- ecdf(sgrain[,"pVal"])
        plot(rpcdf,xlab="rain p-value")
        points(pval.thresh.sig,rpcdf(pval.thresh.sig))
        points(.995,rpcdf(.995))
        legend("top",legend=c(paste0(sum(!unsigr), " oscillate")))
        dev.off()
        #plot(sgrain[,"pVal"],pvs[,1]) # NOTE slight correlation!

    }
    ## use prior permutation calculation as filter
    unsigp <- rep(FALSE, nrow(avg))
    if ( with.permutation!="" ) {
        sgdft <- read.delim(file.path(with.permutation,
                                      paste0(fname,"_fourier.csv")),
                            row.names=1,stringsAsFactors=FALSE)[as.character(sgs[,idcol]),]
        ## FILTER
        sgdft[is.na(sgdft[,"X2_p"]),"X2_p"] <- 1
        unsigp <- sgdft[,"X2_p"] >= pval.thresh.sig

        ## plot
        file.name <- file.path(out.path,paste0(fname,"_filter_permutation"))
        plotdev(file.name,width=4,height=4,type=fig.type,res=300)
        ppcdf <- ecdf(sgdft[,"X2_p"]) ## TODO: this must be argument
        plot(ppcdf,xlab="permutation p-value")
        points(pval.thresh.sig,ppcdf(pval.thresh.sig))
        legend("top",legend=c(paste0(sum(!unsigp), " oscillate")))
        dev.off()
        
    }

    ## FILTER: minimal expressed time-points per segment
    mintpt <- 12
    npoints <- rowSums(avg>0)
    fewpoints <- npoints <= mintpt

    ## plot fewpoints
    file.name <- file.path(out.path,paste0(fname,"_filter_numpoints"))
    plotdev(file.name,width=4,height=4,type=fig.type,res=300)
    npcdf <- ecdf(npoints)
    plot(npcdf,xlab="number of expressed timepoints")
    points(mintpt,npcdf(mintpt),cex=2)
    legend("top",legend=c(paste0(sum(!fewpoints), " expressed")))
    dev.off()

    ## FILTER: no single significant read-count (< pval.thresh.sig)
    unsig <- pvs[,"p.signif"] == 0 
    pcdf <- ecdf(pvs[,"p.signif"])
    file.name <- file.path(out.path,paste0(fname,"_filter_signifraction"))
    plotdev(file.name,width=4,height=4,type=fig.type,res=300)
    plot(pcdf,xlab="fraction of significant reads")
    points(.01,pcdf(.01))
    dev.off()
    
    ## FILTER: total expresssion vs. rain
    minexp <- .05
    tot <-ash(rds[,"r.mean"])
    lowex <- tot<minexp

    ## plot total
    file.name <- file.path(out.path,paste0(fname,"_filter_total"))
    plotdev(file.name,width=4,height=4,type=fig.type,res=300)
    tcdf <- ecdf(tot)
    plot(tcdf,xlab="mean read-count")
    points(minexp,tcdf(minexp))
    legend("right",legend=c(paste0(sum(!lowex), " expressed")))
    dev.off()
    
    ## SHORT
    minlen <- 90
    len <- sgs[,"end"]-sgs[,"start"]+1
    short <- len < minlen             # LONGER THEN 150

    ## plot lengths
    file.name <- file.path(out.path,paste0(fname,"_filter_length"))
    plotdev(file.name,width=4,height=4,type=fig.type,res=300)
    lcdf <- ecdf(len)
    plot(lcdf,xlab="segment length")
    points(minlen,lcdf(minlen))
    legend("right",legend=c(paste0(sum(!short), " long")))
    dev.off()

    ## filter combination
    noise <- lowex | short | fewpoints

    ## SELECT FILTER
    filters <- cbind(lowex=lowex, fewpoints=fewpoints, short=short, 
                     unsig=unsig, unsig.r=unsigr, unsig.p=unsigp,
                     noise=noise)
    
    rmvals <- filters[,cl.filter]
    dat <- avg 
    dat[rmvals,] <- 0 # set to zero, will be removed in processTimeseries
    
    ## TODO: lowly expressed AND non-signficant oscillation LOOKS good!!
    ##table(lowex,unsigr)
    ##rmvals <- lowex&unsigr
    
    if ( verb>0 ) {
        cat(paste("clustered segments\t",sum(!rmvals),"\n"))
        cat(paste("noise segments\t",sum(rmvals),"\n"))
    }

    tset <- processTimeseries(dat,na2zero=TRUE,use.fft=TRUE,
                              smooth.time=smooth.time, trafo=trafo,
                              perm=0, dft.range=dft.range, dc.trafo=dc.trafo,
                              use.snr=use.snr,low.thresh=-Inf)


    ## cluster by flowClust
    testnew <- TRUE ## ALLOWS MULTIPLE EQUAL K!
    if ( testnew ) {
        fcset <- clusterTimeseries2(tset, K=K, method="flowClust",
                                    parameters=c(B=B,tol=tol,lambda=lambda,
                                                 nu=nu,nu.est=nu.est,
                                                 trans=trans, randomStart=0))
    } else {
        fcset <- flowclusterTimeseries(tset, ncpu=ncpu, K=K, selected=fixedK,
                                       B=B, tol=tol, lambda=lambda, merge=merge,
                                       nu=nu, nu.est=nu.est, trans=trans)
    }
    
    ## save all as RData
    ## temporary; until below is fixed
    #if ( save.rdata ) 
    #  save(sgs, rds, phs, pvs, dft, tset, fcset, #sgcls,
    #       file=file.path(out.path,paste0(fname,".RData",sep="")))

    ## get BIC and best clustering
    selected <- selected(fcset) # cluster number of max BIC 
    cls <- fcset$clusters[,selected] # clusters at max BIC
    # Error in fcset$clusters[, selected] : subscript out of bounds

    bic <- fcset$bic  # BIC
    icl <- fcset$icl  # ICL
    max.cli <- fcset$max.cli  # cluster number at max ICL
    max.clb <- fcset$max.clb  # cluster number at max BIC
    
    max.bic <- max(bic,na.rm=TRUE) # max BIC
    max.icl <- max(icl, na.rm=T)   # max ICL

    ## store cluster number, max BIC, numbers of clustered and total segments
    usedk <- fcset$usedk[selected] # NOTE: 
    clnum[type,] <- c(K=usedk, BIC=max.bic,
                      NUMCL=sum(!tset$rm.vals), TOT=nrow(avg))

    ## and add selected global clustering to segment table
    ## write out clusters
    sgcls <- data.frame(ID=sgs[,idcol],sgCL=cls)


    ## flowMerge 
    mselected <- NULL
    if ( merge ) {
        fcset <- mergeCluster(tset, fcset, selected=selected(fcset))
        mselected <- fcset$merged # cluster number of merged clustering
        if ( !is.null(mselected) ) {
            mcls <- fcset$clusters[,mselected] # clusters of merged clustering
            mrg.cl <- fcset$merged.K # cluster number of merged clustering
            ## add to data frame
            sgcls <- cbind.data.frame(sgcls,mCL=mcls)
        }
    }

    ## re-cluster with kmeans
    if ( recluster ) {
        fcset <- reCluster(tset, fcset, select=FALSE)
        rselected <- fcset$reclustered
        sgcls <- cbind.data.frame(sgcls, rCL=fcset$clusters[,rselected])
    }

    file.name <- file.path(out.path,paste(fname,"_clusters",sep=""))
    write.table(sgcls,file=paste(file.name,".csv",sep=""),quote=FALSE,
                sep="\t",col.names=TRUE,row.names=FALSE)


    ## RESORT CLUSTERING
    ## TODO: sort by phase and re-cluster; use vector phs
    ## of average segment phases!
    ## cls.phase <- sapply(cls.srt, function(x)
    ##   phaseDist(phase[cls==x],w=1-rp[cls==x,"pVal"]))
    ## cset$sorting[[bestKcol]] <- cls.srt
    ## ## re-color
    ## cset <- colorClusters(cset)

    
    ## save all as RData
    if ( save.rdata ) 
      save(sgs, rds, phs, pvs, dft, tset, fcset, sgcls, fname,
           file=file.path(out.path,paste0(fname,".RData",sep="")))
    
    ## PLOT CLUSTERING

    if ( verb>0 )
        cat(paste("\tplotting time series clustering\t",time(),"\n"))

    ## re-do tset without removing unsignificant!
    ## "plot-tset"
    pset <- processTimeseries(avg,na2zero=TRUE,use.fft=TRUE,
                              smooth.time=smooth.time, trafo=trafo,
                              perm=0, dft.range=dft.range, dc.trafo=dc.trafo,
                              use.snr=use.snr,low.thresh=-Inf)

    ## plot BIC
    file.name <- file.path(out.path,paste(fname,"_BIC",sep=""))
    plotdev(file.name,width=4,height=4,type=fig.type,res=300)
    par(mai=c(.7,.7,0.1,0.1),mgp=c(1.5,.5,0),tcl=-.3)
    plotBIC(fcset, norm=TRUE)
    dev.off()
     
    ## plot DFT
    file.name <- file.path(out.path,paste(fname,"_DFT",sep=""))
    plotdev(file.name,type=fig.type,res=300,
            width=round(length(dft.range)),height=4)
    par(mfcol=c(2,round(length(dft.range)/2)),
        mai=c(.5,.5,0.1,0),mgp=c(1.5,.5,0),tcl=-.3)
    plotDFT(tset, fcset, cycles=dft.range, pch=1, cex=.5)
    dev.off()

    ## plot best BIC
    file.name <- file.path(out.path,paste0(fname,"_osc_",selected))
    plotdev(file.name,width=4,height=9,type=fig.type,res=300)
    plotClusters(pset,fcset,k=selected,norm="meanzero")
    dev.off()
    ## plot merged
    if ( merge & !is.null(mselected) ) {
        file.name <- file.path(out.path,paste0(fname, "_osc_",mselected))
        plotdev(file.name,width=4,height=9,type=fig.type,res=300)
        plotClusters(pset,fcset,k=mselected,norm="meanzero")
        dev.off()
        ## plot DFT
        file.name <- file.path(out.path,paste0(fname,"_DFT_",mselected))
        plotdev(file.name,type=fig.type,res=300,
                width=round(length(dft.range)),height=4)
        par(mfcol=c(2,round(length(dft.range)/2)),
            mai=c(.5,.5,0.1,0),mgp=c(1.5,.5,0),tcl=-.3)
        tmp<-fcset;tmp$selected <- mselected
        plotDFT(tset, tmp, cycles=dft.range, pch=1, cex=.5)
        dev.off()
    }
    if ( recluster ) {
        file.name <- file.path(out.path,paste(fname, "_osc_",rselected,sep=""))
        plotdev(file.name,width=4,height=9,type=fig.type,res=300)
        plotClusters(pset,fcset,k=rselected,norm="meanzero")
        dev.off()
        ## plot DFT
        file.name <- file.path(out.path,paste0(fname,"_DFT_",rselected))
        plotdev(file.name,type=fig.type,res=300,
                width=round(length(dft.range)),height=4)
        par(mfcol=c(2,round(length(dft.range)/2)),
            mai=c(.5,.5,0.1,0),mgp=c(1.5,.5,0),tcl=-.3)
        tmp<-fcset;tmp$selected <- rselected
        plotDFT(tset, tmp, cycles=dft.range, pch=1, cex=.5)
        dev.off()
    }
 }

## write out summary for analysis over segmentation types!
if ( "clustering" %in% jobs ) {
    if ( verb>0 )
        cat(paste("saving clustering results\t",time(),"\n"))
    clnum <- cbind(ID=rownames(clnum), clnum)
    if ( out.name == "" ) {
        file.name <- "clustering"
    } else {
        file.name <- paste(out.name,"clustering",sep="_")
    }
    file.name <- file.path(out.path,file.name)
    write.table(clnum,file=paste(file.name,".csv",sep=""),quote=FALSE,
                sep="\t",col.names=TRUE,row.names=FALSE)
}

if ( verb>0 )
  cat(paste("DONE\t",time(),"\n"))
if ( !interactive() ) quit(save="no")

## TODO (but elsewhere) :
## PLOT PHASE DIST etc.
## plot phase dist (for diff. weights?),
## read-dist and pval-dist


### TODO - CLASSIFY OVERLAPS
### TODO - CHROMOSOME PLOTS
## get overlaps with cltranscripts and clgenes as in testSegments.R
## classify into upstream - covers - downstream of cluster genes
## calculate enrichment (Jaccard?) per cluster
## potentially: use super-clusters
## plot next to time-course

### MANUAL
##
load("all_20160525/K16_k1_icor3_filt_osc_K14m7.RData")

## segments with 4 peaks!
## mostly D/cd, either full ORF or downstream, sometimes upstream
soi <- seg[which(seg[,"mCL"] == 5),]
