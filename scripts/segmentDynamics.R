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


suppressPackageStartupMessages(library("stringr")) # for 0-padded filenames
suppressPackageStartupMessages(library("flowClust"))
suppressPackageStartupMessages(library("flowMerge"))
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
  make_option(c("--stypes"), type="character", default="", 
              help="sub-set of segments in column 'type'"),
  ## OSCILLATION SETTINGS
  make_option("--pval.thresh", default=1,
              help="phases above this thresh will be set to NA [default %default]"),
  make_option("--phase.weight", action="store_true", default=FALSE,
              help="use weight for calculating average phases of segments: 1-p.value [default %default]"),
  make_option("--pval.thresh.sig", default=1,
              help="pvals above will be counted as non-significant [default %default]"),
  make_option("--read.rng", type="character", default="",
              help="range of time-points for total read-count and Fourier calculations (but not used for clustering!), comma-separated list of integers"),
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
  ##    make_option("--smooth.time.plot", type="integer", default=3, # so far best! 
  ##                help="as smooth.time but only for plotting clusters not used of analysis"),
  ## FLOWCLUST PARAMETERS
  make_option("--ncpu", type="integer", default=1, 
              help="number of available cores for flowClust"),
  ##    make_option("--seed", type="integer", default=1, 
  ##                help="seed for the random number generator before calling flowclust to get stable clustering results (TODO: set.seed, does it work for flowclust?)"),
  make_option(c("--K"), type="character", default="12:20", 
              help="number of clusters to use in flowClust, comma-separated list of integers and colon-separated ranges [default %default]"),
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
  ## OUTPUT OPTIONS
  make_option(c("--jobs"), type="character", default="distribution,timeseries,fourier,clustering", 
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
    if ( lst.args[i]=="numeric" )
        for ( j in 1:length(tmp) ) { # only for numeric modes
            tmp2 <- unlist(strsplit(tmp[[j]], ":"))
            if ( length(tmp2)>1 ) {
                tmp2 <- as.numeric(tmp2)
                tmp[[j]] <- tmp2[1]:tmp2[2]
            }
        }
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



### START 

## generate locally in cwd
dir.create(out.path, showWarnings = FALSE)

if ( ncpu>1 ) # load parallel for flowClust
    suppressPackageStartupMessages(library("parallel"))


### LOAD DATA

## load chromosome index
cf <- read.table(chrfile,sep="\t",header=FALSE)
chrS <- c(0,cumsum(cf[,3])) ## index of chr/pos = chrS[chr] + pos

## load time-series and oscillation data: coor, ts, osc
if ( verb>0 )
    cat(paste("Loading data:",datafile,"\t",time(),"\n"))
load(datafile)

## set missing read range to all!
if ( length(read.rng)==0 ) 
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
    cat(paste("Loading segments:",infile,"\t",time(),"\n"))
segs <- read.table(infile,sep="\t",header=TRUE)

## reduce to requested segment types
if ( length(stypes)>0 ) 
    segs <- segs[as.character(segs[,"type"])%in%stypes,]


if ( fuse.segs ) {
    fuse <- segs[2:nrow(segs),"fuse"]==1
    cat(paste("NOTE: FUSING", sum(fuse), "SEGMENTS, from segment types:\n",
              paste(unique(segs[fuse,"type"]),collapse="; "),"\n"))
    fsegs <- segs[c(TRUE,!fuse),]
    
    fsegs[,"end"] <- segs[c(!fuse,TRUE),"end"]
    segs <- fsegs
}


## replace genome coordinates by continuous index
segs <- coor2index(segs,chrS)

## split by type
lst <- split(segs,segs$type)
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
    phs <- pvs <- rds <- NULL
    if ( "distribution" %in% jobs ) {

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
        ## raw pvalue distribution
        ## NOTE only p.signif is used, fraction of signif. oscill. reads!
        pvs<-t(apply(sgs,1,function(x)
                     pvalDist(pval[x["start"]:x["end"]],pval.thresh.sig)))

        ## total range of expression
        ## NOTE: ONLY TAKING FIRST 20 TIMEPOINTS TO SKIP SHIFT IN THE END?
        rds <-t(apply(sgs,1,function(x)
                      readDist(c(ts[x["start"]:x["end"],read.rng]))))
        
        ## write out phase, pval and read-count distributions
        ## convert back to chromosome coordinates
        sgdst <- data.frame(ID=sgs[,"ID"],rds,phs,pvs)
        file.name <- file.path(out.path,paste(fname,"_dynamics",sep=""))
        write.table(sgdst,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)
    }
    
    if ( !any(c("timeseries","fourier","clustering") %in% jobs) )
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
      cat(paste("segment average time series\t",time(),"\n"))
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
        sgts <- data.frame(ID=sgs[,"ID"], avg)
        file.name <- file.path(out.path,paste(fname,"_timeseries",sep=""))
        write.table(sgts,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)
    }
    

    ## get DFT, use time series processing from segmenTier
    ## this will be written out; re-calculated after filtering
    ## below for clustering
    ## NOTE: using read.rng here as well (above for total read count)
    dft <- NULL
    if ( any(c("fourier","clustering") %in% jobs) ) {

        if ( perm==0 & verb>0 )
          cat(paste("discrete fourier transform\t",time(),"\n",sep=""))
        if ( perm>0 & verb>0 )
          cat(paste("fourier p-values (",perm,
                    ") permutations\t", time(),"\n",sep=""))
        
        tset <- processTimeseries(avg[,read.rng],smooth.time=smooth.time,
                                  trafo=trafo, perm=perm,
                                  dft.range=dft.range, dc.trafo=dc.trafo,
                                  use.snr=TRUE,low.thresh=-Inf, verb=verb)
        
        ## write out phase, pval and DFT from segment averages
        dft <- tset$dft
        if ( perm>0 ) {
            pvl <- tset$pvalues
            colnames(pvl) <- paste(colnames(pvl),"p",sep="_")
            dft <- data.frame(dft,pvl)
        }
        sgdft <- data.frame(ID=sgs[,"ID"], dft)
        file.name <- file.path(out.path,paste(fname,"_fourier",sep=""))
        write.table(sgdft,file=paste(file.name,".csv",sep=""),quote=FALSE,
                    sep="\t",col.names=TRUE,row.names=FALSE)

        ## our use rain
        ## TODO: establish exact period from DO, then use rain
        #res <- rain(t(avg[,read.rng]),period=0.65,deltat=4/60)
    }

    
    if ( !"clustering" %in% jobs ) next
    if ( verb>0 )
        cat(paste("clustering segments\t",time(),"\n"))

    ## CLUSTER DFT OF AVERAGE TIME-SERIES
    ## TODO: instead, collect DFT cluster centers of segments
    ## and use these as centers, or cluster these instead?
    ## OR: trust cluster-segment and take only those reads that
    ## were in major clusters

    ## TODO: use permutation results as filter; optionally load
    ## permutation from previous run!
    dat <- avg 
    ##nosig <- dft[,"X2_p"] >0.1
    ## TODO: use nonsig - no significant pvalue in any fourier component!
    nonsig <- !apply(tset$pvalues,1,function(x) any(x<.05))
    #nonsig[is.na(nonsig)] <- TRUE
    #unsig <- pvs[,"p.signif"] == 0 # NO SINGLE SIGNIFICANT OSCILLATOR
    #lowex <- rds[,"t.0"]>.9        # MINIMAL FRACTION OF READ COUNTS>0
    #len <- sgs[,"end"]-sgs[,"start"]+1    
    #short <- len < 100             # LONGER THEN 150
    rmvals <- nonsig #|unsig #|short #|     # TODO: does short filter help?

    dat[rmvals,] <- 0 # set to zero, will be removed in processTimeseries
    tset <- processTimeseries(dat,smooth.time=smooth.time, trafo=trafo,
                              perm=0, dft.range=dft.range, dc.trafo=dc.trafo,
                              use.snr=TRUE,low.thresh=-Inf)


    ## cluster by flowClust
    fcset <- flowclusterTimeseries(tset, ncpu=ncpu, K=K,
                                   B=B, tol=tol, lambda=lambda,
                                   nu=nu, nu.est=nu.est, trans=trans)
 
    mselected <- fcset$merged # cluster number of merged clustering
    selected <- as.character(fcset$max.clb) # cluster number of max BIC 
    mcls <- fcset$clusters[,mselected] # clusters of merged clustering
    cls <- fcset$clusters[,selected] # clusters at max BIC

    bic <- fcset$bic  # BIC
    icl <- fcset$icl  # ICL
    max.cli <- fcset$max.cli  # cluster number at max ICL
    max.clb <- fcset$max.clb  # cluster number at max BIC
    mrg.cl <- fcset$merged.K # cluster number of merged clustering
    
    max.bic <- max(bic,na.rm=TRUE) # max BIC
    max.icl <- max(icl, na.rm=T)   # max ICL

    ## store cluster number, max BIC, numbers of clustered and total segments
    clnum[type,] <- c(K=selected, BIC=max.bic,
                      NUMCL=sum(!tset$rm.vals), TOT=nrow(avg))

    ## and add selected global clustering to segment table
    ## write out clusters
    sgcls <- data.frame(ID=sgs[,"ID"],sgCL=cls,mCL=mcls)
    file.name <- file.path(out.path,paste(fname,"_clusters",sep=""))
    write.table(sgcls,file=paste(file.name,".csv",sep=""),quote=FALSE,
                sep="\t",col.names=TRUE,row.names=FALSE)

    ## save all as RData
    if ( save.rdata ) 
      save(sgs, rds, phs, pvs, dft, tset, fcset, sgcls,
           file=paste(fname,".RData",sep=""))
    
    ## PLOT CLUSTERING

    if ( verb>0 )
        cat(paste("\tplotting time series clustering\t",time(),"\n"))

    ##if ( smooth.time!=smooth.time.plot )
    ##    tset <- processTimeseries(dat,smooth.time=smooth.time.plot, trafo=trafo,
    ##                              dft.range=dft.range, dc.trafo=dc.trafo,
    ##                              use.snr=TRUE,low.thresh=-Inf)
    
    ## plot time-courses
    ## if no smoothing was done, smooth for plots
    ndat <- tset$ts # avg # 
    navg <- log2(ndat/apply(ndat,1,mean,na.rm=T))

    ## plot BIC
    file.name <- file.path(out.path,paste(fname,"_BIC",sep=""))
    plotdev(file.name,width=4,height=4,type=fig.type,res=300)
    par(mai=c(.7,.7,0.1,0.1),mgp=c(1.5,.5,0),tcl=-.3)
    plot(K, bic, ylim=range(c(bic,icl),na.rm=T),xlab="K",ylab="BIC/ICL")
    lines(K,bic)
    points(K[is.na(bic)], rep(min(bic,na.rm=T),sum(is.na(bic))),
           pch=4, col=2)
    points(max.clb, max.bic,lty=2,pch=4,cex=1.5)
    points(K, icl, col=4,cex=.5)
    points(K[which(icl==max.icl)], max.icl,lty=2,pch=4,cex=1.5,col=4)
    arrows(x0=as.numeric(max.clb),y0=bic[as.character(max.clb)],
           x1=mrg.cl,y1=bic[as.character(max.clb)])
    abline(v=mrg.cl,lty=2)
    legend("right",pch=c(1,1,4,4),lty=c(1,NA,NA,NA),col=c(1,4,2,1),
           legend=c("BIC","ICL","failed","max. BIC"),pt.cex=c(1,.5,1,1.5))
    dev.off()
     

    ## plot merged
    file.name <- file.path(out.path,paste(fname, "_osc_K",mselected,sep=""))
    plotdev(file.name,width=4,height=9,type=fig.type,res=300)
    par(mfcol=c(length(unique(mcls)) - 1,1),mai=c(0,.5,0,.1),cex.lab=1.5,
        mgp=c(1.5,.75,0))
    for ( i in sort(unique(mcls)) ) {
        if ( i== 0 ) next
        matplot(t(navg[mcls==i,]),type="l",lty=1,lwd=.5,col=i,ylim=c(-3,2),
                ylab=sum(mcls==i))
        lines(apply(navg[mcls==i,],2,median,na.rm=TRUE),col="white",lwd=3)
    }
    dev.off()
    ## plot best BIC
    file.name <- file.path(out.path,paste(fname,"_osc_K",selected,sep=""))
    plotdev(file.name,width=4,height=9,type=fig.type,res=300)
    par(mfcol=c(length(unique(cls)) - 1,1),mai=c(0,.5,0,.1),cex.lab=1.5,
        mgp=c(1.5,.75,0))
    for ( i in sort(unique(cls)) ) {
        if ( i== 0 ) next
        matplot(t(navg[cls==i,]),type="l",lty=1,lwd=.5,col=i,ylim=c(-3,2),
                ylab=sum(cls==i))
        lines(apply(navg[cls==i,],2,median,na.rm=TRUE),col="white",lwd=3)
    }
    dev.off()
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
