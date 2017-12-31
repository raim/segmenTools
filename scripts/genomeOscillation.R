#!/usr/bin/Rscript

## reads a genomeData file (chr\tcoor\t\value[1]\tvalue[2]\t...\tvalue[n])
## with time series 1..n for each coordinate,
## and calculates oscillation properties of each time series.

verb <- FALSE 
header <- TRUE ## add head to csv file? useful for piece-wise parallel mode
               ## with later concatenation of results
savecsv <- FALSE ## save DFT/p-values as table files?
savefft <- FALSE ## save DFT for later call?
loadfft <- FALSE ## load previous DFT?
normts <- FALSE ## normalize time-series before FOURIER?
logts <- FALSE ## log2 of time-series before FOURIER?
gauss <- FALSE ## flowClust: gaussian?

args  <- R.utils::commandArgs(excludeReserved=TRUE, asValues=TRUE)
for ( arg in names(args)[2:length(args)] ) {
  ## fix(?) since R 3.0.0 ?
  if ( is.null(args[[arg]]) ) args[arg] <- TRUE
  assign(arg, as.character(args[arg]))
}
## input file
if ( !exists("i",mode="character") ) {
  cat(paste("no data provided, use --i=<file>\n"))
  if ( !interactive() ) quit(save="no")
}
## coordinate columns in input files, 1:length(cc) are taken
## and the value of cc is used as column name
if ( !exists("cc",mode="character") ) cc <- c("chr","coor")
## output file base name (two files are written)
if ( !exists("o",mode="character") ) o <- sub(".RData","",sub(".csv", "", i) )
## start and end values of time-series to use, 0 is default and uses all
if ( !exists("start",mode="character") ) start <- 0
if ( !exists("end",mode="character") ) end <- 0
## DFT normalization type: snr
if ( !exists("norm",mode="character") ) norm <- ""
## number of permutations
if ( !exists("perm",mode="character") ) perm <- 0

## header line? useful for split files in parallelization
if ( exists("nohead",mode="character") ) header <- FALSE
## header suffix for forward strand? used to split values into strands
if ( !exists("frw",mode="character") ) frw <- ""

## cluster? using flowClust, give range of cluster numbers
if ( !exists("cls",mode="character") ) {
  cls <- 0
} else {
  cls <- as.numeric(unlist(strsplit(cls,"\\.\\.")))
}
## chromosome range to to cluster; format: chr,start..end, e.g 1,1..Inf
if ( !exists("rng",mode="character") ) {
  rng <- ""
} else {
  tmp <- unlist(strsplit(rng,","))
  chr <- 1 ## if chrosome is not specified use pos as index!
  if ( length(tmp)==2 ) {
    chr <- as.numeric(tmp[1]) # as.numeric(unlist(strsplit(tmp[2],"\\.\\.")))
    tmp <- tmp[2]
  } 
  pos <- as.numeric(unlist(strsplit(tmp[1],"\\.\\.")))
}
## DFT range to cluster; format: 2,3,6
if ( !exists("dft",mode="character") ) {
  dft <- Inf ## use ALL components
} else {
  dft <- as.numeric(unlist(strsplit(dft,",")))
}

## flowClust: max. num. of EM iterations
if ( !exists("B",mode="character") ) B <- 500 

## parallel? use several cores in flowClust
if ( !exists("ncpu",mode="character") ) ncpu <- "1"

verb <- as.logical(verb)
savecsv <- as.logical(savecsv)
savefft <- as.logical(savefft)
loadfft <- as.logical(loadfft)
normts <- as.logical(normts)
logts <- as.logical(logts)
gauss <- as.logical(gauss)
start <- as.numeric(start)
end <- as.numeric(end)
perm <- as.numeric(perm)
B <- as.numeric(B)

#cls <- as.numeric(cls) # done above, when splitting
ncpu <- as.numeric(ncpu)

get.mesa <- function(x, start=1,frequency=0.25) {

  mesa <- t(apply(x, 1, function(x) {
    x <- ts(x,start=start, frequency=frequency);
    x <- spec.ar(x,method="yule-walker",plot=FALSE); x$spec}))
  mesa
}

get.fft <- function(x, norm="") {
  n <- floor(ncol(x)/2) +1 ## Nyquist-freq
  fft <- t(mvfft(t(x)))[,1:n]
  colnames(fft) <- c("DC",as.character(1:(n-1)))
  if ( norm == "snr" ) fft <- do.norm(fft) ## TODO: implement!
  fft
}

## relative amplitudes ('SNR')
## TODO: implement for complex numbers??
do.norm <- function(am) {
  snr <- am
  for ( i in 2:ncol(am) ) 
    snr[,i] <- am[,i]/rowMeans(abs(am[,-c(1,i)]))
  snr
}
## scale by DC component
do.scale <- function(am) {
  scl <- am[,2:ncol(am)]/am[,1]
}

## fourier permutation
do.perm <- function(x, fft=NULL, perm, norm="", verb=FALSE) {
  #require("permax")
  N <- ncol(x)
  if ( is.null(fft) ) fft <- get.fft(x,norm)
  xam <- abs(fft)/N
  pvl <- matrix(0,nrow=nrow(fft), ncol=ncol(fft))
  dimnames(pvl) <- dimnames(fft)
  ## TODO: use apply and parallel!
  for ( i in 1:perm ) {
    if ( verb & i%%round(perm/10)==0 ) cat(paste(round(i/perm,2)*100,"%, "))
    ## randomize columns and get fourier
    rft <- get.fft(x[,sample(1:ncol(x))],norm)
    ram <- abs(rft)/N
    pvl <- pvl + as.numeric(ram >= xam)
  }
  if ( verb ) cat("\n")
  pvl/perm
}

## OUT-FILE NAMES
## set up file name for DFT of timeseries_dft
dft.name <- o
if ( normts ) dft.name <- paste(dft.name, "_normts",sep="")
if ( logts )  dft.name <- paste(dft.name, "_log2",sep="")
dft.name <- paste(dft.name,"_dft",sep="")
if (norm!="") dft.name <- paste(dft.name, "_", norm,sep="")

## set up file name for clustering files
if ( max(cls)>1 ) {
  ## convert chromosome range for cluster outfile name
  rng.name <- "" 
  if ( rng!="" )
    rng.name <- paste("chr",sub(",","_",sub("\\.\\.","_",rng)),"_",sep="")
  ## convert cluster K range for cluster outfile name
  clsrng <- paste(unique(range(cls)),collapse="-")
  ## was the DFT normalized?
  cl.name <- paste("DFT_",
                   ifelse(norm=="","",paste(norm,"_",sep="")),
                   ifelse(is.infinite(dft[1]), "all",
                          paste(range(dft),collapse="_")),
                   "_K",clsrng,sep="")
  ## construct cluster outfile name
  ## was the time-series normalized?
  ## is flowClust used with Gaussian only?
  cl.name <- paste(rng.name,
                   ifelse(normts, "normts_",""),
                   ifelse(gauss, "gauss_",""),
                   cl.name,sep="")
}


## TODO: use MESA (wrapper for fortran) to obtain amplitude and p-values
## and cut or pad data to get phase via DFT
## TODO - DOMAINS:
## define domains of consistent phase (filtered by amplitude), e.g.,
## use correlation as von-mises distribution
## TODO - TRANSCRIPTS: define transcripts (cufflinks) via amplitudes (filtered
## by consistent phase)?

if ( loadfft ) {
  ## save settings
  settings <- tempfile()
  save.image(settings)
  ## LOAD PREVIOUS FFT: contains only fdat, coor, zs
  outfile <- paste(dft.name,".RData",sep="")
  if ( verb )
    cat(paste("loading DFT:     \t", date(), "\nfile:\t", outfile), "\n",
        file=stdout())
  load(outfile)
  load(settings) # 20160425 - not required since we saved only data?
  unlink(settings)
} else {
  if ( verb )
    cat(paste("reading genome data:     \t", date(), "\n"),file=stdout())
  
  ## load complete data file
  if ( length(grep(".RData$", i)) > 0 )  {
    load(i)
  } else {
    ## TODO: use skip and nrow
    ## use chrS indexing to get line numbers
    ## dat <- read.table(i,sep="\t",header=header,skip=segment[1],nrows=segment[2]-segment[1])

    dat <- read.delim(i,header=header)
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
  }
  
  ## split off coordinates
  coor <- dat[,1:length(cc),drop=FALSE]
  ## add header (for split files in parallel use)
  if ( ! header ) colnames(coor) <- cc
  dat <- dat[,- c(1:length(cc))]
  
  ## filter by specified start:end values
  if ( start == 0 ) start <- 1
  if ( end ==0 ) end <- ncol(dat)
  if ( start != 1 | end != ncol(dat) ) {
    dat <- dat[,start:end]
    o <- paste(o, "_", start,"-",end, sep="")
  }
  
  ## split columns by strand filter (forward/reverse strands)
  if ( frw!="" ) {
    if ( !header ) {
      cat("ERROR: no header but forward strand filter provided!\n")
      quit(save="no")
    }
    nam <- colnames(dat)
    dat <- rbind(dat[,grep(frw,nam)],dat[,grep(frw,nam,invert=T)])
    colnames(dat) <- sub(frw,"",grep(frw,nam,value=T))
    coor <- rbind(cbind(coor,strand=rep(1,nrow(coor))),
                  cbind(coor,strand=rep(-1,nrow(coor))))
    cat("Note: results contain values for both strands!\n")
  }
  
  ## convert to numeric
  nam <- colnames(dat)
  dat <- matrix(as.numeric(dat), ncol=ncol(dat), nrow=nrow(dat))
  colnames(dat) <- nam

  ## TODO: convert to function from here
  
  ## filter zero nucleotides
  zs <- apply(dat,1,sum)==0
  dat <- dat[!zs,,drop=FALSE]

  ## normalize time series?
  if ( normts ) 
    dat <- dat/rowMeans(dat)
  if ( logts ) ## OBSOLETE, NOT USED, will cause -Inf for 0s
    dat <- log2(dat)
  
  if ( verb ) 
    cat(paste("time series normalization:\t",as.character(normts),"\n"),
        file=stdout())
  if ( verb ) 
    cat(paste("time series log2:\t",ifelse(logts,"yes","no"),"\n"),
        file=stdout())
  

  ## result matrix
  N <- floor(ncol(dat)/2) +1
  fdat <- matrix(NA, nrow=nrow(coor), ncol=N)
  colnames(fdat) <- c("DC",as.character(1:(N-1)))
  if ( perm>0 ) pdat <- fdat

  if ( nrow(dat)>0 ) {
    
    ## FOURIER TRANSFORM 
    if ( verb ) {
      cat(paste("discrete fourier transform:\t", date(),"\n"),file=stdout())
      cat(paste("DFT normalization:\t",ifelse(norm=="","none",norm),"\n"),
          file=stdout())
    }
    fft <- get.fft(dat, norm=norm)
    fdat[!zs,] <- fft
    
    ## P-VALUES - takes long
    if ( perm>0 ) {
      if ( verb )
        cat(paste("starting",perm, "permutations:\t",date(),"\n"),file=stdout())
      pvl <- do.perm(dat, fft, perm=perm, norm=norm, verb=verb)
      if ( verb )
        cat(paste("finished",perm,"permutations:\t",date(),"\n"),file=stdout())
      pdat[!zs,] <- pvl
    }  
    if ( verb )
      cat(paste("finished fourier analysis:   \t",date(),"\n"),file=stdout())
    
  }else{
    if ( verb )
      cat(paste("no values, finished at:      \t",date(),"\n"),file=stdout())
  }
 
  ## save fft RData to re-load for clusterings
  if ( savefft ) {
    outfile <- paste(dft.name,".RData",sep="")
    if ( verb )
      cat(paste("saving fourier analysis:   \t",outfile,"\n"),file=stdout())
    #save.image(outfile)
    save('fdat','coor','zs', file=outfile)
  }

  ## save fft results to table files
  if ( savecsv ) {
    outfile <- paste(dft.name,".csv",sep="")
    if ( verb )
      cat(paste("writing results to files:      \t",date(),"\nfile:\t",
                outfile, "\n"), file=stdout())
    fdat <- data.frame(coor,fdat)
    write.table(fdat,outfile,sep="\t",row.names=FALSE,quote=FALSE,na="")
    if ( perm>0 ) {
      outfile <- paste(dft.name,"_pvalues","_p",perm,".csv",sep="")
      if ( verb ) cat(paste("file:\t", outfile, "\n"))
      pdat <- data.frame(coor,pdat)
    write.table(pdat,outfile,sep="\t",row.names=FALSE,quote=FALSE,na="")
    }
  }
}
  
  
## TODO: move this to external file and save memory!?
## OOR:
## 0) pre-filter informative DFT components
## 1) do PCA, and use this as similarity measure
## 2) do flowClust only for subset(s) of data, and assign
##    rest by similarity
## 2a) do flowClust only for p-value < 0.01, cluster rest differently
## TODO: calculate DAG or similarity between clusters and re-sort
## clusters (by phase?)
## TEST: using GDH3 domain, chrI:30000..35000, 20000..40000

if ( max(cls) > 1 ) {
  if ( verb ) {
    cat(paste("clustering:\t", date(),"\n"),file=stdout())
    cat(paste("clusters:\t", paste(unique(range(cls)),collapse=" - "),
              "\n"),file=stdout())
    if ( ncpu>1 ) 
      cat(paste("using", ncpu-1, "additional cores\n"), file=stdout())
  }
  ## set number of cores (-1 since one is used for the main process)
  if ( ncpu>1 ) {
    library("parallel")
    options(cores=ncpu-1)
  } else {
    options(cores=1)
  }
  require("flowClust")

  ## DFT components to use
  if ( is.infinite(dft[1]) ) {
    ## default: don't use DC and Nyquist
    dftCom <- 2:(ncol(fdat)-1) #2:4 # 2:5 #4 #2:12
  }  else {
    ## via command-line, specify all components to use
    dftCom <- dft
  }
  if ( verb )
    cat(paste("components:\t", paste(dftCom,collapse=","),"\n"),file=stdout())

  ## chromosome coordinates of non-zero read-count nt
  chrRange <- coor[!zs,,drop=FALSE]

  ## cluster only range
  ## TODO: both strands are used here! select by coordinate order
  ## TODO: allow use of direct indexing
  chrIdx <- rep(TRUE, sum(!zs))
  if ( rng!="" ) {
    if ( verb ) cat(paste("clustering only range:", rng, "\n"))
    #chr <- 1
    #pos <- c(20000,40000) # c(30000,31000) #
    #pos <- c(1,Inf) # c(30000,31000) #
    ## select strand
    str <- 1;
    if ( pos[1]>pos[2] ) {
      str <- -1
      pos <- rev(pos)
    }
    ## get range
    chrIdx <- (chrRange[,1] == chr & chrRange[,3]==str) & (chrRange[,2] >= pos[1] & chrRange[,2] <= pos[2])
  }

  ## get data, real and imaginary part of the DFT
  clsDat <- fdat[!zs,][chrIdx,]
  clsDat <- cbind(Re(clsDat[,dftCom]),Im(clsDat[,dftCom]))


  ## safe pre-cluster data (for testing, rm later)
  #file.name <- paste("genomeOscillation_",cl.name,"_precluster.RData",sep="")
  #save.image(file=file.name)

  
  ## CLUSTER
  numCls <- min(cls):max(cls)
  ## default flowClust params
  ## B is now command-line param
  #B <- 500 # max. num. of EM iterations # TODO: select lower for genome scan!
  tol <- 1e-5 # tolerance for EM convergence
  lambda <- 1 # intial Box-Cox trafo
  nu <- 4 # Inf for pure Gaussian, initial Box-Cox trafo
  nu.est <- 0 # 0: no, 1: non-specific, 2: cluster-specific estimation of nu
  trans <- 1 # 0: no, 1: non-specific, 2: cluster-specific estimation of lambda
  ## override
  if ( gauss ) { nu <- Inf; trans <- 1}
  fcls <- flowClust::flowClust(clsDat, K=numCls,
                               B=B, tol=tol, lambda=lambda,
                               nu=nu, nu.est=nu.est, trans=trans)
  if ( verb )
    cat(paste("finished:\t", date(),"\n"),file=stdout())

  ## to be safe, store here after long calculation (for testing, rm later)
  #file.name <- paste("genomeOscillation_",cl.name,"_clustered.RData",sep="")
  #save.image(file=file.name)

  ## collect clusterings
  cluster.matrix <- matrix(NA, nrow=nrow(clsDat), ncol=length(numCls))
  colnames(cluster.matrix) <- as.character(numCls)
  bic <- rep(NA, length(numCls))
  names(bic) <- as.character(numCls)
  for ( i in 1:length(fcls) ) {
    if ( length(fcls) > 1 ) fc <- fcls[[i]]
    else fc <- fcls
    cl.num <- as.character(fc@K)
    #cat(paste(i, cl.num, "\n"))
    cluster <- flowClust::Map(fc,rm.outliers=F)
    cluster.matrix[, cl.num] <- cluster
    bic[cl.num] <- fc@BIC
  }
  ## TODO 20160424 : catch failure !
  bic <- bic[!is.na(bic)]
  max.bic <- names(bic)[which(bic==max(bic))]

  ## big clustering results matrix
  rangeClusters <- cbind(chrRange[chrIdx,],cluster.matrix)

  ## save all clustering results to RData file
  ## TODO: save same data as _collectClusterings, to be read by _sortClusters
  # only matrix 'max.clb'
  file.name <- paste("genomeOscillation_",cl.name,"_clustering.RData",sep="")
  save('rangeClusters', 'fcls','bic','clsDat', file=file.name)

  ## clustering data table
  file.name <- paste("genomeOscillation_",cl.name,"_clustering.csv",sep="")
  write.table(rangeClusters, file.name, sep="\t",row.names=F,na="")

  ## BIC Plot
  if ( length(bic) > 1 ) {
    file.name <- paste("genomeOscillation_",cl.name,"_BIC.png",sep="")
    png(file.name)
    plot(names(bic), bic)
    dev.off()
  }
    #str <- rangeClusters[,"strand"]==1
  #png("genomeOscillation_DFT2_5_k25.png")
  #plot(rangeClusters[str,"coor"],rangeClusters[str,max.bic],type="l",axes=F,xlim=c(30000,40000),ylab="cluster", xlab="genome position")       
  #axis(1,at=seq(30000,40000,500),cex=.5);axis(2);
  #dev.off()
  
}

cat(paste("DONE:\t", date(), "\n"))
quit(save="no")
  
  

# End of file
