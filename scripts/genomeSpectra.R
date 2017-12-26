#!/usr/bin/Rscript

## reads a genomeData file (chr\tcoor\t\value[1]\tvalue[2]\t...\tvalue[n])
## with time series 1..n for each coordinate,
## and calculates oscillation properties of each time series.

verb <- FALSE 
plot <- FALSE
save <- FALSE

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
if ( !exists("mx",mode="character") ) { mx <- Inf
} else { mx <- as.numeric(mx) }
## plot range: plot zoomed
if ( !exists("rng",mode="character") ) { rng <- NA
} else { rng <- as.numeric(unlist(strsplit(rng,","))) }
## test range: find maximal amplitude
if ( !exists("test",mode="character") ) { test <- NA
} else { test <- as.numeric(unlist(strsplit(test,","))) }
## log plot? directly passed to plot ("", "xy", "y", "x")
if ( !exists("log",mode="character") ) log <- ""

verb <- as.logical(verb)
plot <- as.logical(plot)
save <- as.logical(save)

## utils
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100) {
  file.name <- paste(file.name, type, sep=".")
  if ( type == "png" )
    png(file.name, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    postscript(file.name, width=width, height=height, paper="special")
  if ( type == "pdf" )
    pdf(file.name, width=width, height=height)
}



## load data, either RData or tab-delimited file
if ( length(grep(".RData$", i)) > 0 )  {
  load(i) ## RData
} else { ## tab-delim file
  dat <- read.table(i,sep="\t",header=TRUE)
  dat <- as.matrix(dat)
  dat[is.na(dat)] <- 0
}

## split data and coordinates
coor <- dat[,cc,drop=FALSE] # is not used anymore, may be skipped!
dat <- dat[,!colnames(dat)%in%cc,drop=FALSE] ## rm coordinates

## TODO: use only some columns
#dat <- dat[,1:11,drop=FALSE]

## PERIODS
nyq <- floor(nrow(dat)/2)
periods <- c(0,nrow(dat)/(1:(nyq-1)))

## DFT
spectra <- matrix(NA,ncol=ncol(dat),nrow=nyq)
colnames(spectra) <- colnames(dat)
for ( i in 1:ncol(dat) ) {
  id <- sub("\\.","_",colnames(dat)[i])
  if ( verb ) cat(paste("DFT of column", id))

  ft <- fft(dat[,i])[1:nyq]
  amps <- abs(ft)/nrow(dat)
  spectra[,i] <- amps

  if ( verb ) cat(".\n")
}


## save results
if ( save ) {
    file.name <- paste(o,"_genomeSpectra.RData",sep="")
    if ( verb ) cat(paste("saving results", file.name, "\n"))
    save('spectra','periods','avgamp',file=file.name)
}

if ( !plot ) quit(save="no")

### DC and average amplitudes:
## TODO: this depends on p/m strand annotation, get rid of this

## average period vs. DC
avgamp <- apply(spectra[2:nrow(spectra),,drop=FALSE],2,mean)
dc <- spectra[1,]

## NOTE: mirror of DC component for forward vs. reverse strand
## also absent and general difference gets smaller without rRNA
file.name <- paste(o,"_genomeSpectrum_DC")
plotdev(file.name,type="png",res=300,width=5,height=3)
par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.5,0))
plot(dc[grep("m$",names(dc))],type="l",col=2,ylim=range(dc),
     ylab="DC component",xlab="sample")
lines(dc[grep("p$",names(dc))],type="l",col=1)
legend("topleft",legend=c("plus","minus"),col=1:2,lty=1)
dev.off()
file.name <- paste(o,"_genomeSpectrum_avgAmp")
plotdev(file.name,type="png",res=300,width=5,height=3)
par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.5,0))
plot(avgamp[grep("m$",names(avgamp))],type="l",col=2,ylim=range(avgamp),
     ylab="average amplitude",xlab="sample")
lines(avgamp[grep("p$",names(avgamp))],type="l",col=1)
legend("topleft",legend=c("plus","minus"),col=1:2,lty=1)
dev.off()

## maximal period in total range plot
if ( is.infinite(mx) ) mx <- periods[2]
totrng <- periods > 0 & periods <= mx

## plot range
if ( is.na(rng[1]) ) rng <- 1:max(periods[periods <= mx])
prdrng <- periods >= rng[1] & periods <= rng[2] ## plot and store range

## range to get maximal amplitude location
if ( is.na(test[1]) ) test <- 1:max(periods[periods <= mx])
testrng <- periods >= test[1] & periods <= test[2]


ylim_rng <- c(0,max(spectra[prdrng,],na.rm=TRUE))
ylim_all <- c(0,max(spectra[totrng,],na.rm=TRUE))
if ( length(grep("y",log))>0 ) {
    ylim_rng[1] <- min(spectra[prdrng,],na.rm=TRUE)
    ylim_all[1] <- min(spectra[totrng,],na.rm=TRUE)
    if ( ylim_rng[1] == 0 ) ylim_rng[1] <- 1e-10
    if ( ylim_all[1] == 0 ) ylim_all[1] <- 1e-10
}

for ( i in 1:ncol(spectra) ) {
    id <- sub("\\.","_",colnames(spectra)[i])
    amps <- spectra[,i]
    
    
    ## plot total range (period <= mx)
    file.name <- paste(o,"_",id,"_genomeSpectrum_all",
                       ifelse(log=="","",paste("_log",log,sep="")),sep="")
    plotdev(file.name,type="png",res=300,width=5,height=3)
    par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.5,0))
    plot(periods[totrng], amps[totrng],type="h",xlab="period [bp]",ylab="amplitude",col="gray",log=log,ylim=ylim_all)
    dev.off()

    
    ## plot of selected range, with maximal period indicated
    if ( !is.na(rng[1]) ) {
        mxpk <- max(amps[testrng])
        peak <- which(amps[testrng]==mxpk)
        peak <- periods[testrng][peak]
        
        file.name <- paste(o,"_",id,"_genomeSpectrum_range",
                           ifelse(log=="","",paste("_log",log,sep="")),sep="")
        if ( verb ) cat(paste("plotting", file.name, "\n"))
        plotdev(file.name,type="png",res=300,width=5,height=3)
        par(mai=c(.5,.5,.1,.1),mgp=c(1.3,.5,0))
        plot(periods[prdrng], amps[prdrng],type="h",xlab="period [bp]",ylab="amplitude",col="gray",log=log,ylim=ylim_rng)
        points(peak,mxpk)
        text(peak,mxpk,labels=paste(round(peak,1),"bp"),pos=4)
        dev.off()
    }
}


quit(save="no")
