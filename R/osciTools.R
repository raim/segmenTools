
#' Discrete Fourier Transformation
#' 
#' A simple wrapper around \code{\link[stats:fft]{mvfft}} to perform
#' Discrete Fourier Transformation, discard redundant components
#' (for real numbers only, the second half of the transform),
#' and name columns by `DC` (direct current, a term from electrical
#' engineering) and the numbers of cycles in the data set reflected
#' by each component. The last column is the Nyquist frequency.
#' @param x data to be transformed
#' @export
get_fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n,drop=FALSE]
    if ( n==1 )
        colnames(fft) <- "DC" # is this necessary or better stop()?
    else
        colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}

## TODO: simple routine for windowed FFT, see yeastSeq2016 DO analysis
fft_window <- function(x, win, Tmax, Tmin, res=1) {

    ## calculate window size <-> Tmax (one is required)
    if ( missing(win) & !missing(Tmax) )
        win <- Tmax/res
    else if ( missing(Tmax) & !missing(win) )
        Tmax <- win*res

    ## call fourier in windows, via apply
    last <- length(x)-win+1
    dft <- sapply(1:last, function(i) fft(x[i:(i+win-1)]))
    dft <- dft[1:(floor(win/2) +1),] # rm mirror
    dft <- dft/win # scale amplitudes
    dc <- abs(dft[1,]) # DC component
    dft <- dft[2:nrow(dft),] # rm DC component
    Tosc <- win/(1:(nrow(dft))) # period in index time unit
    Tcenter <- 1:last + win/2 # center of fourier window in index time unit

    ## convert by resolution info
    Tosc <- Tosc*res # period in hours
    Tcenter <- Tcenter*res  # center of window in hours
    
    ## cut to minimal period
    rmp <- Tosc<Tmin
    dft <- dft[!rmp,]
    Tosc <- Tosc[!rmp]

    ## TODO: clarify, Tcenter?
    list(time=Tcenter, period=Tosc, DC=dc, DFT=dft)
}

#' convert phase (degree) to polar coordinates
#' @param phase phase in degree
#' @param amplitude amplitude
#' @param unit "degree" or "rad"
#' @seealso \code{\link{complex2degree}}
#' @export
degree2complex <- function(phase, amplitude=1, unit="degree"){
    ## to radian
    if ( unit=="degree" )
        phase <- phase*pi/180
    xy.coords(x=amplitude*cos(phase),
              y=amplitude*sin(phase))
}


#' convert complex numbers to degrees
#'
#' convert complex number representation of polar coordinates
#' (such as provided by \code{\link[stats:fft]{mvfft}} or segmenTools'
#' wrapper \code{\link{get_fft}}
#' to degrees, and map to 0-360 representing the first cycle
#' in an experiment
#' @param x a vector or matrix of complex numbers
#' @seealso \code{\link{get_fft}}
#' @export
complex2degree <- function(x) {

    ## atan2 gives phase shifts relative to a cosine,
    ## from -pi (delay) to pi (advance) wrt peak;
    ## to map phase on actual time, we need to change the sign ...
    ph <- - atan2(Im(x),Re(x)) * 180/pi
    ## and shift the negative phases ahead
    ph <- ifelse(ph<=  0,ph + 360, ph)
    ph <- ifelse(ph> 360,ph - 360, ph) # TODO: required?
    ph
}

#' circular density
#'
#' calculates kernel density estimates for circular data (`phase angles').
#' It is a convenience wrapper around the
#' \code{\link[circular:density.circular]{density.circular}}
#' function from package \code{circular} that (a) coerces results to numeric
#' types to avoid automatic handling of the S3 class circular data
#' (eg. polar plots), (b) allows to scale densities by the absolute number
#' (argument \code{freq}), and (c) can alternatively invoke
#' \code{\link[stats:density]{density}} from the \code{stats} package
#' (see option \code{high.mem}) for long phase vectors (`genome-wide')
#' due to high memory usage of
#' \code{\link[circular:density.circular]{density.circular}}.
#' NOTE: results are equivalent only for the scaled version freq=TRUE;
#' TODO: find out why densities differ by a factor of ~100 between the 2
#' @param x phases in degrees and of class numeric
#' @param bw the smoothing bandwith parameter
#' @param n the number of equally spaced points at which density
#' will be estimated
#' @param freq if TRUE densities \code{y} will be scaled by the total
#' number of measurement \code{N} as \code{N * y/sum(y)}; the resulting
#' density function will integrate to \code{N}
#' @param units phase angle units, 'degrees' or 'radians'
#' @param high.mem use \code{\link[stats:density]{density}} to calculate
#' kernel densities, based on copying data from -360 to 0 and 360 to
#' 720 degrees
#' @param na.rm remove NA values before calculation
#' @param ... further arguments to
#' \code{\link[circular:density.circular]{density.circular}}
#' or \code{\link[stats:density]{density}} 
#' @seealso \code{\link[circular:density.circular]{density.circular}},
#' \code{\link[stats:density]{density}} 
#'@export
circ.dens <- function(x, bw=18, n=512, freq=FALSE, units="degrees", high.mem=FALSE, na.rm=TRUE, ...) {

    if ( na.rm ) x <- x[!is.na(x)]

    ##x <- 2*pi * x/360 # in radian # TODO allow/use radian?
    if ( units=="degrees" )
        maxc <- 360
    else if ( units=="radians" )
        maxc <- 2*pi
    
    if ( high.mem ) { # use stats:density
        x <- c(x-maxc, x, x+maxc)
        cd <- stats::density(x, bw=bw, n=n, from=0, to=maxc, ...)
    } else { # use circular:density
        ph <- circular::circular(x, units=units)
        cd <- circular::density.circular(ph, bw=bw, n=n, ...)
    }
    ## results
    xx <- as.numeric(cd$x)

    ## normalize to total density 1 and multiply with N
    ## NOTE: after normalization results are comparable!
    if ( freq ) 
        cd$y <- length(x) * cd$y/sum(cd$y)
    
    cd <- list(x=xx, y=cd$y, bw=bw)
    cd
}


#' calculate peak phase
#' 
#' calculate phases via a Discrete Fourier Transformation, but starting
#' at the first time point, such that phase 0/360 indicates peak
#' at time 0, and phase 180 indicates peak at half of the cycle
#' TODO: implement systematic phase error for time series that
#' do not cover multiples of full cycles; phase error related to
#' pi/resolution, see comments in testPhase
#' @param x a time-series matrix with columns as time points,
#' or a \code{timeseries} object as returned by
#' \code{\link[segmenTier]{processTimeseries}}
#' @param cycles optional integer vector for the DFT components
#' (number of cycles in the data) for which phases are to be calculated;
#' all are returned if \code{cycle} is not specified
#' @param degrees logical to indicate whether phases should be reported
#' in degrees (0 - 360) or radians (0 - 2*pi)
#' @seealso \code{\link{get_fft}} and \code{\link{complex2degree}}
#' @export
calculatePhase <- function(x, cycles, degrees=TRUE) {

    if ( inherits(x, "timeseries") ) # segmenTier class !
        dft <- x$dft
    else dft <- get_fft(data.matrix(x))

    ## all cycles:
    if ( missing(cycles) )
        cycles <- 2:ncol(dft)
    else cycles <- cycles +1

    ## get requested cycle numbers
    phase <- dft[,cycles,drop=FALSE]
    ## convert complex number representation to degrees
    phase <- complex2degree(phase)
    
    ## convert back to radian if requested
    if ( !degrees ) phase <- phase/360 * 2*pi
    
    phase
}

#' tests phase recovery by \code{\link{calculatePhase}}
#'
#' generates cosine waves with fixed periods but random
#' amplitudes and systematic phase shifts and sampled
#' at varying resolution to test phase recovery performance
#' via Discrete Fourier Transformation in \code{\link{calculatePhase}}.
#' TODO: allow random phases, noise on period and general
#' addition of noise
#' @param n number of cosines to generate
#' @param cyc numbers of cycles to generate
#' @param T period
#' @param res resolution: samples per period
#' @param xlim x-axis range to show, defaults to one period \code{+/- T/res}
#' @seealso \code{\link{calculatePhase}}
#' @export
testPhase <- function(n=5, cyc=4, T=24, res=6, xlim=c(-T/res,T+T/res)) {

    ## time vector
    ## NOTE: omitting last timepoint (-T/res)
    time <- seq(0,cyc*T-T/res,T/res)

    ## generate cosine waves where phase 0 means peak at origin t0
    phases <- seq(0,2*pi,length.out=n) # phase shift from time 0 in pi
    amps <- abs(stats::rnorm(n,1,.3)) # amplitude around 1
    
    y <- matrix(NA,ncol=length(time),nrow=length(phases))
    
    for ( i in 1:nrow(y) )
        y[i,] <- amps[i]*cos(2*pi*time/T - phases[i])
     
    ## phase(dat,cyc) vs. phases*180/pi
    ## TODO: calculate amplitudes and draw full reconstructed sine!
    ## TODO: calculate phase for time series that are not multiples of
    ## full cycles; phase error for e.g. one timepoint to much
    ## seems to be pi/res
    rephases <- calculatePhase(y,cyc)[,1]
    
    matplot(time,t(y),type="l",col=1:nrow(y),xlim=xlim,lty=1)
    axis(3,at=time,labels=NA)
    mtext("samples", 3, 1)
    abline(v=seq(0,cyc*T,T),lty=3,col="gray");abline(h=0,lty=3,col="gray") 
    abline(v=rephases*T/360,col=1:nrow(y)) ## plot re-covered phases
    tmp <- cbind(orig=phases,pred=rephases*2*pi/360)
    ##plot(tmp,ylim=c(0,2*pi),xlim=c(0,2*pi),xaxs="i",yaxs="i");abline(a=0,b=1)
}

## time-series processing
## function taken from FLUSH.LVS.bundle/R/normalize_LVS.R and via tataProject.R
## http://www.meb.ki.se/~yudpaw/
## http://www.ncbi.nlm.nih.gov/pubmed/18318917
#' LVS/LOS normalization
#'
#' This function normalizes the columns of a matrix by a set of
#' invariant data rows (\code{LOS}, a "least-oscillating set") to be
#' defined by the user (option \code{los.id}). For each column a LOESS
#' fit, local fitting of a 2nd order polynomial with default
#' parameters of the R \code{stats} function
#' \code{\link[stats:loess]{loess}} is calculated between the values
#' of the reference set (\code{LOS}) and their median or median
#' (option \code{norm.ref} over all columns.  This is a modified
#' version of the function `normalize_LVS` from package
#' \code{FLUSH.LVS} (Calza et al. 2008), as used for normalization by
#' a "least-oscillating" set of "microarray" probes from a
#' transcriptome time-series of "yeast respiratory oscillations" by
#' Machne & Murray, 2012. Alternatively to LOESS but un-tested, the 
#' data can also be normalized by smoothing or simple min/max scaling.
#' Also un-tested and highly experimental: to avoid cutting data at
#' highest values of the \code{LOS} set, one can \code{add.maxima} of each
#' column to the LOESS fits, or after the fit and normalization reset
#' (\code{reset.extrema}) or re-scale (\code{extend.max}) maximal values.
#' @param x a data matrix to be normalized
#' @param los.id a numerical or logical vector, providing the indices
#' or a TRUE/FALSE vector for a set of "least-oscillating" or "least-variant"
#' data rows
#' @param norm.ref take the "median" or the "mean" of the reference
#' data set (\code{los.id}) for normalization
#' @param norm.type use "loess", "smooth" or "scale" normalization; see Details
#' @param add.maximum untested, see Details
#' @param reset.extrema untested, see Details
#' @param extend.max untested, see Details
#' @references
#' Calza S, Valentini D, Pawitan Y (2008)
#' Normalization of oligonucleotide arrays based on the least-variant
#' set of genes. BMC Bioinformatics 9:140.
#'
#' Machne R, Murray DB (2012)
#' The yin and yang of yeast transcription: elements of a global feedback
#' system between metabolism and chromatin. PLoS One 7(6):e37906
#' 
#' @export
normalize.los <- function(x, los.id,norm.ref=c("median","mean"),
                          norm.type=c("loess","smooth","scale"),
                          add.maximum=FALSE, reset.extrema=FALSE,
                          extend.max=FALSE) {
    
    if(missing(los.id))
        stop("must provide index (logical or numerical) for LOS genes")
    
    if(!is.logical(los.id) && !is.numeric(los.id))
        stop("los.id must be a vector either logical or numeric")
    
    ## generate reference data: mean/median for each reference probe
    norm.ref <- match.arg(norm.ref)
    ref.data <- switch(norm.ref,
                       "median"=apply(x,1,median),
                       "mean"=apply(x,1,mean))
    
    out <- matrix(NA, ncol=ncol(x), nrow=nrow(x))

  ## SIMPLE SCALING
  if ( norm.type[1]=="scale" ) {
    out <- t(t(x) * apply(t(ref.data[los]/x[los,]),1,mean))
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    return(out)       
  }
    
  for ( i in 1:ncol(x) ) {
    los <- los.id
    ## add all individual array maxima to extend fitting curve
    if ( add.maximum ) {
      if ( is.logical(los)) los[x[,i]>=max(x[los,i])]<-TRUE
      if ( is.numeric(los)) los<-c(los,which(x[,i]>=max(x[los,i])))
    }
    if ( norm.type[1]=="loess" ) {
      ## fit polynomial of degree 2 for reference data points
      sm <- loess(x[los,i]~ref.data[los], degree = 2)

      ## sm$x IS ref.data[los] , i.e. time-series medians/means
      ## sm$fitted IS x[los,i] curve??
            
      ## linearly interpolate fitted curve at raw data points,
      ## rule=2: extrapolate at maxima to last value (rule=1:
      ## extrapolate as NA)
      a <- approx(x=sm$fitted, y=sm$x, xout=x[,i], rule = 2)$y
    } else if ( norm.type[1]=="smooth") {
      ## generate a smoothed spline for reference data points
      sm <- smooth.spline(y=x[los,i] , x= ref.data[los])
      ## as above
      a <- approx(x=sm$y, y=sm$x, xout=x[,i], rule = 2)$y
    }
    
    ## HANDLE DATA BEYONG APPROX RANGE, x>=max(los)
    if ( reset.extrema ) {
      
      ## TODO: instead use simple scaling with highest LOS value
      ##max.los <- which(los & x[,i] == max(x[los,i]))
      ##a[x[,i]>=max(x[los,i])] <- x[x[,i]>=max(x[los,i]),i] * a[max.los]/x[max.los,i]))
      
      ## TODO: was >= and <= correct(er) ??
      a[x[,i] > max(x[los,i])] <- x[x[,i] > max(x[los,i]),i]
      a[x[,i] < min(x[los,i])] <- x[x[,i] < min(x[los,i]),i]  
    } else if ( extend.max ) {
      ## use simple scaling with highest LOS value
      max.los <- which(los & x[,i] == max(x[los,i]))
      a[  x[,i] > max(x[los,i]) ] <-
        x[x[,i] > max(x[los,i]),i] * a[max.los]/x[max.los,i]
      ## TODO: x[,i]<=min(x[los,i]) !!??
      
    }
    out[,i] <- a
  }
    
  colnames(out) <- colnames(x)
  rownames(out) <- rownames(x)
  
  return(out)
}
