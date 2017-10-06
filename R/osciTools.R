
#' Discrete Fourier Transformation
#' 
#' A simple wrapper around \code{\link[stats]{mvfft}} to perform
#' Discrete Fourier Transformation, discard redundant components
#' (for real numbers only, the second half of the transform),
#' and name columns by `DC` (direct current, a term from electrical
#' engineering) and the numbers of cycles in the data set reflected
#' by each component. The last column is the Nyquist frequency.
#' @param x data to be transformed
#' @export
get_fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n]
    if ( n==1 )
        colnames(fft) <- "DC" # is this necessary or better stop()?
    else
        colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}

#' convert complex numbers to degrees
#'
#' convert complex number representation of polar coordinates
#' (such as provided by \code{\link[stats]{mvfft}} or segmenTools'
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

#' calculate peak phase
#' 
#' calculate phases via a Discrete Fourier Transformation, but starting
#' at the first time point, such that phase 0/360 indicates peak
#' at time 0, and phase 180 indicates peak at half of the cycle
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

    if ( class(x)=="timeseries" ) # segmenTier class !
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
#' @param xlim x-axis range to show, defaults to one period
#' @seealso \code{\link{calculatePhase}}
#' @export
testPhase <- function(n=5, cyc=4, T=24, res=6, xlim=c(-T/res,T+T/res)) {

    ## time vector
    time <- seq(0,cyc*T,T/res)

    ## generate cosine waves where phase 0 means peak at origin t0
    phases <- seq(0,pi,length.out=n) # phase shift from time 0 in pi
    amps <- abs(stats::rnorm(n,1,.3)) # amplitude around 1
    
    y <- matrix(NA,ncol=length(time),nrow=length(phases))
    
    for ( i in 1:nrow(y) )
        y[i,] <- amps[i]*cos(2*pi*time/T - phases[i])
     
    ## phase(dat,cyc) vs. phases*180/pi
    ## TODO: why is this slightly shifted (.75 degree) wrt original phases?
    rephases <- calculatePhase(y,cyc)[,1]

    matplot(time,t(y),type="l",col=1:nrow(y),xlim=xlim,lty=1)
    abline(v=seq(0,cyc*T,T),lty=3,col="gray");abline(h=0,lty=3,col="gray") 
    abline(v=rephases*T/360,col=1:nrow(y)) ## plot re-covered phases
}

