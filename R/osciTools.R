
#' Discrete Fourier Transformation
#' 
#' perform Discrete Fourier Transformation using \code{mvfft} from
#' the \code{stats} package, and returning the non-redundant (for
#' real numbers) first half of the transform, i.e., from the DC
#' (direct current) component to the Nyquist frequency
#' @param x data to be transformed
#' @export
get.fft <- function(x) {
    n <- floor(ncol(x)/2) +1 ## Nyquist-freq
    fft <- t(stats::mvfft(t(x)))[,1:n]
    if ( n==1 )
        colnames(fft) <- "DC"
    else
        colnames(fft) <- c("DC",as.character(1:(n-1)))
    fft
}

#' calculate peak phase
#' 
#' calculate phases via a Discrete Fourier Transformation, but starting
#' at the first time point, such that phase 0/360 indicates peak
#' at time 0, and phase 180 indicates peak at half of the cycle
#' @param x a time-series matrix with columns as time points,
#' or a \code{timeseries} object as returned by \code{\link{processTimeseries}}
#' @param cycles optional integer vector for the DFT components
#' (number of cycles in the data) for which phases are to be calculated;
#' all are returned if \code{cycle} is not specified
#' @param degrees logical to indicate whether phases should be reported
#' in degrees (0 - 360) or radians (0 - 2*pi)
#' @export
phase <- function(x, cycles, degrees=TRUE) {

    if ( class(x)!="timeseries" ) {
        if ( class(x)!="matrix" )
            x <- data.matrix(x)
        x <- processTimeseries(x, use.fft=TRUE, use.snr=FALSE)
    }
    if ( missing(cycles) )
        cycles <- 2:ncol(x$dft)
    else cycles <- cycles +1

    phase <- x$dft[,cycles,drop=FALSE]
    ## TODO: check whether sign is correct?
    ## atan2 gives phase shifts relative to a cosine, from -pi (delay)
    ## to pi (advance) wrt peak
    phase <- atan2(Im(phase),Re(phase)) * 180/pi
    ## atan2: negative phase means shift to right
    ## positive phase means shift to left
    ## re-define to indicate peak time wrt time t0=0 to t=Period
    ## 1) switch sign, positive means later!
    phase <- - phase 
    ## adjust phase angles,
    ## such that 0 is at the beginning of data
    phase <- ifelse(phase <=  0, phase + 360, phase)
    phase <- ifelse(phase > 360, phase - 360, phase)
    if ( !degrees ) phase <- phase/360 * 2*pi
    phase
}

#' tests phase recovery by \code{\link{phase}}
#'
#' @param n number of cosines to generate
#' @param cyc numbers of cycles to generate
#' @param T period
#' @param res resolution: samples per period
#' @param xlim x-axis range to show, defaults to one period
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
    rephases <- phase(y,cyc)[,1]

    matplot(time,t(y),type="l",col=1:nrow(y),xlim=xlim,lty=1)
    abline(v=seq(0,cyc*T,T),lty=3,col="gray");abline(h=0,lty=3,col="gray") 
    abline(v=rephases*T/360,col=1:nrow(y)) ## plot re-covered phases
}

