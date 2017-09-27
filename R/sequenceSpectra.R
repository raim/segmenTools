
## This is a code-base initiated by Robert Lehmann for
## the publication Lehmann, Machne, Herzel in Nucl Acids Res 2014.
## It is used to calculate periodicities of dinucleotide patterns
## in genomes.

################################################################
# takes DNAString sequence and a set of patterns, returns a sequence of 0s, with 1s marking pattern occurrences
# seq - DNAString input sequence
# patterns - list of characters describing the hit patterns
make.hit.seq <- function(seq, patterns) {
	genome.dn.pos <- lapply(patterns, function(x) matchPattern(x, seq, fixed=FALSE))
	names(genome.dn.pos) <- patterns
	#add hits for patterns in sequence
	pat.pos <- rep(0,length(seq))
	for(p in patterns) 
          pat.pos[ unlist( start(genome.dn.pos[[p]]) ) ] <- 1
	pat.pos
} 
#######################################
# calculates autocorrelation function on 01-series for the specified range
# ts - series consisting of [0,1] marking pattern hit positions
# d.range - range in which to compute the autocorrelation function
acf.hph <- function(ts, d.range) {
	N <- length(ts)
	Naa <- rep(NA, length(d.range))
	#assuming that the sequence consists only of 0 and 1, sum(ts) yields the total number of pattern hits 
	pa <- sum(ts)/N
	#iter over distances
	for(d in d.range){
		j <- 1:(N-d)
		Naa[d] <- sum(ts[j] == 1 & ts[j+d] == 1)
	}
	## calculate autocorrelation
	Naa <- Naa/(N-d.range) - pa^2
	Naa
}
################################################################
# calculates the power spectrum for the provided signal
power.spect.mean <- function(signal, Periods, norm.spect=FALSE, s) {
	mi <- complex(real = 0, imaginary = -1)
	L <- length(signal)
        if ( missing(s) )
          s <- c(1:(L))

         
	amplitudes <- sapply(Periods, function(p) abs(sum(signal * exp(mi * s * (2*pi)/p))))
	names(amplitudes) <- Periods
	
	if(norm.spect) { # normalized spectra if specified
		mean.amp <- mean(amplitudes)
		if(mean.amp != 0)
			amplitudes <- amplitudes /  mean.amp
	}
	amplitudes
}

################################################################
# power spectrum after Mrazek
power.spect.m <- function(signal, Periods, norm.spect=FALSE, s) {
	mi <- complex(real = 0, imaginary = -1)
	L <- length(signal)
        if ( missing(s) )
          s <- c(1:(L))

        #s <- acf2spectrum.range[1]:acf2spectrum.range[2]
	amplitudes <- sapply(Periods,
                             function(p)
                             abs(sum(signal * exp(mi * s * (2*pi)/p))))
	names(amplitudes) <- Periods
	
	if(norm.spect) # normalized spectra if specified
          amplitudes <- ((max(Periods)-min(Periods)+1) * amplitudes) / sum(amplitudes)
        amplitudes
}
################################################################
# power spectrum after Kravatskaya
power.spect.k <- function(signal, Periods, norm.spect=FALSE) {
	mi	<- complex(real = 0, imaginary = -1)
	L	<- length(signal)
	n	<- 0:(L-1)
	if(is.null(Periods)) {
		qn	<- 2 * pi * n / L
		names(qn) <- 1/n * 100
		names(qn)[1] <- 0
		qn <- qn[order(names(qn),decreasing=TRUE)]
	} else {
		qn	<- 2 * pi * 1/Periods
		names(qn) <- Periods
	}
	m 	<- 1:L
	p.spect <- sapply(qn, function(p)
				abs( L^(-1/2) * sum( signal * exp(mi * p * m)) )
	)
	names(p.spect) <- names(qn)
	if(norm.spect) { # normalize if requested
		Nh <- sum(signal) # number of hits in sequence
		if(Nh != 0) # when no pattern hit is found: avoid division by zero 
			p.spect <- p.spect /  (Nh * (L - Nh) / L * (L-1))
	}
	p.spect
}

#######################################
#smooth vector function
moving.avg <- function(x,n=5, na.rm=TRUE){
	res = filter(x,rep(1/n,n), sides=2)
	if(na.rm)
			res[-which(is.na(res))]
	else
		res
}
#######################################
#smooth matrix function
moving.avg.mat <- function(x,dim=1, n=5, na.rm=TRUE){
	if(n==0) {
		warning('n=0 would provide empty martix, returing input mat instead...')
		return(x)
	}
	if(dim==1) { # for some funy reason, apply transposes matrix when dim=1
		res = t(apply(x, dim, moving.avg, n=n, na.rm=FALSE))
		colnames(res) = colnames(x)
		if(na.rm)
			res = res[,-which(is.na(res[1,]))]
	} else {
		res = apply(x, dim, moving.avg, n=n, na.rm=FALSE)
		rownames(res) = rownames(x)
		if(na.rm)
			res = res[-which(is.na(res[,1])),]
	}
	res
}
################################################################
# calculates the periodicity spectrum of the provided genomic sequence 
seq.periodicity.spectrum <- function(in.seq, patterns, smooth.acf=TRUE, 
		acf.range, acf2spectrum.range, use.powerspect=TRUE, powerspect.periods, norm.spectrum='n') {
  if(!class(norm.spectrum) == 'character') {
    warning('invalid value for "norm.spectrum", found class ',class(norm.spectrum),', required character!\n')
    return()
  }
  ##find pattern hits
  ww.pos <- make.hit.seq(in.seq, patterns)
  ## calc. autocorrelation function on hit seq and take specified section
  if(smooth.acf)
    acf.section <- moving.avg(acf.hph(ww.pos, 1:acf.range),n=3, na.rm=TRUE)[acf2spectrum.range[1]:acf2spectrum.range[2]]
  else
    acf.section <- acf.hph(ww.pos, 1:acf.range)[acf2spectrum.range[1]:acf2spectrum.range[2]]
  ## calulcate amplitude spectrum of section of autocorrelation function
  if(!use.powerspect) {
    amps <- amplitude.spect(acf.section)
    if(!is.null(powerspect.periods))
      return(amps[which(as.numeric(names(amps)) >= min(powerspect.periods) &  
                        as.numeric(names(amps))>= max(powerspect.periods))])
    else
      return(amps)
    
  } else { 
    switch(norm.spectrum,
           'n'={return(power.spect.m(signal=acf.section, Periods=powerspect.periods, norm.spect=FALSE, s=acf2spectrum.range[1]:acf2spectrum.range[2]))},
           'm'={return(power.spect.m(signal=acf.section, Periods=powerspect.periods, norm.spect=TRUE, s=acf2spectrum.range[1]:acf2spectrum.range[2]))},
           'k'={return(power.spect.k(signal=acf.section, Periods=powerspect.periods, norm.spect=TRUE))},
           'mean'={return(power.spect.mean(signal=acf.section, Periods=powerspect.periods, norm.spect=TRUE, s=acf2spectrum.range[1]:acf2spectrum.range[2]))}
           )
  }
}
