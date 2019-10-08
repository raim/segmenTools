#' segmenTools : analysis of genome segmentations by segmenTier
#'@author Rainer Machne \email{raim@tbi.univie.ac.at}
#'@docType package
#'@name segmenTools
#'@section Dependencies: basic (\code{stats}, \code{graphics}, \code{grDevices}), clustering, \code{flowClust}, \code{flowMerge}
#'@importFrom utils write.table read.delim read.table
#'@importFrom graphics image axis par plot matplot points lines legend arrows strheight strwidth text abline hist spineplot polygon mtext layout
#'@importFrom grDevices png dev.off rainbow gray xy.coords rgb col2rgb  colorRampPalette densCols gray.colors
#'@importFrom stats mvfft ecdf loess predict qt quantile runmed sd var phyper heatmap rnorm kmeans approx fft smooth.spline median na.omit
##@bibliography /home/raim/ref/tata.bib
##@importFrom segmenTier clusterCor_c 
NULL # this just ends the global package documentation


### GENOME SEGMENT UTILS & ANALYSIS

### DATA STAT & TRANSFORMATION UTILS


#' asinh data transformation
#'
#' asinh trafo, an alternative to log transformation that has less
#' (compressing) effects on the extreme values (low and high values),
#' and naturally handles negative numbers and 0
#' @param x data to be transformed
#' @return log(x+sqrt(x^2+1))
#' @export
ash <- function(x) log(x+sqrt(x^2+1))
#' log trafo handling zeros by adding 1
#' @param x data to be transformed
#' @return log(x+1)
#' @export
log_1 <- function(x) log(x+1)

#' mean-0 normalization
#' @param x data to be transformed
#' @return t(apply(x,1,scale))
#' @export
meanzero <- function(x) t(apply(x,1,scale))

#' log2 ratio normalization
#'
#' @param x data to be transformed
#' @param na.rm remove NA values for mean calculation
#' @return log2(x/apply(x,1,mean))
#' @export
lg2r <- function(x,na.rm=TRUE) {
	y=as.matrix(log2(x/apply(x,1,mean,na.rm=na.rm)))
	if(sum(is.infinite(y))>0)
        warning("generate -Inf -> convert to NaN")
	y[!is.finite(y)]=NaN
	return(y)
	}

## box-cox trafo for negative values (Bickel and Doksum 1981)
## as used in flowclust
bc <- function(x,lambda) (sign(x)*abs(x)^lambda-1)/lambda
## amplitude box-cox trafo for complex polar coordinates 
bcdft <- function(x, lambda) {
    if ( class(x)=="matrix" )
        return(apply(x,2, bcdft, lambda))
    ## Box-Cox transform amplitude
    y <- bc(abs(x), lambda)
    ## amplitude scaling factor
    sf <- (y-min(y,na.rm=T))/abs(x)
    x*sf
}
## show effect of Box-Cox transformation
## on amplitude of complex DFT components
testbcdft <- function(tset, lambda=.5, cycle=3, col) {
    ## bsp. mit ts/csets 
    yft <- bcdft(tset$dft,lambda)

    ## colors
    if ( missing(col) )
        col <- rep("#00000077",nrow(tset$dat))
    else if ( class(col)=="clustering" )
        col <- clusterColors(col)

    #png("test_amplitude_boxcox.png",units="in",res=100,width=4.7,height=9)
    par(mfcol=c(2,1),mai=c(.5,.7,.01,.01), mgp=c(1,.25,0))
    plot(tset$dft[,cycle+1],cex=.5, col=col,
         xlab=paste0("Re_",cycle),ylab=paste0("Im_",cycle))
    abline(v=0,col=1,lwd=2)
    abline(h=0,col=1,lwd=2)
    plot(yft[,cycle+1],cex=.5, col=col,
         xlab=paste0("Re_",cycle),ylab=paste0("Im_",cycle))
    abline(v=0,col=1,lwd=2)
    abline(h=0,col=1,lwd=2)
    #dev.off()
}
## show effect of Box-Cox transformation in
## separate DFT components (dat is tset$dat)
testbc <- function(tset, lambda=.9, cycle=2, col) {

    dat <- tset$dat
    xy <- dat[,c(paste0("Re_",cycle),paste0("Im_",cycle))]
    
    ## colors
    if ( missing(col) )
        col <- rep("#00000077",nrow(tset$dat))
    else if ( class(col)=="clustering" )
        col <- clusterColors(col)

    par(mfcol=c(2,1),mai=c(.5,.7,.01,.01),mgp=c(1,.25,0))
    plot(xy,cex=.5,col=col)
    abline(v=0,col=1,lwd=2)
    abline(h=0,col=1,lwd=2)
    plot(bc(xy,lambda),cex=.5,col=col)
    abline(v=bc(0,lambda),col=2,lwd=2)
    abline(h=bc(0,lambda),col=2,lwd=2)
}

#' moving average using \code{\link[stats]{filter}}
#' @param x data vector along which a moving average will be calculated
#' @param n moving average window size
#' @param circular logical see help of function \code{\link[stats]{filter}}
#' @export
ma <- function(x, n=5, circular=FALSE) {
    stats::filter(x,rep(1/n,n), sides=2, circular=circular)
}

#' calculate 95% confidence intervals for the given
#' data vector using a t-distribution
#' @param x data vector for which the 95% confidence interval
#' will be calculated
#' @param na.rm remove NA values before calculation
#' @export
ci95 <- function(x,na.rm=FALSE) {
    if ( na.rm ) x <- x[!is.na(x)]
    n <- length(x)
    if ( n<2 ) return(NA)
    error <- qt(0.975, df=n-1) * sd(x)/sqrt(n)
    return(error)
}



### PLOT UTILS
#' Switch between plot devices
#' @param file.name file name without suffix (.png, etc)
#' @param type plot type: png, jpeg, eps, pdf, tiff or svg
#' @param width figure width in inches
#' @param height figure height in inches
#' @param res resolution in ppi (pixels per inch), only for 'png' and 'tiff'
#' @param bg background color
#' @export
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100,
                    bg="white") {
  file.name <- paste(file.name, type, sep=".")
  if ( type == "png" )
    grDevices::png(file.name, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    grDevices::postscript(file.name, width=height, height=width,paper="special",
                          horizontal = FALSE, onefile = FALSE)
  if ( type == "pdf" )
    grDevices::pdf(file.name, width=width, height=height, bg=bg)
  if ( type == "tiff" )
    grDevices::tiff(file.name, width=width, height=height, units="in", res=res)
  if ( type == "svg" )
    grDevices::svg(file.name, width=width, height=height)
  if ( type == "jpeg" )
    grDevices::jpeg(file.name, width=width, height=height, units="in", res=res)
}

##

#' 2D density heatmap plot
#'
#' Uses base R's \code{\link[grDevices:densCols]{densCols}} to
#' calculate local densities at each point in a scatterplot, and then
#' replaces them by a colored scheme, copied from Josh O'Brien posted
#' at
#' https://stackoverflow.com/questions/17093935/r-scatter-plot-symbol-color-represents-number-of-overlapping-points/17096661#17096661
#' @param x x-coordinates
#' @param y y-coordinates
#' @param pch \code{pch} argument to plot
#' @param nbin number of bins for both dimensions, can be a single
#'     number for both dimensions, or separate numbers, see argument
#'     \code{nbins} to function
#'     \code{\link[grDevices:densCols]{densCols}} and \code{gridsize}
#'     to \code{\link[KernSmooth:bkde2D]{bkde2D}}
#' @param colf color map function used to create a color gradient,
#'     eg. \code{\link[grDevices:colorRampPalette]{colorRampPalette}}
#'     or \code{viridis}
#' @param ... arguments to plot (hint:cex can be useful)
#' @return the plotted data.frame, including local densities
#' @export
dense2d <- function(x, y, pch=20, nbin=c(128,128), 
                    colf=grDevices::colorRampPalette(c("#000099","#00FEFF",
                                                      "#45FE4F", "#FCFF00",
                                                      "#FF9400", "#FF3100")),
                    ...) {
    
    df <- data.frame(x, y)
    
    ## TODO: nbin = 128 in densCols, vs. 256 colors in colorRampPalette
    ## KernSmooth::bkde2D(cbind(x, y), bandwidth=c(20,20), gridsize=c(128,128))
    ## Use densCols() output to get density at each point
    xcol <- grDevices::densCols(x, y, nbin=nbin,
                                colramp=grDevices::colorRampPalette(c("black",
                                                                      "white")))
    ## black - white, rgb components equal
    ## each corresponds to density, convert to number between 1 and 256
    df$dens <- col2rgb(xcol)[1,] + 1L
  
    ## Map densities to colors
    ncol <- 256 # TODO: fixed due to RGB range? 
    cols <- colf(ncol)
    df$col <- cols[df$dens]
    
    ## Plot it, reordering rows so that densest points are plotted on top
    plot(y~x, data=df[order(df$dens),], pch=pch, col=col, ...)
    
    ## normalize to 1, for legend?
    df$dens <- df$dens/ncol
    
    invisible(df)
}

#' Map a distribution to a color scheme
#'
#' Returns a list containing colors and breaks, to
#' be used in \code{image} plots, and the cut data and colors
#' for each datum. Note that the function is specifically used
#' for color scheme selection in the \code{genomeBrowser} project.
#' @param x data vector
#' @param mn minimal value to cut data at, if missing
#' \code{\link[stats:quantile]{quantile}} will
#' be used with parameter \code{q} or \code{q[1]}
#' @param mx maximal value to cut data at, if missing the quantile
#' \code{q[2]} or, if only one {q} is provided, \code{1-q[1]} will be used
#' @param q parameter \code{probs} for quantile selection of
#' mn/mx data cuts, here a single value or two values; see \code{mn}
#' and \code{mx} arguments
#' @param colf color map function used to create a color gradient,
#' eg. \code{\link[grDevices:gray.colors]{gray.colors}} or
#' \code{viridis}
#' @param n number of colors, first argument to \code{colf}
#' @param plot plot a legend, ie. a histogram of the full data and an
#' \code{image} of the color scheme below
#' @param xlab x-axis label to use in legend plot
#' @param heights relative heights of histogram and image plots
#' @param xlim optional limits to x-axis
#' @param mai \code{par("mai")} plot parameter for plot margins
#' @param ... further arguments to \code{colf}, e.g. \code{start} and
#' \code{end} in \code{\link[grDevices:gray.colors]{gray.colors}} 
#' @export
selectColors <- function(x, mn, mx, q=.1, colf=grDevices::gray.colors,  n=100,
                         plot=TRUE, xlab="score", xlim, heights=c(.75,.25),
                         mai=c(.75,.75,.1,.1), ...) {

    ## full data range
    rng <- range(x, na.rm=TRUE)
    brk <- seq(rng[1], rng[2], length.out=n+1)

    ## cut data
    if ( length(q)==1) q <- c(q, 1-q)
    if ( missing(mn) ) mn <- quantile(x, q=q[1], na.rm=TRUE) 
    if ( missing(mx) ) mx <- quantile(x, q=q[2], na.rm=TRUE) 
    x.cut <- x
    x.cut[which(x.cut>mx)] <- mx
    x.cut[which(x.cut<mn)] <- mn

    ## colors & breaks for cut data
    cols <- colf(n, ...)
    x.cols <- cols[1+(n-1)*(x.cut-min(x.cut))/(max(x.cut)-min(x.cut))]
    cbrk <- seq(mn, mx, length.out=n+1)

    if ( missing(xlim) )
        xlim <- rng
    
    ## plot legend
    if ( plot ) {

        lbrk <- seq(min(brk),max(brk),length.out=length(cols)+1)
        cdat <- lbrk
        cdat[cdat>mx] <- mx
        cdat[cdat<mn] <- mn
        
        layout(t(t(1:2)), heights=heights, widths=1)
        mai[1] <- mai[1]/5
        par(xaxs="i", yaxs="i", mai=mai)
        hist(x.cut,breaks=brk, border=2, col=2, xlim=xlim, axes=FALSE, main=NA)
        hist(x,breaks=brk, add=TRUE)
        abline(v=c(mn,mx), col=2, lty=2)
        axis(2)
        axis(1, labels=NA)
        mai[1] <- mai[1]*5
        par(xaxs="i", mai=mai)
        image_matrix(x=lbrk, z=t(cdat), breaks=cbrk,
                     col=cols, axis=1, ylab="color", xlab=xlab, xlim=xlim)
    }
    list(breaks=cbrk, col=cols, x.cut=x.cut, x.col=x.cols, xlim=xlim)
}

#' plot multiple cumulative distribution functions
#' 
#' plot multiple cumulative distribution functions of overlap statistics, as
#' provided by function \code{\link{getOverlapStats}}
#' @param x x-values for which cumulative distribution functions of overlap characteristics are plotted
#' @param CDF named list of cumulative distribution functions
#' @param type character indicating the type of the overlap CDF
#' @param range draw the range of values either as "polygon" or "lines"
#' @param col lines color vector
#' @param lty line type vector
#' @param h horizontal cut-off lines
#' @param v vertical cut-off lines
#' @param ylab y-axis label
#' @param ylim y-axis limits
#' @param ... further parameters to plot
#' @export
plot_cdfLst <- function(x=seq(0,2,.05), CDF, type="rcdf", range="polygon",
                        col, lty, h=c(.2,.8), v=c(0.8,1.2),
                        ylab="cum.dist.func.", ylim=c(0,1), ...) {


    ## group by colors and plot mean and ci95 as polygon!
    cls <- sort(unique(col))
    if ( length(cls)<length(col) ) {
        #plot_cdfGroups(x, )
        cdfmat <- matrix(NA, nrow=length(x), ncol=length(CDF))
        for ( i in 1:length(CDF) )
            cdfmat[,i] <- CDF[[i]][[type]](x)
        ## calculate meana and ci95
        cdfmn <- matrix(NA,nrow=length(x), ncol=length(cls))
        colnames(cdfmn) <- cls
        cdfci <- cdflo <- cdfhi <- cdfmn
        for ( cl in as.character(cls) ) {
            cdfmn[,cl] <- apply(cdfmat[,col==cl],1,mean)
            cdfci[,cl] <- apply(cdfmat[,col==cl],1,ci95)
            cdflo[,cl] <- apply(cdfmat[,col==cl],1,min)
            cdfhi[,cl] <- apply(cdfmat[,col==cl],1,max)
        }
        plot(x,rep(1,length(x)),col=NA,xlim=range(x),ylim=ylim,
             ylab=ylab,...)
        abline(v=v,lty=2)
        abline(h=h,lty=2)
        abline(h=0:1, lty=2, col="gray",lwd=.75)
        for ( cl in as.character(cls) ) {
            if ( range=="polygon" ) {
                px <-c(x,rev(x))
                ##py <- c(cdfmn[,cl]+cdfci[,cl],rev(cdfmn[,cl]-cdfci[,cl]))
                py <- c(cdfhi[,cl],rev(cdflo[,cl]))
                polygon(px, py, col=sub("FF$","77",cl),border=NA)
            } else if ( range=="lines" ) {
                lines(x,cdfhi[,cl], col=cl,lwd=1, lty=2)
                lines(x,cdflo[,cl], col=cl,lwd=1, lty=2)
            }
            lines(x,cdfmn[,cl], col=cl,lwd=3, lty=1)
        }
    } else {    
        plot(x, CDF[[1]][[type]](x),type="l",
             col=NA,main=NA, ylim=ylim, ylab=ylab, ...)
        abline(v=v,lty=2)
        abline(h=h,lty=2)
        abline(h=0:1, lty=2, col="gray",lwd=.75)
        for ( i in 1:length(CDF) ) 
          lines(x,CDF[[i]][[type]](x),col=col[i],lty=lty[i])
    }
}



### SEGMENT UTILS

#' which segment covers a position
#' @param pos chromosome position in continuous index
#' @param seg segment (genomic interval)
#' @export
whichSegment <- function(pos, seg) 
    which(seg[,"start"]<= pos & seg[,"end"]>=pos)

#' splits segmenTier classes into a table
#' 
#' specifically tailored to segment strings in segmenTier;
#' splits a list of strings by a separator and constructs
#' a table from the entries.
#' @param sgtypes a list of strings, that is converted to a table
#' of classes based on a string separator (\code{sep})
#' @param sep the separator for classes in the string
#' @param gsep a separator within classes
#' to be used as classification ID (column header of the produces
#' class table)
#' @export
getSegmentClassTable <- function(sgtypes, sep="_", gsep=":") {
    ## get classes
    cllst <- strsplit(sgtypes,"_")
    if ( length(unique(unlist(lapply(cllst, length))))>1 )
        cat(paste("WARNING: inconsistent segment classes\n"))
    ## get class ids
    clid <- unique(unlist(lapply(cllst, function(x) sub(":.*","",x))))

    ## fill data.frame
    cltab <- data.frame(matrix(NA, nrow=length(sgtypes), ncol=length(clid)))
    colnames(cltab) <- clid
    rownames(cltab) <- sgtypes
    for ( i in 1:length(cllst) ) {
        tmp <- strsplit(cllst[[i]],":")
        class <- unlist(lapply(tmp, function(x) x[2]))
        names(class) <- unlist(lapply(tmp, function(x) x[1]))
        cltab[i, names(class)] <- class
    }
    ## convert numeric to numeric
    for ( i in 1:ncol(cltab) ) {
        col <- cltab[,i]
        if ( suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
            cltab[,i] <- as.numeric(col)
        } else {
            cltab[,i] <- as.factor(col)
        }
    }
    ## rm non-varied
    rmc <- apply(cltab, 2, function(x) length(unique(x))==1)
    cltab <- cltab[,!rmc,drop=FALSE]
    
    cltab
}

#' splits segmenTier segment class strings into classes
#' 
#' specifically tailored to segment strings in segmenTier
#' @param sgtypes a list of strings, that is converted to a table
#' of classes based on a string separator (\code{sep})
#' @param sep the separator for classes in the string
#' @param gsep currently not used, a separator within classes
#' to be used as classification ID (column header of the produces
#' class table)
#' @export
getSegmentClasses <- function(sgtypes, sep="_", gsep=":") {
    sgclasses <- sort(unique(unlist(strsplit(sgtypes,"_"))))
    ## filter classes that are not varied
    rm <- NULL
    for ( class in sgclasses ) 
        if ( length(grep(class,sgtypes))==length(sgtypes) )
            rm <- c(rm, which(sgclasses==class))
    if ( !is.null(rm) )
        sgclasses <- sgclasses[-rm]
    sgclasses
}

### WRAPPERS of segmentOverlap for specific purposes
## calculate overlap between two sets of genome segments, e.g.
## test a segmentation of RNA-seq data vs. known features or transcripts
## NOTE: all chromosome coordinates must be mapped to a continuous index
## via coor2index
## NOTE: mod from $TATADIR/yeast/scripts/analyzeSeq2013_utils.R

## wrapper around \code{\link{segmentOverlap}}, used to 
## annotate the query set by a column in the target set
## @param query the query set of segments (genomic intervals)
## @param target the target set of segments (genomic intervals)
## @param col column names to copy from target to it's best match in query
annotateQuery <- function(query, target, col) {
    ## TODO: fix this, make real annotate query
}

## is query upstream of target?
## TODO:
## 0) coor2index w/o strand: ordered coors (start > end), remember strand
## 1) sort by start
## 2) for each start, consider all abs(end-start) <= maxD (cut at chrS)
## 3) calculate distances a=start-start, b=end-start, c=end-end, d=start-end
## NOTE: d = a + b + c [b is negative for overlapping]
## 4) classify into left, covers/equals, right
## 5) strands: translate to upstream, downstream, overlaps/antisense 
## 6) report characteristic distances (and jaccard for overlaps?)
orientation <- function(q, t, coors, frw.id=c(1,"+")) {
    chr <- coors[q,"chr"]
    if ( chr != coors[t,"chr"] ) return(NA)
    frw <- coors[q,"strand"] %in% frw.id
    str <- coors[q,"strand"] == coors[t,"strand"]
    ss <- coors[q,"start"] - coors[t,"start"]
    ee <- coors[q,"end"] - coors[t,"end"]
    se <- coors[q,"start"] - coors[t,"end"]
    es <- coors[q,"end"] - coors[t,"start"]

    ## multiply with ifelse(frw, 1, -1) ?

    ## if frw:
    ## q upstream of t: str & ss < 0 & es<0 -> - es
    ## q downstream t:  str & ss > 0 & se>0 -> - se 
    ## q/t convergent: !str & ee < 0 -> ee
    ## q/t divergent:  !str & ss > 0 -> ss
    
    if ( str ) {
        
    } else {
    }
}

## TODO: pairwise genome segment annotation
## TODO: convert query to target according to rules
## TODO: 20180402
## functions upstream/downstream/overlaps
## bandSparse(nque, ntar, ...?) matrices ss, ee, es, se, Jaccard
## -> rules to calculate bandSparse matrices:
## upstream, downstream, divergent, convergent, each with characteristic
## distance
##
## pairwise segment anntation, a wrapper around
## \code{\link{annotateTarget}} allowing for high-level
## inferences on relative positions of genome segments
## to each other
segmentPairs <- function(query, qcol="ID", chrS, distance, verb=1,
                         rules=c("divergent","convergent","antisense"),
                         strands=list(frw=c("1","+"),rev=c("-1","-"))) {

 }

#' annotate target segments by overlapping query segments
#' 
#' wrapper around \code{\link{segmentOverlap}}, used to 
#' annotate the target set by a column in the query set. See
#' Details for required input formats.
#'
#' Wrapper around  \code{\link{segmentOverlap}}, used to 
#' annotate the target set by a column in the query set. Note that
#' coordinates must already be `indexed', see \code{?\link{segmentOverlap}}
#' and \code{?\link{coor2index}}. The default options (\code{only.best=TRUE},
#' \code{collapse=TRUE} yield a result matrix with 1 row for each
#' \code{target} (\code{nrow(results)==nrow(target)}, annotated by
#' only the best ranking overlapping \code{query} segment(s), ie. the
#' the \code{query} segment(s) with the highest Jaccard index. If multiple
#' \code{query} segments have the same Jaccard index, they are collapsed
#' into ;-separated lists. If \code{collapse=FALSE}, each overlapping
#' \code{query} segment will have its own row, and the result matrix
#' is longer (\code{nrow(results)>=nrow(target)}). If \code{only.best=FALSE},
#' all overlapping \code{query} segments will be reported, and optionally
#' sorted by rank (\code{sort}).
#' Additionaly a minimal Jaccard filter can be applied
#' (option \code{minJaccard}), and only \code{query} segments with a
#' higher Jaccard index are reported.
#' @param query the query set of segments (genomic intervals)
#' @param target the target set of segments (genomic intervals)
#' @param qcol column names to copy from target to it's best match in query,
#' defaults to all columns
#' @param tcol column names of target to include in the returned table,
#' defaults to none
#' @param prefix column name prefix for the copied columns
#' @param details set to \code{TRUE} to add details of the used
#' match, i.e., union and intersect, and relative position of the query
#' to the matching targets
#' @param only.best only consider the top-ranking query hit
#' @param sort sort hits by rank
#' @param minJaccard minimal Jaccard index (intersect/union) threshold;
#' overlaps below this value will NOT be reported; NOTE that this will
#' set argument \code{details} to \code{TRUE}
#' @param collapse  if \code{TRUE} multiple query hits are collapsed
#' into a single row, with ;-separated strings in the respective fields.
#' @param msgfile file pointer for progress messages and warnings, defaults to
#' stdout, useful when using in context of command line pipes
#' @export
annotateTarget <- function(query, target, qcol=colnames(query), tcol,
                           prefix, details=FALSE, only.best=TRUE, sort=TRUE,
                           minJaccard, collapse=TRUE, msgfile=stdout()) {

    ## activate "details" for jaccard filter
    if ( !missing(minJaccard) )
        details <- TRUE

    ## force sorting with only.best option!
    if ( only.best ) sort <- TRUE
    
    ## TODO: use details flag to also bind details of overlap (left/right)
    #cltr <- annotateQuery(query, target, qcol)
    cltr <- segmentOverlap(query=query, target=target,
                           collapse=FALSE,
                           add.na=TRUE, sort=sort, untie=FALSE,
                           details=details, msgfile=msgfile)

    ## bind query column to overlap table
    cltr <- cbind(cltr, query[cltr[,"query"],qcol,drop=FALSE])

    ## REDUCE TO TOP-RANKING QUERY HIT
    best <- cltr
    if ( only.best )
      best <- cltr[cltr[,"qrank"]==1,,drop=FALSE]
    ## APPLY THRESHOLD
    if ( !missing(minJaccard) ) {
        rmJ <- which(best[,"intersect"]/best[,"union"] < minJaccard)
        best[rmJ,2:ncol(best)] <- NA
        ## TODO 20190124: avoids NA;NA in results!
    }
    
    ## collapse or remove duplicated with qrank==1
    ## TODO: is this necessary, or can we just use collapse=TRUE
    ## in call to segment overlap and strsplit required columns?
    if ( sum(duplicated(best[,"target"])) )
        cat(paste("handling", sum(duplicated(best[,"target"])), "duplicated:",
                  ifelse(collapse,"collapse","split"), "\n"),file=msgfile)
    dups <- rdups <- which(duplicated(best[,"target"]))
    if ( collapse ) {
        ## de-factor before pasting!! (20180423)
        for ( k in 1:ncol(best) )
            if ( is.factor(best[,k]) )
                best[,k] <- as.character(best[,k])
        while ( length(rdups)>0 ) {
            d <- rdups[1]
            id <- best[d,"target"]
            idx <- which(best[,"target"]==id) ## find all 
            ovl <- best[idx,]
            covl <- rep(NA, ncol(ovl))
            names(covl) <- colnames(ovl)
            for ( j in 1:ncol(ovl) )
                covl[j] <- paste(ovl[,j],collapse=";")
            covl["target"] <- ovl[1,"target"]
            covl["tlen"] <- ovl[1,"tlen"]
            if ( only.best ) 
                covl["qrank"] <- ovl[1,"qrank"]
            best[idx[1],] <- covl
            rdups <- rdups[!rdups%in%idx]
        }
    
        if ( length(dups)>0 )
            best <- best[-dups,]
    }
    
    ## e.g., details = c("qpos", "intersect", "union")
    if ( details )
        qcol <- c(qcol, "qpos", "qrank", "qlen", "intersect", "union")
    
    ## get and optionally rename requested columns
    ## TODO: why was match used? use dplyr::join !
    #idx <- match(best[,"target"],1:nrow(target))
    addcol <- best[,c("target",qcol), drop=FALSE]

    ## add prefix to column ids
    if ( !missing(prefix) )
        if ( prefix!="" )
            colnames(addcol) <- paste(prefix,"_",colnames(addcol),sep="")

    ## filter target columns
    if ( !missing(tcol) ) {
        target <- target[best[,"target"],tcol,drop=FALSE]
        colnames(target) <- tcol
        addcol <- cbind(target, addcol)
    }
    addcol
    ## add index to column ID if already present
    ## TODO: this is from clustering in segmenTier, adapt
    ## and make smarter: count existing indexes
    ##if ( any(duplicated(colnames(target))) ) {
    ##    sel <- paste(colnames(target),".1",sep="")
    ##    cnt <- 2
    ##    while( sum(duplicated(sel)) ) {
    ##        sel[duplicated(sel)] <- sub("\\..*",paste(".",cnt,sep=""),
    ##                                    sel[duplicated(sel)])
    ##       cnt <- cnt+1
    ##    }
    ##    colnames(target) <- sub("\\.1$","",sel)
    ##}
}

#' Collect statistics from from \code{\link{segmentOverlap}}
#' 
#' Collect statistics from from \code{\link{segmentOverlap}},
#' a nested list of overlap statistics in lists
#' and vectors for a given segment \code{type}
#' @param ovlStatLst list of overlap statistics, where individual
#' entries come from \code{\link{segmentOverlap}}
#' @param type name of the overlap statistics list first level
#' @export
collectOvlStats <- function(ovlStatLst, type) {

    lst <- ovlStatLst[[type]]
    ## rm empty
    len <- unlist(lapply(lst, function(x) length(x) ))
    if ( any(len==1) )
      cat(paste("removing",sum(len==1),"empty tests:",
                paste(names(lst)[len==1],collapse="; "),"\n"))
    lst <- lst[len>1]
    
    ## collect stats
    CDF <- lapply(lst, function(x) if ( length(x)>1 ) x$CDF)
    class(CDF) <- "cdfLst"
    DIST <- lapply(lst, function(x) if ( length(x)>1 ) x$DIST)
    height <- matrix(unlist(lapply(lst,
                                   function(x) x$height)),ncol=2,byrow=TRUE)
    jaccard <- unlist(lapply(lst,
                             function(x) x$jaccard))
    j.prcnt <- unlist(lapply(lst,
                             function(x) x$j.prcnt))
    j.cutoff <- unlist(lapply(lst,
                             function(x) x$j.cutoff))
    #j.prcnt<- covlStats$j.prcnt# percent of targets covered with J>threshold
    hitnum <- unlist(lapply(lst,
                                function(x) x$hitnum))
    numhit <- unlist(lapply(lst,
                            function(x) x$numhit))
    qnum <- unlist(lapply(lst,
                          function(x) x$qnum))
    tnum <- unique(unlist(lapply(lst,
                                 function(x) x$tnum)))
    nms <- names(CDF)
    names(jaccard) <- names(hitnum) <- names(numhit) <-
        names(qnum) <- names(CDF) <- names(DIST) <- NULL
    if ( length(tnum)>1 ) # temporary until sure that it works
        stop("ERROR,",type,"tnum should be unique for types")
    res <- list(CDF=CDF, DIST=DIST,
                jaccard=jaccard, j.prcnt=j.prcnt, j.cutoff=j.cutoff,
                height=height, hitnum=hitnum, numhit=numhit,
                qnum=qnum, tnum=tnum, nms=nms)
    class(res) <- "overlapStatLst"
    res
}

## TODO
plotOverlap <- function(ovlstats,type="rcdf",file.name) {
    if ( !missing(file.name) )
        png(file.name,width=400,height=400)
    ## TODO
   if ( !missing(file.name) ) dev.off()   
}


#' Jaccard vs. Ratio Segment Overlap Plots
#'
#' Generates a plot of Jaccard indices vs. query/target overlap Ratios,
#' optionally decorated with threshold lines and counts.
#' @param jaccard Jaccard Indices
#' @param ratio query/target length ratio
#' @param symm plot symmetric version, ie. for ratio>1 plot 1/ratio with
#' inverted y-axis
#' @param nbin \code{nbin} parameter for \code{\link{dense2d}}
#' @param minn minimal number of available points to use \code{\link{dense2d}},
#' a conventional scatter plot (\code{\link{points}}) will be used if less
#' than \code{minn} data are available
#' @param rlim ratio (y-)axis limis
#' @param jlim Jaccard (x-)axis limits
#' @param decorate decorate the plot with Jaccard threshold and quartal
#' numbers
#' @param j.thresh Jaccard threshold, line will be drawn at this point,
#' and counts reported for values below and above
#' @param yd fraction of plot height to use for quartal number plot
#' decoration
#' @param did optional data ID to be used in the plot legend
#' @param tot optional total number of targets to indicate non-overlapping
#' targets in the plot legend
#' @param ... further arguments to \code{\link{dense2d}} or \code{\link{points}}
#' @export
jrplot <- function(jaccard, ratio, symm=TRUE, nbin=512,  minn=10, 
                   rlim=c(0,2), jlim=c(0,1.05), decorate=TRUE,
                   j.thresh=0.5, yd=.1,did="", tot, ...) {

    ## legend text
    leg.txt <- ifelse(did=="","",paste0(did,":\n"))
    leg.txt <- paste0(leg.txt, sum(!is.na(jaccard)))
    if (!missing(tot) )
        if ( !is.na(tot) )
            leg.txt <- paste0(leg.txt,"/",tot)

    ## scale to symmetric plot!
    if ( symm )
        ratio[which(ratio>1)] <- 2 - 1/ratio[which(ratio>1)]
    
    plot(jaccard,ratio, ylim=rlim, xlim=jlim, col=NA,
         xlab=expression("Jaccard Index,"~J*"="*I/U),
         ylab=NA, axes=FALSE)


    ## draw min/max lines and threshold for J vs. R
    if ( symm )
        abline(a=2, b=-1, col="black", lwd=1, lty=2)
    else {
        j <- seq(0,2,.01)
        lines(j,1/j, col="black", lwd=1, lty=2)
    }
    abline(a=0,b=1, col="black", lwd=1, lty=2)
    abline(v=j.thresh, col="red", lwd=2, lty=2)

    ## axes
    ##abline(v=1, col="red", lwd=2, lty=2)
    abline(h=1)
    axis(1)
    ## construct axes for symmetric plot
    if ( symm ) {
        axis(2, labels=NA, col.ticks=NA) # empty line
        ## long tick at symmetry axis 1
        mgp <- par("mgp")
        mgp1 <- mgp
        mgp1[2] <- mgp1[2]*2
        axis(2, at=1, labels=1, cex=1.1, tcl=2*par("tcl"), mgp=mgp1)
        ## axis for Q/T
        qt <- pretty(seq(0,1,.1))
        axis(2, at=qt[-length(qt)])
        ## axis for T/Q
        xt <- pretty(seq(1,2,.1))
        xt <- xt[-1]
        axis(2, at=xt, labels=2-xt)
        mtext(expression(R*"="*Q/T), 2, mgp[1], adj=.25)
        mtext(expression(1/R*"="*T/Q), 2, mgp[1], adj=.75)
    } else {
        axis(2)
        mtext(expression("Length Ratio,"~R*"="*Q/T), 2, mgp[1])
    }
    #box()
    df <- NULL
    if ( length(na.omit(jaccard))<minn )
        points(jaccard,ratio, pch=20, col="#00000077", ...)
    else {
        par(new=TRUE)
        df <- dense2d(jaccard,ratio, ylim=rlim, xlim=jlim, nbin=nbin,
                      xlab=NA, ylab=NA, axes=FALSE, ...)
    }

    ## plot decoration
    if ( decorate ) {

        ## total number legend
        if ( FALSE ) # TODO: make optional or do externally?
            legend("top",leg.txt, bg="#FFFFFFAA",box.col=NA, cex=1.5)
        else
            mtext(sub("\\n"," ", leg.txt), 3, mgp[1]/2, cex=1.5)
        
        ## quartal numbers
        
        ## above Jaccard threshold
        good <- ratio >= 1 & jaccard>j.thresh
        bad <-  ratio <  1 & jaccard>j.thresh
        
        ylim <- par("usr")[3:4]
        ydff  <- yd*diff(ylim)
        ylow  <- ylim[1]+ydff
        yhigh <- ylim[2]-ydff
        
        text(j.thresh, yhigh-ydff,
             bquote(Q>T*":"~.(sum(good,na.rm=TRUE))),pos=4)
        text(j.thresh, ylow +ydff,
             bquote(Q<T*":"~.(sum(bad, na.rm=TRUE))),pos=4)
        
        ## below jaccard threshold
        long  <- ratio>=1 & jaccard<=j.thresh
        short <- ratio< 1 & jaccard<=j.thresh
        shorts <- sum(short, na.rm=TRUE)
        longs  <- sum(long,  na.rm=TRUE)
        
        text(j.thresh, yhigh, bquote(Q>T*":"~.(longs)),pos=2)
        text(j.thresh, ylow,  bquote(Q<T*":"~.(shorts)),pos=2)
        
        ## this was only useful for !symm
        ##if ( any(ratio>rlim[2],na.rm=TRUE) )
        ##    text(1/rlim[2],rlim[2],paste0(sum(ratio>rlim[2],na.rm=TRUE),
        ##                                  ": R>", rlim[2]),
        ##         pos=4, cex=.7)
    }
    
    ## pass on dense2d return for potential density legend
    invisible(df)
}

### MAIN ALGORITHM - FINDS OVERLAPS BETWEEN TWO SETS
### CHROMOSOMAL SEGMENTS

#' Overlaps between two sets of chromosome segments
#' 
#' Simple sweeping algorithm to find overlapping intervals. Both
#' intervals must be in continuous index (coor2index) and start<end,
#' with strand information implied in the continuous index position.
#' see \code{bedtools} \code{intersects}, the Binary Interval
#' Search (BITS) algorithm, and the \code{NClist} algorithm (implemented
#' in packages IRanges and GenomicRanges) for similar tools.
#' It loops over targets and collects all queries that match, and
#' optionally adds NA lines for non-matched targets.
#'
#' If \code{details==TRUE}, a column \code{qpos} will indicate the
#' relative position of the query to the target (e.g. `inside' means
#' that the query is `inside' the target, `left' means at lower
#' coordinates then the target). The positions `right' and `left'
#' can later can be converted to `upstream/downstream/etc.' by
#' \code{\link{index2coor}} which re-introduces strand information.
#' To analyze relative positions between features on distinct
#' strands, the target or query has to be mapped to the opposite
#' strand by \code{\link{switchStrand}} before passing it to
#' \code{segmentOverlap}.
#' @param query the query set of segments (genomic intervals)
#' @param target the target set of segments (genomic intervals)
#' @param details add details on the relative positions of the query
#' wrt the target (covers,left,inside,right); note, that  keywords 'left'
#' and 'right' will be converted to 'upstream' or 'downstream' by
#' index2coor, accounting for strand info.
#' @param distance not yet implemented: maximal distance for neighbor count
#' @param add.na add an empty overlap line for targets without overlaps
#' @param untie if several queries have equal top rank of 1, the rank
#' can be  replaced by simple order, such that only the first query
#' will have rank 1
#' @param collapse if \code{TRUE} multiple query hits are collapsed
#' into a single row, with ;-separated strings in the respective fields.
#' If both \code{add.na} and \code{collapse} are set to TRUE the resulting
#' overlap matrix will be of the same dimension as the input target, i.e.,
#' it will contain one row for each target.
#' @param sort sort query hits by their \code{rank}
#' @param msgfile file pointer for progress messages and warnings, defaults to
#' stdout, useful when using in context of command line pipes
#' @export
segmentOverlap <- function(query, target, details=FALSE, distance,
                           add.na=FALSE, untie=FALSE, collapse=FALSE,
                           sort=FALSE, msgfile=stdout()) {

    ## get target and query ID - only required for messages
    if ( "ID" %in% colnames(target) ) {
        tids <- as.character(target[,"ID"])
    } else if ( !is.null(rownames(target)) ){
        tids <- rownames(target)
    } else {
        tids <- 1:nrow(target)
    }
    if ( "ID" %in% colnames(query) ) {
        qids <- as.character(query[,"ID"])
    } else if ( !is.null(rownames(query)) ){
        qids <- rownames(query)
    } else {
        qids <- 1:nrow(query)
    }

    overlaps <- NULL
    for ( k in 1:nrow(target) ) {

        ## overlapping 
        idx <- which(query[,"end"] >= target[k,"start"])
        idx <- idx[which(query[idx,"start"] <= target[k,"end"])]

        ## non-overlapping neighbors
        ## TODO: implement and use below!?
        if ( !missing(distance) ) {
            #nidx <- which(target[k,"start"] - query[,"end"] <= distance
            ## left neighbor
            ##nleft <- query[,"end"] <= target[k,"start"]
            dst <- target[k,"start"] - query[,"end"]
            nleft <- dst <= distance & dst > 0
                    
            ## right neighbor
            ##nrigt <- query[,"start"] >= target[k,"end"]
            dst <- query[,"start"] - target[k,"end"] 
            nrigt <- dst >= distance & dst > 0
        }
      
        ## detailed overlap structure
        ## TODO: rm redundant code above
        ## TODO: allow upstream and downstream ranges!
        if ( details ) {
            ## 1: query ends right of target start (ALL overlaps)
            RhLt <- query[,"end"]   >= target[k,"start"]
            ## 2: query starts left of target end (ALL overlaps)
            LlRt <- query[,"start"] <= target[k,"end"]
            ## 3: query starts right of target start (INSIDE or RIGHT)
            LhLt <- query[,"start"] >= target[k,"start"]
            ## 4: query ends left of target end (INSIDE or LEFT)
            RlRt <- query[,"end"]   <= target[k,"end"]
            ## 5: query ends right of target end (INSIDE, COVERS or RIGHT)
            ##    ~ !RlRt, but incl. end
            RhRt <- query[,"end"]   >= target[k,"end"]  
            ## 6: query starts left of target start (INSIDE, COVERS or LEFT)
            ##    ~ !LhLt, but incl. start
            LlLt <- query[,"start"] <= target[k,"start"]
            ## TODO: divergent/convergent/tandem neighbors
            ##    -> "neighbors", can later (in index2coor) be classified into
            ##    convergent/divergent (different strand in --antisense search,
            ##    and tandem (same strand, in --sense search)
            
            ## RULES
            ## combinations:
            overlp <- LlRt & RhLt # general overlap: should contain all below
            inside <- LhLt & RlRt # is inside target
            covers <- LlLt & RhRt # covers target
            inside[inside & covers] <- FALSE ## EXACT HIT
            left   <- !LhLt & !RhRt & RhLt # left overlapping
            right  <- !RlRt & !LlLt & LlRt # right overlapping

            ## debugging: these should NOT happen
            if ( sum(overlp,na.rm=T)!=length(idx) )
                stop(k,"ERROR, overlap vs. idx")
            if ( sum(overlp,na.rm=T) !=
                sum((covers|inside)|(left|right),na.rm=T) )
                stop(k,"ERROR, overlap vs. detail")
            if ( any(rowSums(cbind(inside,covers,right,left),na.rm=T)>1) )
                stop(k,"ERROR, too many hits")
        }

        ## add empty target if requested
        if ( length(idx)==0 ) {
            if ( add.na ) {
                result<-c(target=k,query=NA,intersect=0, union=0, qlen=0,tlen=0)
                if ( details ) result <- c(result,qrank=1,qpos=NA)
                else result <- c(result,qrank=1)
                overlaps <- rbind(overlaps,result)
            }
            next
        }
        ## calculate overlap %
        tgrng <- target[k,"start"]:target[k,"end"]
        tlen <- length(tgrng)
        ovl <- NULL
        tps <- NULL
        for ( id in idx ) {

            ## query coors
            qurng <- query[id,"start"]:query[id,"end"]
            qlen <- length(qurng)

            ## INTERSECT/UNION
            intersect <- length(intersect(tgrng,qurng))
            union <- length(union(qurng,tgrng))

            ## store overlaps of target k with query
            ## values to calculate ratios in \code{getOverlapStats(ovl)}:
            ## Jaccard measure: J(query,target) = intersect/union
            ## ratio:           ratio = qlen/tlen
            ## target coverage: tgcov=intersect/tlen
            ## query coverage:  qucov=intersect/qlen
            result <- c(target=k, query=id, 
                        intersect=intersect, union=union, 
                        qlen=qlen, tlen=tlen)
            ## relative query position wrt target
            if ( details ) {
                qpos <- NULL #factor(levels=c("left","covers","inside","right"))
                ## TODO: number positions from left to right
                if ( id %in% which(left) )   qpos <- c(qpos,"left")
                if ( id %in% which(covers) ) qpos <- c(qpos,"covers")
                if ( id %in% which(inside) ) qpos <- c(qpos,"inside")
                if ( id %in% which(right) )  qpos <- c(qpos,"right")

                ## debugging: these should NOT happen
                if ( is.null(qpos) )
                    stop(k, "ERROR: must have one qpos")
                if ( length(qpos)>1 )
                    stop(k, "ERROR: can only have one qpos")
                tps <- c(tps,qpos)
            }
            ovl <- rbind(ovl, result)
        }
        ## RANK: 1 for best match!
        ## rank by Jaccard=intersect/union
        jc <- ovl[,"intersect"]/ovl[,"union"]
        qrank <- rank( -unlist(jc), ties.method="min") 
        ## attempt to untangle ties by second criterium: query coverage
        ##if ( sum(duplicates(qrank))>0 ) {
        if ( sum(qrank==1)>1 ) {
            cat(paste("WARNING:",sum(qrank==1),
                      "with qrank 1 (max. jaccard=intersect/union) for target",
                      tids[k],":", paste(qids[idx[qrank==1]],collapse=";")),
                file=msgfile)
            if ( untie ) {
                cat(paste(" - un-tieing by order"),file=msgfile)
                ## NOTE: assignment of first hit to rank 1
                qrank[order(qrank)] <- 1:length(qrank)
                cat(paste(":",paste(qids[idx[qrank==1]],collapse=";")),
                    file=msgfile)
            }
            cat(paste(".\n"),file=msgfile)
        }
        
        rownames(ovl) <- NULL
        if ( details )  ovl <- data.frame(ovl,qrank=qrank,qpos=tps,
                                          stringsAsFactors=FALSE)
        else ovl <- data.frame(ovl,qrank=qrank)

        ## sort by qrank, top-hit (qrank==1) first
        ## TODO: test this
        if ( sort ) 
            ovl <- ovl[order(ovl[,"qrank"]),]
        
        ## collapse all query hits into ;-separated lists
        if ( collapse ) {
            covl <- rep(NA, ncol(ovl))
            names(covl) <- colnames(ovl)
            for ( j in 1:ncol(ovl) )
                covl[j] <- paste(ovl[,j],collapse=";")
            covl["target"] <- k
            covl["tlen"] <- tlen
            ovl <- covl
        }
        
        overlaps <- rbind(overlaps, ovl)
    }
    rownames(overlaps) <- NULL

    ## TODO - add target rank
    
    data.frame(overlaps, stringsAsFactors=FALSE)
}

#' Statistics of overlaps between two segment sets. 
#' @param ovl overlap table from \code{\link{segmentOverlap}}
#' @param ovlth threshold fraction of overlap of
#' target and query to be counted as 'good' hit
#' @param minj the minimal Jaccard index above which the fraction
#' is reported as \code{JaccPrcnt}
#' @param minf fraction cutoff, the Jaccard index reported for this
#' fraction of test sets is reported as \code{JaccCutoff}
#' @param hrng lower and upper thresholds of the ratio CDF (query length/target
#' length); the fraction of 'best' hits within this range will be reported
#' as 'height'
#' @param tnum number of targets in the original call, required for total CDF
#' ('acdf')
#' @param qnum number queries in the original call, just passed to results
#' @param tid ID for the target set, just passed on to results and used in plot
#' @param qid ID for the query set, just passed on to results and used in plot
#' @export
getOverlapStats <- function(ovl, ovlth=.8, minj=0.8, minf=0.2, hrng=c(.8,1.2), tnum=NA, qnum=NA, qid=NA, tid=NA) {

    ## returning NULL if no query had been found
    if ( !any(!is.na(ovl[,"query"])) )
      return(NULL)
    
    ## input ovl is from \code{segmentOverlap(ovl)}
    ## Calculate ratios:
    ## Jaccard measure: J(query,target) = intersect/union
    ## ratio:           ratio = qlen/tlen
    ## target coverage: tgcov=intersect/tlen
    ## query coverage:  qucov=intersect/qlen
    jaccard <- ovl[,"intersect"]/ovl[,"union"] 
    ratio <- ovl[,"qlen"]/ovl[,"tlen"]
    tgcov <- ovl[,"intersect"]/ovl[,"tlen"] # handle NA correctly?
    qucov <- ovl[,"intersect"]/ovl[,"qlen"] # -"-


    ## ALL HITS - total target coverage
    ## CDF of all hits, i.e. summing target coverages from multiple segments
    ## ?? todo: get total length and calculate ratio
    ## TODO: required?
    alltrg <- NULL
    if ( !is.na(tnum) ) {
        ## calculate total coverage of each target
        alltrg <- unlist(lapply(1:tnum,function(x)
            sum(tgcov[ovl[,"target"]==x])))
        alltrg[is.na(alltrg)] <- 0
        acdf <- ecdf(alltrg)
    }

    ## BEST HITS: rank==1 - hit with maximal coverage of target 
    ## TODO: rm duplicates properly, i.e. also multiple qrank==1!
    rnk <- ovl[,"qrank"]

    ## BEST HIT - total Jaccard
    J <- sum(ovl[rnk==1,"intersect"])/sum(ovl[rnk==1,"union"])

    ## total Jaccard
    ## todo: total Jaccard requires to account for unmachted qlen/tlen
    ##Jtot <- sum(ovl[,"intersect"])/sum(ovl[,"union"])
    
    ## 'GOOD' HIT: minimal mutual coverages > threshold
    ##             only used to calculate sum below
    sgt <- tgcov[rnk==1] >= ovlth & qucov[rnk==1]>= ovlth


    ## best hit target coverage distributions
    ## TODO: required?
    tgcovdist <- hist(tgcov[rnk==1],breaks=0:20/20,plot=FALSE)
    
    ## cumulative dist function of target coverage

    ## CDF of best hit target coverage
    ## TODO: required?
    tcdf <- ecdf(tgcov[rnk==1]) 

    ## relative best hit segment/test length ratio
    ## CDF of relative best hit ratio, i.e. w/o non-covered targets
    rat <- ratio[rnk==1]
    rrcdf <- ecdf(rat)   ## NA not counted!


   
    ## absolute best hit segment/test length ratio
    ## CDF of absolute best hit ratio, i.e., including non-covered targets
    ## TODO: align this with tnum
    if ( !is.na(tnum) ) {
        rat[is.na(rat)] <-0  # NOTE: difference only for add.na=TRUE
        if ( length(rat)<tnum ) # todo: why?
            rat <- c(rat,rep(0,tnum-length(rat) ))
        rcdf <- ecdf(rat)   
    }

    ## HEIGHTS: target recovery fraction within threshold - get from CDF
    height <- c(rcdf(hrng[1]),rcdf(hrng[2]))
    
    ## relative jaccard (intersect/union)
    ## CDF of relative jaccard ratio, i.e. w/o non-covered targets
    jac <- jaccard[rnk==1]
    rjcdf <- ecdf(jac)   ## NA not counted!

    ## absolute jaccard segment/test length ratio
    ## CDF of absolute jaccard ratio, i.e., including non-covered targets
    ## TODO: align this with tnum
    if ( !is.na(tnum) ) {
        jac[is.na(jac)] <-0  # NOTE: difference only for add.na=TRUE
        if ( length(jac)<tnum ) # todo: why?
            jac <- c(jac,rep(0,tnum-length(jac) ))
        jcdf <- ecdf(jac)
        j.prcnt <- 1-jcdf(minj) # fraction above jaccard cutoff
        j.cutoff <- quantile(jac,minf) # jaccard index at CDF cutoff
    }
   
    ## SINGLE OPTIMIZATION MEASURES
    ## distrubtion of queries per target
    tab <- table(ovl[,"target"])
    brks <- max(100, max(tab))
    tgnumdist <- hist(table(ovl[,"target"]), breaks=brks, plot=FALSE)
    ## TODO is mean a good measure?
    numhit <- mean(tab)            ## MINIMIZE: ~ too much fragmentation 
    hitnum <- sum(sgt, na.rm=TRUE) ## MAXIMIZE: ~ # of 'good' hits!

    ## LIST OF CDFs
    CDF <- list(rcdf=rcdf,rrcdf=rrcdf,tcdf=tcdf,acdf=acdf,
                  jcdf=jcdf,rjcdf=rjcdf,
                  qid=qid, tid=tid) # required when results are resorted
    ## LIST OF DISTRIBUTIONS
    DIST <- list(tgcovdist=tgcovdist, # best hit target coverage distribution
                 tgnumdist=tgnumdist) # num. of hits per target distribution
                
    ## TODO: other values? 
    ## TODO: which of these is actually useful except rcdf and rrcdf?
    res <- list(CDF=CDF, # list of CDFs of different measures
                DIST=DIST, # list of distributions
                height=height, # fraction of best hit ratios within hrng
                jaccard=J,     # MAX: jaccard measure of best hits
                j.prcnt=j.prcnt, # % of best hits with j > minj
                j.cutoff=j.cutoff, # jaccard index at CDF(minf) 
                hitnum=hitnum, # MAX: num. of hits with minimal mutual coverage
                numhit=numhit, # MIN: avg. num. of hits per target
                qnum=qnum, tnum=tnum,# number of queries/targets 
                qid=qid, tid=tid)
    class(res) <- "overlapStats"
    return(res)
}


#' Jaccard-Index overlap test for classes of segments (genomic intervals)
#'
#' calculates the Jaccard index, including a simple permutation test, between
#' different classes in a query and a target set of segments
#' (genomic intervals), where coordinates have been converted to a
#' continuous index over all chromosomes with \code{\link{coor2index}}.
#' Note, that this ignores chromosome borders!
#' @details Reports the Jaccard index (\code{J=intersect/union)}) between
#' two distinct sets of segments (matrix "jaccard" in the results object),
#' and the relative intersect lengths, i.e., intersect divided by
#' the total target length (matrix "intersect.target") or the total
#' query length (matrix "intersect.query").
#' If argument \code{perm>0}, a simple permutation is performed,
#' sampling randomly from
#' all inter-segment distances and segment lengths, and ignoring optional
#' query sub-classifications (argument \code{qclass}). Note, that chromosome
#' ends are ignored. The total length of the query range (genome length, for
#' both strands, if both are used) can be passed in argument \code{total},
#' and if missing the start of the first segment is also used as the distance
#' of the final segment to the query range end.
#' The results of the permutation test (argument \code{perm>0}) can be
#' plotted directly with \code{\link{plotOverlaps}}, where option \code{text}
#' allows to plot either the jaccard index or the relative intersect sizes.
#' @param query query set of segments
#' @param target target set of segments
#' @param qclass column name which holds a sub-classification (clustering) of
#' the query segments, omit or pass empty string ("") to use all
#' @param tclass column name which holds a sub-classification (clustering) of
#' the target segments, omit or pass empty string ("") to use all
#' @param total total length of the query range (genome length), if missing
#' the start of the first segment is also used as end
#' @param perm number of permutations to perform
#' @param verb integer level of verbosity, 0: no messages, 1: show messages
#' @export
segmentJaccard <- function(query, target, qclass, tclass, total, perm=0, verb=1) {
    if ( missing(qclass) ) qclass <- ""
    if ( missing(tclass) ) tclass <- ""

    ## query classes
    if ( qclass=="" ) {
        qcls <- as.factor(rep("query", nrow(query)))
        ## TODO: use this to fix perm use with qclass==""
        ## see below in permutation
        #query <- cbind(query,TMPCLASS=qcls)
        #qclass <- "TMPCLASS"
    } else {
        qcls <- as.factor(query[,qclass])
    }
    qcls.srt <- sort(unique(qcls))
    qN <- length(qcls.srt)
    
    ## target classes
    if ( tclass=="" ) {
        tcls <- as.factor(rep("target", nrow(target)))
    } else {
        tcls <- as.factor(target[,tclass])
    }
    tcls.srt <- sort(unique(tcls))
    tN <- length(tcls.srt)
    
    ## get full ranges for all query classes
    qcls.rng <- rep(list(NA), qN)
    names(qcls.rng) <- qcls.srt
    for ( i in 1:qN ) {
        cl <- qcls.srt[i]
        idx <- qcls==cl
        rng <- apply(as.matrix(query[idx,c("start","end")]), 1,
                     function(x) x["start"]:x["end"])
        names(rng) <- NULL
        qcls.rng[[cl]] <- unlist(rng)
    }
    
    ## get full ranges for all target classes
    tcls.rng <- rep(list(NA), tN)
    names(tcls.rng) <- tcls.srt
    for ( i in 1:tN ) {
        cl <- tcls.srt[i]
        idx <- tcls==cl
        rng <- apply(as.matrix(target[idx,c("start","end")]), 1,
                     function(x) x["start"]:x["end"])
        names(rng) <- NULL
        tcls.rng[[cl]] <- unlist(rng)
    }

    ## get intersect/union of all query:target class pairs
    J.real <- matrix(NA, nrow=qN, ncol=tN)
    colnames(J.real) <- tcls.srt
    rownames(J.real) <- qcls.srt
    ## additional data
    I.target <- I.query <- J.real # intersect/target
    for ( i in 1:qN ) {
        for ( j in 1:tN ) { 
            is <- length(intersect(qcls.rng[[i]],tcls.rng[[j]]))
            un <- length(union(qcls.rng[[i]],tcls.rng[[j]]))
            J.real[i,j] <- is/un
            I.target[i,j] <- is/length(tcls.rng[[j]]) 
            I.query[i,j] <- is/length(qcls.rng[[i]]) 
        }
    }
    #J.real

   
    if ( perm>0 ) {

        ## randomize queries

        J.pval <- J.real
        J.pval[] <- 0

        ## sort query
        query <- query[order(query$start),]
        for ( i in 1:perm ) {

            if ( verb>0 )
                cat(paste(i/perm," "))

            rquery <- randomSegments(query, qclass=qclass, total=total)

            ## TODO: test cluster length distribution?
            
            J.rnd <- segmentJaccard(rquery, target, qclass="type", tclass,
                                    perm=0)
            J.pval <- J.pval + as.numeric(J.rnd$jaccard >= J.real)
        }
        cat(paste("\n"))

        ## p-value
        J.pval <- J.pval/perm
    }

    ## results
    ovl <- list()
    ovl$jaccard <- J.real
    ovl$intersect.target <- I.target
    ovl$intersect.query <- I.query
    if ( perm>0 ) 
        ovl$p.value <- J.pval
    
    return(ovl)
}

#' randomize locations of input segments
#' 
#' randomizes the locations of input segments, while maintaining
#' the length distributions, optionally different distributions
#' of different segment classes. Coordinates have been converted to a
#' continuous index over all chromosomes with \code{\link{coor2index}}.
#' Note, that this ignores chromosome borders!
#' TODO: allow to pass chromosome length index chrS to account
#' for chromosome borders.
#'
#' The function is also used for permutation tests in
#' \code{\link{segmentJaccard}}.
#' @param query query set of segments to be randomized
#' @param qclass column name which holds a sub-classification (clustering) of
#' the query segments, omit or pass empty string ("") to use all
#' @param total total length of the query range (genome length), if missing
#' the start of the first segment is also used as end.
#' @export
randomSegments <- function(query, qclass, total) {
    
    if ( missing(qclass) ) qclass <- ""

    ## ORDER!
    query <- query[order(query$start),]

    ## query classes
    if ( qclass=="" ) {
        qcls <- as.factor(rep("query", nrow(query)))
        query <- cbind(query,TMPCLASS=qcls)
        qclass <- "TMPCLASS"
    } 
    qcls <- as.factor(query[,qclass])

    ## segment lengths
    sglen <- apply(query[,c("start","end")], 1, diff)
    sglen <- sglen+1
    sgcls <- as.factor(query[,qclass])
    
    ## inter-segment lengths
    ## TODO: handle nrow==1
    qnum <- nrow(query)
    islen <- apply(cbind(query[1:(qnum-1),"end"],
                         query[2:qnum,"start"]), 1, diff)
    ## NOTE: the start of the first real query segment
    ## is also used as the distance of the final segment
    ## unless a total length is explicitly provided as argument
    ## len; total =  2*sum(chrL)
    if ( !missing(total) ) end <- total - query[qnum,"end"]+1
    else {
        end <- query[1,"start"]
        total <- max(query[,c("start","end")])+end -1
    }
    islen <- c(query[1,"start"], islen, end)
    islen <- islen-1
    
    ## debug check whether total lengt is reproduced
    ## TODO: allow this check, but requires chrL to be known!
    if ( sum(islen)+sum(sglen) != total )
        stop()
    if ( FALSE )
        if ( any(islen<0) )
            cat(paste("NOTE: overlapping segments in randomized query set!\n"))
    
    ## RANDOMIZATION
    ## sample segment lengths
    ridx <- sample(1:qnum)
    rcls <- sgcls[ridx] # store cluster to keep cluster length dist!
    rsglen <- sglen[ridx]
    ## sample intersegment lengths
    rislen <- sample(islen)

    ## TODO: avoid start==0 which can happen for directly adjacent segments!
    
    ## construct randomized segmentation
    tot <- length(rsglen)+length(rislen)
    cumlen <- rep(NA,tot)
    cumlen[seq(1,tot,2)] <- rislen
    cumlen[seq(2,tot,2)] <- rsglen
    cumlen <- cumsum(cumlen)
    
    rquery <- data.frame(start=cumlen[seq(1,tot-1,2)],
                         end=cumlen[seq(2,tot-1,2)],
                         type=rcls)
    rquery
}

## FILL UP GENOME ADD MISSING SEGMENTS
#' Get inter-segments
#'
#' finds all ranges not present in the input segments, e.g. non-coding
#' regions in a list of coding regions
#' @param seg input segments with indexed coordinates (see
#' \code{\link{coor2index}})
#' @param chrS a chromosome index, indicating at wich positions
#' chromosomes start; this is required for handling chromosome ends
#' and forward and reverse strand values
#' @param indexed boolean value indicating whether the coordinates
#' in argument \code{seg} have already been indexed; if not
#' \code{\link{coor2index}}) will be applied and the returned inter-segments
#' will be mapped back to chromosome coordinates by \code{\link{index2coor}}
#' @param expand if \code{TRUE} the returned inter-segments will contain
#' the same columns as the input; otherwise only coordinates are returned
#'@export
fillGenome <- function(seg, chrS, indexed=TRUE, expand=TRUE) {

    ## find all non-covered segments

    ## apply coor2index
    if ( !indexed )
      seg <- coor2index(seg,chrS)
    ## expand segments to full range
    rng <- unique(unlist(apply(seg, 1, function(x) x["start"]:x["end"] )))
    ## find all not present
    all <- 1:(2*max(chrS))
    mss <- all[!all%in%rng]
    ## collapse adjacent to segments
    dff <- which(diff(mss)>1) 
    starts <- c(mss[1],mss[dff+1])
    ends <- c(mss[dff],mss[length(mss)])
    iseg <- cbind(start=starts,end=ends)
    ## cut at chromosome ends
    iseg <- splitsegs(iseg,chrS)

    ## map back to chromosome coordinates
    if ( !indexed )
      iseg <- index2coor(iseg,chrS)
    else iseg <- cbind(chr=1,iseg,strand=idx2str(iseg,chrS)[,1])
    
    ## expand to same matrix as input
    if ( expand ) {
        imat <- as.data.frame(matrix(NA, ncol=ncol(seg), nrow=nrow(iseg)))
        colnames(imat) <- colnames(seg)
        imat[,colnames(iseg)] <- iseg
        iseg <- imat
    }
    iseg
}

### SEGMENT READ STATISTICS

#' Calculates phase distributions
#' 
#' calculates the circular mean and R statistics, copied from package
#' \code{CircStats} and after
#' http://en.wikipedia.org/wiki/Directional_statistics#Measures_of_location_and_spread and
#' https://en.wikipedia.org/wiki/Mean_of_circular_quantities
#' TODO: recognize bimodal?
#' @param phs vector of phases in degrees
#' @param w optional weights for weighted mean
#' @param degrees phases are in degrees (0-360), set FALSE for radians
#' @export
phaseDist <- function(phs, w, degrees=TRUE) {
    
    avg <- c(mean=NA, r=NA, na=sum(is.na(phs))/length(phs))

    ## removing NA
    if ( !missing(w) ) w <- w[!is.na(phs)]
    if ( degrees ) 
        phs <- phs[!is.na(phs)] * pi/180 # convert degree to radian
    
    if ( length(phs)>0 ) {
        if ( !missing(w) ) {
            c <- sum(cos(phs)*w)
            s <- sum(sin(phs)*w)
            n <- sum(w)
        }else {
            c <- sum(cos(phs))
            s <- sum(sin(phs))
            n <- length(phs)
        }
        r <- sqrt(c^2 + s^2)/n
        mean <- atan2(s, c)
        if ( degrees ) {
            mean <- mean*180/pi # convert back to degree
            mean <- ifelse(mean<=  0,mean + 360, mean)
        }
        avg[1:2] <- c(mean=mean, r=r)
    }
    if ( !missing(w) ) names(avg) <- paste("w",names(avg),sep=".")
    avg
}
#' cluster p-value summary
#' 
#' calculates the mean of p-values and number of significant
#' p-values (<\code{threshold}) and number of NA values
#' @param pvs vector of p-values
#' @param threshold p-value threshold to be counted as signficant
#' @export
pvalDist <- function(pvs,threshold=0.01) {
  pvs[pvs==1] <- NA # convert back to NA
  avg <- c(p.mean=mean(pvs,na.rm=TRUE),
           p.signif=sum(pvs<threshold,na.rm=TRUE)/length(pvs),
           p.na=sum(is.na(pvs))/length(pvs))
  avg
}

#' distribution of read-counts, i.e., the mean, var, min&max and the number
#' of NA values
#' @param rds vector of read-counts, e.g., within a genomic interval (segment)
#' @export
readDist <- function(rds) {
    avg <- c(r.tot=sum(rds,na.rm=TRUE),
             r.len=length(rds),
             r.nrm=sum(rds,na.rm=TRUE)/length(rds),
             r.mean=mean(rds,na.rm=TRUE),
             r.var=var(rds,na.rm=TRUE),
             r.min=min(rds,na.rm=TRUE),
             r.max=max(rds,na.rm=TRUE),
             r.0=sum(rds==0,na.rm=TRUE)/length(rds),
             r.na=sum(is.na(rds))/length(rds))
  avg
}

#' Average read-counts of segments
#' 
#' Average read-count time-series of a genomic interval (segment). Note
#' that the diverse attempts to smooth or filter read-counts before
#' taking the averages did not give good results, and we remained using
#' avg="mean" and no other transformations.
#' @param rds matrix of read-count time-series
#' @param avg function to use as average statistics, \code{mean} or
#' \code{median}
#' @param mean.ratio calculate the mean ratio over time of each read-count
#' time-series before taking the total average
#' @param rm.extreme remove extreme value reads,
#' i.e., where the temporal mean is \code{<0.1} or \code{>0.9} of all means
#' @param endcut fraction of the ends the genomic interval to remove before
#' taking the total average
#' @param k parameter \code{k} for the running median
#' function \code{\link[stats]{runmed}}; the running median will be calculated
#' only if \code{k>1}
#' @param endrule parameter  \code{endrule} for the running median
#' function \code{\link[stats]{runmed}}
#' @export
segmentAverage <- function(rds, avg="mean", mean.ratio=FALSE,
                           rm.extreme=FALSE,
                           endcut=0, k=1, endrule="median") {

    ## not used - doesnt work well
    if ( mean.ratio ) {
        ## store mean
        mn <- mean(rds)
        ## mean-ratio: scale each read time-series by its mean
        rds <- rds/apply(rds,1,mean)
    }

    ## not used: cut ends
    if ( endcut>0 ) {
        N <- nrow(rds)
        end <- round(N*endcut)
        idx <- (1+end):(N-end)
        rds <- rds[idx,]
    }
    ## not used: smooth each position
    if ( k>1 )
      rds <- t(apply(rds,1,runmed,k=k,endrule=endrule))

    ## not used: remove extrema
    if ( rm.extreme ) {
        rmn <- apply(rds,1,mean)
        qu <- quantile(rmn,c(.1,.9))
        rm <- rmn<qu[1] | rmn>qu[2]
        rds[rm,] <- NA
    }

    ## not used: loess - not really working!
    if ( avg=="loess" ) { 
      y <- c(rds)
      X <- rep(1:ncol(rds),each=nrow(rds))
      data <- data.frame(x=X,y=y)
      fit <- loess(y ~ X, data=data, degree = 2) # slow!
      #plot(y ~ x, data=data,pch=19,cex=0.1)
      #j <- order(data$x)
      #lines(data$x[j],fit$fitted[j],col="red",lwd=3)
      return(predict(fit,1:ncol(rds)))
    }
    ##matplot(t(rds),type="l",lty=1,lwd=.5,col="black",log="")
    ##lines(apply(rds,2,mean,na.rm=T),col=2,lwd=3)

    ## not used (see above): ... and scale by real mean
    if ( mean.ratio ) avg <- avg*mn
    
    ## return average (mean or median)
    ## median not used
    apply(rds,2,get(avg,mode="function"), na.rm=TRUE)
}

## CLUSTER ANALYSIS

## copied from gplots to avoid dependency just for this
my.colorpanel <- function (n, low, mid, high) 
{
    odd <- function (x) x%%2 == 1 # from gtools
    if (missing(mid) || missing(high)) {
        low <- col2rgb(low)
        if (missing(high)) 
            high <- col2rgb(mid)
        else high <- col2rgb(high)
        red <- seq(low[1, 1], high[1, 1], length = n)/255
        green <- seq(low[3, 1], high[3, 1], length = n)/255
        blue <- seq(low[2, 1], high[2, 1], length = n)/255
    }
    else {
        isodd <- odd(n)
        if (isodd) {
            n <- n + 1
        }
        low <- col2rgb(low)
        mid <- col2rgb(mid)
        high <- col2rgb(high)
        lower <- floor(n/2)
        upper <- n - lower
        red <- c(seq(low[1, 1], mid[1, 1], length = lower), seq(mid[1, 
            1], high[1, 1], length = upper))/255
        green <- c(seq(low[3, 1], mid[3, 1], length = lower), 
            seq(mid[3, 1], high[3, 1], length = upper))/255
        blue <- c(seq(low[2, 1], mid[2, 1], length = lower), 
            seq(mid[2, 1], high[2, 1], length = upper))/255
        if (isodd) {
            red <- red[-(lower + 1)]
            green <- green[-(lower + 1)]
            blue <- blue[-(lower + 1)]
        }
    }
    rgb(red, blue, green)
}



### TODO - NOT WORKING CODE

## more general approach (faster/better?) of segmentOverlap
## loops over segments instead!
## Here, we join query and target segments and get the order
## of all start and end sides. Then, we go through queries. If a query
## has the same order of start and end sites, then it does NOT overlap
## with any target. [problem: order of equal positions?]
## TODO - NOT WORKING: atm this does NOT find all overlaps!
## requires additional sort of query/target start vs. target/query ends
## to find overlapping - not sure whether this would be any time improvement
## see e.g BIT algo http://europepmc.org/articles/PMC3530906
segmentOverlap.v2 <- function(query, target, details=FALSE, add.na=FALSE) {
    cols <- c("start","end")
    all<-rbind(query[,cols],target[cols]);
    n <- nrow(query)
    m <- nrow(all)
    stord <- order(all[,"start"])
    enord <- order(all[,"end"])
    ## search and classify all insert between queries 1:n
    overlaps <- data.frame()
    ## go through queries and get overlapping targets
    all <- NULL
    for ( j in 1:n ) {
        i <- which(stord == j)
        k <- which(enord == j)
        if ( i<k ) {
            ovl <- stord[(i+1):k] - n
            #ovl2 <- enord[(i-1):k] - n
            ## TODO: instead search end of query in target starts
            ## and start of query in target ends
        
            ## details
            ## scan each target
            qu <- query[j,"start"]:query[j,"end"]
            res <- NULL
            for ( t in ovl ) {
                tg <- target[t,"start"]:target[t,"end"]
                is <- intersect(tg,qu) # intersection
                len <- length(is)   
                totlen <- diff(range(c(qu,tg))) +1 # union
                ## store overlaps of target k with query
                result<-c(query=j, target=t, intersect=len,
                          ratio=round(length(qu)/length(tg),3),#query/target
                          jaccard=round(len/totlen,3), #intersect/union (Jaccard)
                          tgcov=round(len/length(tg),3), #intersect/target
                          qucov=round(len/length(qu),3)) #intersect/query
                ## is isect left/right/inside/over target?
                if ( details ) {}
                res <- rbind(res,result)
            }
            if ( !is.null(res) ) 
                all <- rbind(all,res)
            ## TODO: this doesn't cover all!! diff. stord and enord
            ovl <- cbind(query=rep(j,length(ovl)),target=ovl)
            overlaps <- rbind(overlaps,ovl)
        }        
    }
    #overlaps[order(overlaps[,"target"]),]
    rownames(all) <- NULL
    all[order(all[,"target"]),]
}

### PRE-SEGMENTATION

#' Pre-segmentation of the time-series
#' 
#' pre-segmentation of the time-series into chunks that can
#' be handled by segmenTier.
#' @param ts the time-series of readcounts for the complete chromosome,
#' rows are chromosomal positions and columns are time-points; reverse
#' strand rows at the bottom of the matrix. Option \code{chrS} can be
#' used to handle chromosome ends and to optionally (\code{map2chrom})
#' map the resulting primary segment coordinates to chromosome coordinates.
#' @param avg the broad moving average of read-count presence
#' (number of time-points with >0 reads) for a first broad segmentation
#' @param minrd the minimal number of time-points with reads in the broad
#' moving average used as cutoff between segments
#' @param favg as \code{avg}, but a narrower moving average used in
#' end scanning that can result in fusing back segments w/o good separation
#' @param border string indicating whether to "expand" or "trim" borders,
#' using the finer moving average in \code{favg}
#' @param minds minimum distance between two segments (will be fused
#' if distance is smaller)
#' @param rmlen minimum segment length for removal (shorter segments
#' will be dropped)
#' @param minsg minimum segment length for fusion with (shorter) neighbor
#' @param chrS a chromosome index, indicating at wich positions
#' chromosomes start; this is required for handling chromosome ends
#' and forward and reverse strand values, but can be omitted.
#' @param map2chrom if true, argument \code{chrS} is required to map
#' the segment coordinates to chromosomal coordinates
#' @param seg.path a directory path where individual segments' data will
#' be written to as tab-delimited .csv files; no files will be written if
#' \code{seg.path} is not provided.
#' @param plot.borders logical indicating whether plots of the scanned
#' borders should be generate (in directory at \code{fig.path})
#' @param fig.path a directory path for plots of the segment end scanning;
#' no figures will be plotted if \code{fig.path} is not provided.
#' @param fig.type image type, "png" or "pdf"
#' @param verb integer level of verbosity, 0: no messages, 1: show messages
#' @export
presegment <- function(ts, avg=1000, minrd=8,
                       favg=100, border=c("expand","trim"),
                       minds=250, rmlen=250, minsg=5e3,
                       chrS, map2chrom=FALSE,
                       seg.path, plot.borders=FALSE, fig.path, fig.type="png",
                       verb=1) {

    if ( verb> 0 )
        cat(paste("Calculating total read-counts and moving averages...\n"))

    ## plot borders fig path
    if ( missing(fig.path) )
        fig.path <- getwd() # plot to working directory

    ## cut at chromosome ends depends on presence of chrS index
    cutChromosomes <- !missing(chrS)
    
    ## total time series
    numts <- rowSums(ts > 0) ## number timepoints with reads
    
    ## moving averages of read-count presence 
    avgts <- ma(numts,n=avg,circular=TRUE) # long mov.avg, initial split
    if ( favg>1 )
        avgfn <- ma(numts,n=favg,circular=TRUE) # short mov.avg, end extension
    else avgfn <- numts
    
    ## main primary segment definition!
    ## TODO: should be >= for "minimal"
    segs <- avgts > minrd # smoothed total read number is larger then threshold
    ## set chromosome ends to FALSE as well
    if ( cutChromosomes ) {
        chrends <- sort(c(chrS[2:length(chrS)],(chrS+1)[2:length(chrS)-1]))
        if ( idx2str(length(segs),chrS)==-1 ) # add reverse strand splits
            chrends <- c(chrends, chrends + max(chrS))
        segs[chrends] <- FALSE
    }
    
    ## areas below expression threshold
    empty <- which(!segs) 

    ## distance between empty areas
    emptycoor <- empty
    if ( cutChromosomes ) # accounts for chromosome ends - TODO: could be faster
        emptycoor <- idx2coor(empty, chrS)[,"coor"]
    distn <- diff(emptycoor) 

    ## PRE-SEGMENTS: ALL WHERE DISTANCE BETWEEN EMPTY AREAS IS >1
    start <- empty[which(distn>1)]+1
    end <- start + distn[which(distn>1)]-2

    ## fuse close segments, < minds
    ## NOTE: this will fuse chromosome ends again
    ## so perhaps the effort above is not required
    close <- start[2:length(end)] - end[2:length(end)-1] < minds

    if ( verb>0 )
      cat(paste("Fusing close segments, distance <",minds,
                "bp\t", sum(close), "\n",sep=""))

    start <- start[c(TRUE,!close)]
    end <- end[c(!close,TRUE)]
    ## remove too small segments
    small <- end-start < rmlen
    if ( verb>0 )
      cat(paste("Removing small segments, <",rmlen,
                "bp\t",sum(small),"\n",sep=""))
    start <- start[!small]
    end <- end[!small]
    primseg <- cbind(start,end)


    ## (2) expand ends in both directions until mov.avg. (n=10) of signal is 0
    ## TODO: analyze gradients and minima, and add to appropriate segments
    ## TODO: border expand | trim
    if ( border[1]=="expand" ) {
        if ( verb>0 )
            cat(paste("Scanning borders ...\n"))
        fused <- 0
        for ( sg in 2:nrow(primseg) ) {
            rng <- primseg[sg-1,2]:primseg[sg,1]
            ## scan from both sides, and if overlapping
            ## take minimum between
            k <- j <- min(rng)
            i <- max(rng)
            while ( k<(i-1) ) { # expand left segment to right
                if ( avgfn[k]==0 ) break
                k <- k+1
            }
            while ( i>(j+1) ) { # expand right segment to left
                if ( avgfn[i]==0 ) break
                i <- i-1
            }
            primseg[sg-1,2] <- k
            primseg[sg,1] <- i
            
            if ( i <= k )  {
                if (  verb > 1 )
                    cat(paste("segment #",sg,i-k,
                              "overlap, will be fused with",sg-1,"\n"))
                fused <- fused +1
            }
            
            if ( !plot.borders ) next
    
            ## plot borders
            bord <- range(rng)
            rng<- max(rng[1]-1000,1):(rng[length(rng)]+1000)
            file.name <-file.path(fig.path,
                                  ifelse(i<=k,
                                         paste("fused_",fused,sep=""),
                                         paste("border_",sg-1-fused,sep="")))
            plotdev(file.name,width=4,height=4,type=fig.type)
            plot(rng,numts[rng],type="l",ylim=c(-2,24),
                 main=ifelse(i<=k,"fuse",""))
            lines(rng,avgts[rng],col=3)
            lines(rng,avgfn[rng],col=2);
            graphics::abline(h=8,col=3)
            if ( bord[1]==bord[2] ) {
                #cat(paste(sg, "borders equal\n"))
                graphics::abline(v=bord[1],col=3)
            }else
                graphics::arrows(x0=bord[1],x1=bord[2],y0=-2,y1=-2,col=3)
            if ( k==i ) {
                #cat(paste(sg, "k==i equal\n"))
                graphics::abline(v=k,col=3)
            } else
                graphics::arrows(x0=k,x1=i,y0=-1,y1=-1,col=2)
            dev.off()
        }
    
        ## (3) fuse primary segments with distance <=1, incl.
        ## those where ends where swapped in end extension
        start <- primseg[,1]
        end <- primseg[,2]
        close <- start[2:length(end)] - end[2:length(end)-1] < 2
        
        if ( verb>0 )
            cat(paste("Fusing close segments\t", sum(close), "\n",sep=""))
    
        start <- start[c(TRUE,!close)]
        end <- end[c(!close,TRUE)]

        ## NOTE: primseg v4 - primseg v3 can be reproduced with minsg==1
        ## TODO: account for chromosome ends!
        ## (4) recursively fuse short (<minsg) to shorter neighbor
        if ( verb>0 )
            cat(paste("Fusing small segments, <",minsg,
                      "bp\t", sum(end - start +1 < minsg), "\n",sep=""))
    
        while( sum(end - start +1 < minsg) ) {
            sglen <- end - start +1
            idx <- which.min(sglen)
            if ( idx==1 ) # first: fise with next
                prev <- FALSE 
            else if ( idx==length(start) ) # last: only previous possible
                prev <- TRUE
            else # take shorter neighbor
                prev <- sglen[idx-1]<sglen[idx+1]
            ## fuse
            if ( prev ) {
                start <- start[-idx]
                end <- end[-(idx-1)]
            } else {
                start <- start[-(idx+1)]
                end <- end[-idx]
            }
        }
    } else if ( border=="trim" ) {
        ## use favg to trim ends
        if ( verb>0 )
            cat(paste("Trimming segment ends",nrow(primseg),"\n"))
        for ( sg in 1:nrow(primseg) ) {
            rng <- primseg[sg,1]:primseg[sg,2]
            idx <- which(avgfn[rng]>0)
            if ( length(idx) ) {
                primseg[sg,1] <- rng[idx[1]]
                primseg[sg,2] <- rng[idx[length(idx)]]
            } else primseg[sg,] <- NA
        }
        ## rm segments below threshold
        if ( verb>0 )
            cat(paste("\tremoving",sum(is.na(primseg[,1])),"\n"))
        primseg <- primseg[!is.na(primseg[,1]),]
        start <- primseg[,1]
        end <- primseg[,2]
    }


    ## TODO: split too long segments
    ## while ( sglen>maxsg ) {}

    ## TODO: assign small intersegment to smaller adjacent segments

    ## (5) split chromosome ends!
    ## TODO: why is multiple chromosome end handling required?
    ## TODO: attach small telomeric segments to next?
    ## TODO: doesnt work at chrXVI ??
    ## get chromosomes of starts and ends via chrS
    if ( cutChromosomes ) {

        seg <- splitsegs(cbind(start,end),chrS,verb=verb)
        start <- seg[,"start"]
        end <- seg[,"end"]

        ## remove too small segments again
        small <- end-start < rmlen
        if ( verb>0 )
          cat(paste("Removing small segments, <",rmlen,
                    "bp\t", sum(small), "\n",sep=""))
        start <- start[!small]
        end <- end[!small]
    }
   
    ## TODO (6) split too long segments!

    
    primseg <- cbind(start,end)  ## DONE - PRIMARY SEGMENT v4

    ## map back to original chromosome coordinates
    if ( cutChromosomes & map2chrom ) {
        #primseg <-cbind(start=primseg[,1], end=primseg[,2])
        primseg <- index2coor(primseg,chrS)
    }
    ## write out data for each segment, if requested
    if ( !missing(seg.path) ) {
        if ( verb>0 )
            cat(paste("Writing segment data to single files.\n"))
        writeSegments(data=ts, segments=primseg, name="primseg", path=seg.path)
    }
    primseg
}

#' splits segments that span chromosome ends
#'
#' Finds segments that span chromosomes ends and splits those
#' in two segments on each covered chromosome. The input must
#' contain columns "start" and "end"; these will be modified for
#' chromosome-spanning segments. All other entries in the matrix
#' will be copied, unless an "idcol" is specified, which will receive the
#' suffix "_2" for one of two copies.
#' @param segs a matrix of segment start and end coordinates given
#' in columns named "start" and "end"
#' @param chrS a chromosome index, indicating at wich positions
#' chromosomes start; this is required for handling chromosome ends
#' and forward and reverse strand values
#' @param idcol column holding segment IDs; IDs of split segments will
#' receive the suffix "_2" for once copy#' 
#' @param verb integer level of verbosity, 0: no messages, 1: show messages
#' @export
splitsegs <- function(segs, chrS, idcol, verb=0) {

    start <- segs[,"start"]
    end <- segs[,"end"]

    if ( any(end<start) )
        stop("splitsegs requires ordered start<end coordinates")

    ## remember all columns
    col.srt <- colnames(segs)

    ## copy IDs, split IDs will receive a suffix
    if ( missing(idcol) ) idcol <- NULL
    if ( !is.null(idcol) ) ids <- as.character(segs[,idcol])
    else ids <- as.character(1:nrow(segs))

    ## copy all other columns
    othercols <- col.srt[!col.srt%in%c("start","end",idcol)]
    if ( length(othercols)==0 ) othercols <- NULL
        
    ## get strand
    schr <- idx2chr(start,chrS) # forward strand
    echr <- idx2chr(end,chrS) # reverse strand
    splt <- which(echr!=schr) # which are spanning chromosome ends?
    
    if ( verb>0 )
        cat(paste("Splitting chromosome ends",length(splt),"\n"))
    
    ## split chromosome-spanning segments, and fuse with rest
    old <- cbind(start[-splt], end[-splt])
    old.ids <- ids[-splt]
    if ( !is.null(othercols) )
        old.data <- segs[-splt,othercols,drop=FALSE]
    
    ## construct split coordinates
    str <- idx2str(start[splt],chrS)
    ## (str==-1)*max(chrS) adds minus strand to end
    new<-rbind(cbind(start[splt],chrS[schr[splt]+1] + (str==-1)*max(chrS)),
               cbind(chrS[schr[splt]+1]+1 + (str==-1)*max(chrS) ,end[splt]))
    ## add suffix to ID
    new.ids <- paste0(rep(ids[splt],2),rep(c("","_2"),each=length(splt)))

    ## copy all other columns
    if ( !is.null(othercols) ) { 
        new.data <- rbind(segs[splt,othercols,drop=FALSE],
                          segs[splt,othercols,drop=FALSE])
    }
        
    ## bind
    seg <- rbind(old,new)
    ids <- c(old.ids,new.ids)
    if ( !is.null(othercols) ) data <- rbind(old.data, new.data)

    start <- seg[,1] # re-assign start/end of segments
    end <- seg[,2]
    
    ## re-order
    srt <- order(start)
    ids <- ids[srt]
    end <- end[srt]
    start <- start[srt]
    if ( !is.null(othercols) ) data <- data[srt,]

    ## collate results and reproduce input data structure
    res <- data.frame(start=start, end=end,ID=ids, stringsAsFactors=FALSE)
    if ( !is.null(othercols) )
        res <- cbind(res, data)
    colnames(res) <- c("start","end", idcol, othercols)
    res <- res[,col.srt] # re-sort as original
    res
}


#' write out segment data to individual files
#' 
#' writing out data for segments to files, used for writing the primary
#' segments by \code{\link{presegment}} to single files that are
#' then further processed by segmenTier. This avoids loading the
#' complete data set and allows parallel usage of \code{segmenTier}.
#' @param data a data matrix to which coordinates in \code{segments}
#' refer to
#' @param segments a matrix that must contain "start" and "end" (columns)
#' of segments; these coordinates will be extracted from \code{data}
#' and written to individual files, numbered by the row number in segments.
#' @param path optional output path where files will be written, if not supplied
#' files will end up in the current working directory (`getwd`)
#' @param name name used as prefix in file names 
#' @export
writeSegments <- function(data, segments, path, name="segment") {

    for ( i in 1:nrow(segments) ) {
        rng <- segments[i,"start"]:segments[i,"end"]
        ##if ( length(rng) < 100 )
        ##    cat(paste("short segment",i, ":", length(rng),"\n"))
        dat <- data[rng,]
        id <- ifelse("ID"%in%colnames(segments),
                     as.character(segments[i,"ID"]), i)
        file.name <- paste(name, "_",id,".csv",sep="")
        if ( !missing(path) )
            file.name <- file.path(path, file.name)
        write.table(dat,file.name,row.names=FALSE,sep="\t",quote=FALSE)
    }
 }
