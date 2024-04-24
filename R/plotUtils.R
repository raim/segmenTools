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


#' Add sub-figure label to corners of plots
#'
#' copied from code by January Weiner 3rd at
#' \url{https://www.r-bloggers.com/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/},
#' originally published at
#' \url{https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/}
#' @param text figure label
#' @param region add to a figure, plot or device
#' @param pos position key, as in \code{legend}
#' @param cex label size
#' @param ... additional parameters for call to \code{text}
#' @export
figlabel <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)

  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}



#' obselete, use \code{\link{num2col}} 
#'
#' Maps values to colors.
#' A simple function copied from 
## from https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r to map numeric values to colors
#' @param x numeric values
#' @param cols a vector of colors
#' @param limits min/max cut-off values of x
#' @export
val2col <-function(x, cols,limits=NULL) {
    
    if(is.null(limits)) limits=range(x)

    cols[findInterval(x,seq(limits[1],limits[2],
                            length.out=length(cols)+1), all.inside=TRUE)]
}

#' convert numeric values to color range:
#' @param x a numeric vector
#' @param limits optional data limits for min/max color, every x
#'     lower/higher will get the extreme colors
#' @param q auto-select limits by these two quantiles (argument
#'     \code{prob} in function \code{\link[stats:quantile]{quantile}}
#' @param pal color palette, alternatively \code{colf} and \code{n}
#'     can be supplied
#' @param colf color palette function
#' @param n number of different colors
#' @export
num2col <- function(x, limits, q, pal, colf=viridis::viridis, n=100){
    if ( missing(pal) ) pal <- colf(n)
    if ( missing(limits) ) limits <- range(x, na.rm=TRUE)
    if ( !missing(q) ) {
        if ( length(q)!=2 )
            stop("argument q must be of length 2")
        limits <- quantile(x, probs=q)
    }
    pal[findInterval(x,seq(limits[1],limits[2],
                           length.out=length(pal)+1), all.inside=TRUE)]
}

#' scatter plot with correlation statistics
#'
#' @param x the x coordinates of the points in the plot.
#' @param y the y coordinates of the points in the plot.
#' @param outliers vector of indices or logical vector of x,y values to
#' exclude from correlation analysis.
#' @param na.rm remove NA values from x and y.
#' @param cor.method method to calculate correlation and p-value via
#'     \code{\link[stats:cor.test]{cor.test}}.
#' @param line.methods regression line methods, where \code{"ols"} is
#'     a linear regression (ordinary least squares) with R's
#'     \code{\link[stats:lm]{lm}} and \code{"tls"} (total least squares)
#'     is calculated manually via \code{\link[base:eigen]{eigen}} and
#'     \code{\link[stats:cov]{cov}}.
#' @param line.col colors of the regression lines drawn with
#'     \code{line.methods}.
## @param line.type line types of the regression lines drawn with
##     \code{line.methods}.
#' @param circular NOT fully implement, treat data as circular.
#' @param cor.legend plot correlation, p-value and slope (TLS) as a legend.
#' @param signif number of digits to be shown for p-values in the plot
#'     legend (argument \code{digits} to \code{\link{signif}}.
#' @param round number of decimal places to be shown for correlation
#'     values in the plot legend (argument \code{digits} to
#'     \code{\link{round}}.
#' @param density indicate local point densities by using
#'     \link{dense2d} instead of the base \code{\link{plot}} function.
#' @param pch point symbols.
#' @param cex point size.
#' @param legpos position of the legend.
## @param col color(s) of plotted points, if \code{density=FALSE}, or
##     color palette function if \code{density=TRUE} (argument
##     \code{colf} to \code{\link{dense2d}}).
#' @param ... further arguments to the plotting function;
#'     \code{\link{plot}} or, if \code{density=TRUE},
#'     \code{\link{dense2d}}.
#' @export
plotCor <- function(x, y, outliers,
                    cor.method=c("pearson", "kendall", "spearman"),
                    line.methods=c("ols","tls"),
                    na.rm=TRUE, circular=FALSE,
                    cor.legend=TRUE, line.col=c(1,2), pch=20, cex=1,
                    legpos,
                    signif=1, round=1, density=TRUE, ...) {

    xy <- data.frame(x=x, y=y)

    xyo <- NULL
    if ( !missing(outliers) ) {
        if ( is.logical(outliers) )
            outliers <- which(outliers)
        xyo <- xy[outliers,]
        xy <- xy[-outliers,]
    }
    ## clean data from NA
    if ( na.rm ) {
        rmna <- apply(xy,1, function(x) any(is.na(x)))
        if ( sum(rmna)>0 )
            warning("removing ", sum(rmna), " of ",nrow(xy),
                    " rows with NA values")
        xy <- xy[!rmna,]
    }

    
    ## line fit and r-squared
    lfit <- tls <- beta <- crt <- NULL

    if ( circular ) {
        xyc <- as.data.frame(circular(xy, units="radians", type="angles"))

        ## TODO: implement circular correlation and line fit!
        ## TODO: convert to circular first
        ## TODO: bpnreg
        lfit <- circular::lm.circular(y=xyc$y, x=xyc$x, type="c-c", order=2)
        ## add polynomial fit
        fline <- lfit$fitted
        fline[fline>pi] <- fline[fline>pi] - 2*pi
        ##points(as.numeric(df$x), fline, col=1, pch=19,cex=.3)

        crt <- circular::cor.circular(x=xyc$x, y=xyc$y, test=TRUE)
        cr <- round(crt$cor, round)
        pv <- signif(crt$p.value, signif)
    } else {
        ## lin.reg (OLS)
        if ( "ols" %in% line.methods )
            lfit <- lm(y ~ x, data=xy)
        ## TLS, via eigen(cov)
        if ( "tls" %in% line.methods ) {
            v <- eigen(cov(xy))$vectors
            beta <- v[2,1]/v[1,1] # slope
            alpha <- mean(xy$y) - beta*mean(xy$x) # intercept
            tls <- list(alpha=alpha, beta=beta)
            ## for legend
            sl <- round(beta, round)
        }
        ## correlation
        if ( cor.method[1]!="" ) {
            crt <- cor.test(xy$x, xy$y,
                            method=cor.method, use="pairwise.complete")
            ## for legend
            cr <- round(crt$estimate, round)
            pv <- signif(crt$p.value, signif)
        }
    }
    
    if ( density )
        dense2d(xy$x, xy$y, circular=circular, cex=cex, pch=pch, ...)
    else
        plot(xy$x, xy$y, pch=pch, cex=cex, ...)
    if ( !circular ) {
        if ( "ols"%in%line.methods )
            abline(lfit, col=line.col[1])
        if ( "tls"%in%line.methods )
            abline(a=tls$alpha, b=tls$beta, col=line.col[2])
    } else {
        points(as.numeric(xy$x), fline, col=1, pch=19, cex=.3)
    }

    if ( !missing(outliers) )
        points(xyo$x, xyo$y, col=2, pch=4, cex=cex, ...)

    ## TODO: replace line legend by TLS and OLS
    if ( cor.legend ) {
        ## construct legend from available data
        leg <- lty <- lcol <- lpch <- NULL
        if ( !circular ) {
            if ( cor.method[1]!="" ) {
                leg <- c(leg,
                         as.expression(bquote(r == .(cr))),
                         as.expression(bquote(p == .(pv))))
                lty <- c(lty, NA, NA)
                lcol <- c(lcol, NA, NA)
                lpch <- c(lpch, NA, NA)
                
            }
            if ( "ols" %in% line.methods ) {
                leg <- c(leg, expression(OLS))
                lty <- c(lty, 1)
                lcol <- c(lcol, line.col[1])
                lpch <- c(lpch, NA)
            }
            if ( "tls" %in% line.methods ) {
                leg <- c(leg, expression(TLS))
                lty <- c(lty, 1)
                lcol <- c(lcol, line.col[2])
                lpch <- c(lpch, NA)
            }
        } else { # circular
            leg <- c(as.expression(bquote(r == .(cr))),
                     as.expression(bquote(p == .(pv))))
        }
        if ( !missing(outliers) ) {
            leg <- c(leg, expression(excluded))
            lty <- c(lty, NA)
            lcol <- c(lcol, 2)
            lpch <- c(lpch, 4) 
        }
        if ( length(leg)>0 )
            if ( missing(legpos) )
                legpos <- ifelse(cr<0,"topright","topleft") 
            legend(legpos,
                   legend=leg,
                   col=lcol, lty=lty, pch=lpch, seg.len=.5, bty="n",
                   y.intersp=.75, x.intersp=0.75)
    }
    invisible(list(xy=xy, cor=crt, fit=lfit, tls=tls))
}

#' 2D density heatmap plot
#'
#' Uses base R's \code{\link[grDevices:densCols]{densCols}} to
#' calculate local densities at each point in a scatterplot, and then
#' replaces them by a colored scheme, copied from Josh O'Brien posted
#' at
#' https://stackoverflow.com/questions/17093935/r-scatter-plot-symbol-color-represents-number-of-overlapping-points/17096661#17096661
#' @param x x-coordinates
#' @param y y-coordinates
#' @param xlim, numeric vectors of length 2, giving the x 
#'     coordinate ranges.
#' @param ylim as \code{xlim} for y.
#' @param pch \code{pch} argument to plot
#' @param nbin number of bins for both dimensions, can be a single.
#'     number for both dimensions, or separate numbers, see argument
#'     \code{nbins} to function
#'     \code{\link[grDevices:densCols]{densCols}} and \code{gridsize}
#'     to \code{\link[KernSmooth:bkde2D]{bkde2D}}.
#' @param circular treat x and y as circular coordinates in radian.
#' @param colf color map function used to create a color gradient,
#'     eg. \code{\link[grDevices:colorRampPalette]{colorRampPalette}}
#'     or \code{viridis}.
#' @param ... arguments to plot (hint:cex can be useful).
#' @return the plotted data.frame, including local densities.
#' @export
dense2d <- function(x, y, pch=20, nbin=c(128,128), circular=FALSE,
                    xlim, ylim,
                    colf=viridis::viridis,
                    ...) {

    ## NOTE, previous default color function was
    ## grDevices::colorRampPalette(c("#000099","#00FEFF",
    ##                               "#45FE4F", "#FCFF00",
    ##                               "#FF9400", "#FF3100"))
    
    if ( circular ) {
        ## add data +/- 2*pi to get circular density
        xa <- c(x-2*pi,x,x+2*pi,
                x-2*pi,x,x+2*pi,
                x-2*pi,x,x+2*pi)
        ya <- c(y,y,y,
                y-2*pi,y-2*pi,y-2*pi,
                y+2*pi,y+2*pi,y+2*pi)
        if ( missing(xlim) ) xlim <- c(-pi,pi)
        if ( missing(ylim) ) ylim <- c(-pi,pi)
        df <- data.frame(x=xa, y=ya)
        
    } else {
        df <- data.frame(x=x, y=y)
        if ( missing(xlim) ) xlim <- range(x, na.rm=TRUE)
        if ( missing(ylim) ) ylim <- range(y, na.rm=TRUE)
    }
    
    ## TODO: nbin = 128 in densCols, vs. 256 colors in colorRampPalette
    ## KernSmooth::bkde2D(cbind(x, y), bandwidth=c(20,20), gridsize=c(128,128))
    ## Use densCols() output to get density at each point
    xcol <- grDevices::densCols(df$x, df$y, nbin=nbin, 
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
    plot(y~x, data=df[order(df$dens),], pch=pch, col=col,
         xlim=xlim, ylim=ylim, ...)
    
    ## normalize to 1, for legend?
    df$dens <- df$dens/ncol
    
    invisible(df)
}
