### PLOT UTILS

#' calculate plot device dimensions for any matrix or an overlap object
#' 
#' @export
sizeOverlaps <- function(x, mai=par('mai'), w=.25, h=.25) {

    ## catch our genera
    if ( inherits(x, "clusterOverlaps") )
        x <- x$p.value
    
    dww <- sum(mai[c(2,4)]) + w*ncol(x)
    dhh <- sum(mai[c(1,3)]) + h*nrow(x)
    list(W=dww, H=dhh)
}

#' Sanitize file names.
#' @param string to sanitize
#'@export
fsanitize <- function(str) 
    gsub("\\(","o",gsub("\\]", "c", fs::path_sanitize(str)))



#' Switch between plot devices
#' @param file.name file name without suffix (.png, etc)
#' @param type plot type: png, jpeg, eps, pdf, tiff or svg
#' @param width figure width in inches
#' @param height figure height in inches
#' @param res resolution in ppi (pixels per inch), only for 'png' and
#'     'tiff'
#' @param bg background color
#' @param mai optional par("mai") setting, used for automatic plot
#'     dimensions
#' @param overlap an object of class clusterOverlaps from which plot
#'     dimensions are calculated; useful for plotOverlaps and dotplot
#' @export
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100,
                    bg="white", mai, overlap, w=.25, h=.2, verb=0) {


    ## calculate dimensions from overlap object (clusterOverlaps)
    if ( !missing(overlap) ) {
        if ( missing(mai) ) mai <- par("mai")
        dims <- sizeOverlaps(overlap, mai = mai, w=w, h=h)
        width <- dims$W
        height <- dims$H
        if ( verb >0 )
            cat(paste("width and height from clusterOverlaps object\n"))
    }
    
    file.name <- paste(file.name, type, sep=".")
    ##if ( type=="cairopdf" )
    ##    file.name <- sub("\\.cairopdf$", ".pdf", file.name)
    if ( type == "png" )
        grDevices::png(file.name, width=width, height=height, units="in",
                       res=res, bg=bg)
    if ( type == "eps" )
        grDevices::postscript(file.name, width=height, height=width,paper="special",
                              horizontal = FALSE, onefile = FALSE)
    ## NOTE 20240425: added cairopdf to allow fonts
    ## this may change pdf output of older scripts!!
    if ( type == "pdf" ) 
        grDevices::cairo_pdf(file.name, width=width, height=height, bg=bg)
    if ( type == "tiff" )
        grDevices::tiff(file.name, width=width, height=height, units="in",
                        res=res, bg=bg)
  if ( type == "svg" )
      grDevices::svg(file.name, width=width, height=height)
    if ( type == "jpeg" )
        grDevices::jpeg(file.name, width=width, height=height, units="in",
                        res=res, bg=bg)
    
    ## optionally setup mai
    if ( !missing(mai) )
        par(mai=mai)

    invisible(list(W=width, H=height))
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
#' @param shadow use \code{\link{shadowtext}} to get background colors
#' @param ... additional parameters for call to \code{text}
#' @export
figlabel <- function(text, region="figure", pos="topleft", cex=NULL,
                     shadow=FALSE, ...) {

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

    if ( shadow ) shadowtext(x1, y1, text, cex=cex, ...)
    else text(x1, y1, text, cex=cex, ...)
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

#' Scatter plot with correlation statistics.
#'
#' @param x the x coordinates of the points in the plot.
#' @param y the y coordinates of the points in the plot.
#' @param outliers vector of indices or logical vector of x,y values
#'     to exclude from correlation analysis.
#' @param classes classification of all values (length(x)), for which
#'     correlations are calculated separateley.
#' @param log a character string which contains ‘"x"’ if the x axis is
#'     to be logarithmic, ‘"y"’ if the y axis is to be logarithmic
#'     (base 10) and ‘"xy"’ or ‘"yx"’ if both axes are to be
#'     logarithmic.
#' @param na.rm remove NA values from x and y.
#' @param verb verbosity level: report removed NAs?
#' @param cor.method method to calculate correlation and p-value via
#'     \code{\link[stats:cor.test]{cor.test}}.
#' @param line.methods regression line methods, where \code{"ols"} is
#'     a linear regression (ordinary least squares) with R's
#'     \code{\link[stats:lm]{lm}} and \code{"tls"} (total least
#'     squares) is calculated manually via
#'     \code{\link[base:eigen]{eigen}} and
#'     \code{\link[stats:cov]{cov}}.
#' @param line.col colors of the regression lines drawn with
#'     \code{line.methods}.  ## @param line.type line types of the
#'     regression lines drawn with ## \code{line.methods}.
#' @param circular treat data as circular, NOT fully
#'     implemented/tested.
#' @param circular.jitter hack for equispaced phases (for
#'     circular=TRUE), adds jitter to x or y (arguments
#'     circular.jitter=c('x','y')) to avoid NA in cor.circular.
#' @param title plot correlation parameters (as in legend) on the top
#'     of the plot.
#' @param cor.legend plot correlation, p-value and slope (TLS) as a
#'     legend.
#' @param line.legend plot line fit method as a legend.
#' @param signif number of digits to be shown for p-values in the plot
#'     legend (argument \code{digits} to \code{\link{signif}}.
#' @param round number of decimal places to be shown for correlation
#'     values in the plot legend (argument \code{digits} to
#'     \code{\link{round}}.
#' @param density indicate local point densities by using
#'     \link{dense2d} instead of the base \code{\link{plot}} function.
#' @param col point color if \code{density==FALSE}.
#' @param pch point symbols.
#' @param cex point size.
#' @param legpos position of the legend.
#' @param legcex font size (cex) of the legend.
#' @param legbg background color of the legend, default: transparent
#'     via alpha=0.  ## @param col color(s) of plotted points, if
#'     \code{density=FALSE}, or ## color palette function if
#'     \code{density=TRUE} (argument ## \code{colf} to
#'     \code{\link{dense2d}}).
#' @param ... further arguments to the plotting function;
#'     \code{\link{plot}} or, if \code{density=TRUE},
#'     \code{\link{dense2d}}.
#' @export
plotCor <- function(x, y, outliers, classes,
                    cor.method=c("pearson", "kendall", "spearman"),
                    line.methods=c("ols","tls"),
                    log="",
                    na.rm=TRUE, verb=0,
                    circular=FALSE, circular.jitter='',
                    cor.legend=TRUE, line.legend=FALSE,
                    title=FALSE, xlab, ylab,
                    line.col=c(1,2), pch=20, cex=1, axes=TRUE,
                    legpos, legbg="#FFFFFF99", legcex=1,
                    signif=1, round=2, density=TRUE, col, ...) {

    
    xy <- data.frame(x=x, y=y)

    ## default xlabel: variable name
    if ( missing(xlab) ) xlab <- deparse(substitute(x))
    if ( missing(ylab) ) ylab <- deparse(substitute(y))

    ## store eliminated data
    kept <- 1:nrow(xy)

    logaxis <- FALSE
    if ( length(grep("x", log))>0 ) {
        xy$x <- log10(xy$x)
        logaxis <- axes
        axes <- FALSE
    }
    if ( length(grep("y", log))>0 ) {
        xy$y <- log10(xy$y)
        logaxis <- axes | logaxis
        axes <- FALSE
    }

    xyo <- NULL
    if ( !missing(outliers) ) {
        if ( is.logical(outliers) )
            outliers <- which(outliers)
        if ( length(outliers)>0 ) {
            xyo <- xy[outliers,]
            xy <- xy[-outliers,]
            kept <- kept[-outliers]
            if ( !missing(col) ) col <- col[-outliers]
        }
    }
    ## clean data from NA
    if ( na.rm ) {
        rmna <- apply(xy,1, function(x) any(is.na(x)))
        if ( sum(rmna)>0 & verb>0 )
            cat(paste("removing ", sum(rmna), " of ",nrow(xy),
                      " rows with NA values\n"))
        xy <- xy[!rmna,]
        kept <- kept[!rmna]
        if ( !missing(col) )
           col <- col[!rmna]

        rminf <- apply(xy, 1, function(x) any(is.infinite(x)))
        if ( sum(rminf)>0 )
            cat(paste("removing ", sum(rminf), " of ",nrow(xy),
                      " rows with INF values\n"))
        xy <- xy[!rminf,]
        kept <- kept[!rminf]
        if ( !missing(col) )
            col <- col[!rminf]
    }

    
    ## line fit and r-squared
    lfit <- tls <- beta <- crt <- NULL

    pv <- NA
    if ( circular ) {

        ## circular correlation and line fit!

        xyc <- xy
        for ( var in circular.jitter )
          xyc[[var]] <- jitter(xyc[[var]])
  
        xyc <- as.data.frame(circular::circular(xyc,
                                                units="radians",
                                                type="angles"))
        
 
        ## TODO: use bpnreg ?
        ## TODO: allow only 1D circular+ type='c-l': maximum
        ## likelihood regression model proposed by Fisher and Lee (1992),
        ## with linear variable x:
        ## lfit <- lm.circular(x=y, y=x, type='c-l', init=0)

        ## lm.circular.cc: 
        lfit <- circular::lm.circular(y=xyc$y, x=xyc$x, type="c-c", order=2)
        ## add polynomial fit
        fline <- lfit$fitted
        fline[fline>pi] <- fline[fline>pi] - 2*pi
        ##points(as.numeric(df$x), fline, col=1, pch=19,cex=.3)

        crt <- circular::cor.circular(x=xyc$x, y=xyc$y, test=TRUE)
        cr <- round(crt$cor, round)
        pv <- signif(crt$p.value, signif)
    } else {

        if ( !missing(classes) ) {}
        
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

            ## TODO: use .Machine$double.xmin (2.2e-308) or
            ## .Machine$double.eps (2.2e-16) as cutoff and plot < in legend 
            
        }
    }

    ## LEGEND EXPRESSIONS
    ## p-value legend expression
    nexpr <- cexpr <- pexpr <- texpr <- "NA"
    nexpr <- bquote(n == .(nrow(xy)))
    if ( !is.na(cr) )
        cexpr <- bquote(r == .(cr))
    if ( !is.na(pv) )
        if ( pv < .Machine$double.xmin ) {
            pexpr <- bquote(p < .(signif(.Machine$double.xmin,2)))
            texpr <- bquote(n == .(nrow(xy))*","~
                                r == .(cr)*","~
                                    p < .(signif(.Machine$double.xmin,2)))
        } else {
            pexpr <- bquote(p == .(pv))
            texpr <- bquote(n == .(nrow(xy))*","~
                                r == .(cr)*","~
                                    p == .(pv))
        }

    ## DATA
    if ( density )
        dense2d(xy$x, xy$y, circular=circular, cex=cex, pch=pch, axes=axes, ...)
    else {
        if ( missing(col) ) col <- "#000000AA"
        plot(xy$x, xy$y, pch=pch, cex=cex, col=col, axes=axes,
             xlab = xlab, ylab = ylab, ...)
    }
    ## FIT LINES
    if ( !circular ) {
        if ( "ols"%in%line.methods )
            abline(lfit, col=line.col[1])
        if ( "tls"%in%line.methods )
            abline(a=tls$alpha, b=tls$beta, col=line.col[2])
    } else {
        points(as.numeric(xy$x), fline, col=1, pch=19, cex=.3)
    }

    ## indicated outliers
    if ( !missing(outliers) )
        points(xyo$x, xyo$y, col=2, pch=4, cex=cex, ...)

    ## axes
    if ( logaxis ) {
        if ( length(grep("x", log))>0 ) {
            logaxis(1)
        } else axis(1)
        if ( length(grep("y", log))>0 ) {
            logaxis(2)
        } else axis(2)
    }

    ## TODO: replace line legend by TLS and OLS
    if ( cor.legend ) {
        ## construct legend from available data
        leg <- lty <- lcol <- lpch <- NULL
        if ( !circular ) {
            if ( cor.method[1]!="" ) {
                leg <- c(leg,
                         as.expression(nexpr),
                         as.expression(cexpr),
                         as.expression(pexpr))
                lty <- c(lty, NA, NA, NA)
                lcol <- c(lcol, NA, NA, NA)
                lpch <- c(lpch, NA, NA, NA)
                
            }
            if ( line.legend ) {
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
            }
        } else { # circular
            leg <- c(as.expression(cexpr),
                     as.expression(pexpr))
        }
        if ( !missing(outliers) ) {
            if ( length(outliers)>0 ) {
                leg <- c(leg, expression(excluded))
                lty <- c(lty, NA)
                lcol <- c(lcol, 2)
                lpch <- c(lpch, 4) 
            }
        }

        if ( all(is.na(lty)) ) lty <- NULL
        if ( length(leg)>0 )
            if ( missing(legpos) )
                legpos <- ifelse(cr<0,"topright","topleft") 
        legend(legpos,
               legend=leg,
               col=lcol, lty=lty,
               pch=lpch, seg.len=.5,
               box.col=NA, bg=legbg, cex=legcex,
               y.intersp=.75, x.intersp=0.75)
    }
    if ( title ) {
        titl <- as.expression(texpr)
        mtext(titl,3,0, cex=.9)
    }
        
    invisible(list(xy=xy, n=nrow(xy), kept=kept, r=cr, p=pv, 
                   cor=crt, fit=lfit, tls=tls))
}

#' Plot linear regressions for subsets (clusters) of data
#'
#' @param x the x coordinates of the points in the plot.
#' @param y the y coordinates of the points in the plot.
#' @param col point colors; if not provided cls.col will be used.
#' @param cls a vector of cluster labels for the xy data.
#' @param cls.col a named vector of cluster label colors.
#' @param cls.srt a sorting and sub-selection of cluster labels.
#' @param alpha color alpha level applied to cluster colors of the
#'     plotted points.
#' @param method method used to calculate the correlations.
#' @param legend add a legend showing correlations for each cluster.
#' @param leg.pos position of the legend.
#' @param verb verbosity level.
#' @param ... arguments to \link{plotCor}.
#' @export
plotCorMulti <- function(x, y, col, cls, cls.col, cls.srt,
                         alpha = 0.5, method = 'pearson', line='ols',
                         legend=TRUE, leg.pos='bottomright', verb = 1, ...) {


     if ( missing(cls.srt) ) 
         cls.srt <- unique(cls)

    if ( missing(col) ) {
        col <- add_alphas(cls.col[cls], alpha=alpha)
        col[!cls%in%cls.srt] <- NA
    } 

    ## filter all non-plotted cohorts
    keep <- cls%in%cls.srt
    if ( any(!keep) ) {
        x <- x[keep]
        y <- y[keep]
        col <- col[keep]        
        cls <- cls[keep]
    }

    
    fit <- plotCor(x, y, 
                   cor.legend=FALSE, title=FALSE,
                   density=FALSE, line.methods='',
                   col=col, ...)
    
    ## filter clustering my eliminated data
    cls <- cls[fit$kept]
    
    ## correlations for each cluster
    xy <- fit$xy
    isn <- !is.na(xy$x) & !is.na(xy$y) ##& is.finite(xy$x) & is.finite(xy$y)


    fits <- list()
    for ( cl in cls.srt ) {

        set <- cls==cl & isn
        
        if ( verb>0 )
            cat(paste('fitting', cl))
        if ( sum(set)<2 ) {
            if ( verb>0 )
                cat(paste(':\tnot enough data\n'))
            next
        } else if ( verb>0 )
            cat(paste('\n'))

        ## correlation
        cr <- cor.test(xy$y[set], xy$x[set], method = method)
        ## linear regression (ordinary least squares)
        ft <- lm(xy$y[set] ~ xy$x[set])
        ## total least squares
        v <- eigen(cov(xy[set,]))$vectors
        beta <- v[2,1]/v[1,1] # slope
        alpha <- mean(xy$y[set]) - beta*mean(xy$x[set]) # intercept
        fits[[cl]] <- list(fit=ft, cor=cr)
        
        if ( 'ols'%in%line ) {
            abline(ft, col=cls.col[cl], lwd=2)
        } 
        if ( 'tls'%in%line ) {
            abline(a=alpha, b=beta, col=cls.col[cl], lwd=2)
        }
   

    }
    if ( legend ) {
        rhos <- sapply(fits, function(x) unname(x$cor$estimate))
        legend(leg.pos, legend=paste0(names(rhos),": ", round(rhos,2)),
               seg.len=.75, 
               col=cls.col[names(rhos)], lty=1, bg='#ffffff77', box.col=NA)
    }
    invisible(fits)
}

#' Nicely formated date axis over three date ranges.
#' @param side axis side.
#' @param xlim x-axis limits (as.POSIXct), used to calculate  appropriate ticks.
#' @export
dateaxis <- function(side=1, xlim) {

    require(lubridate)


    ## from/to 1st of January
    if ( missing(xlim) ) {
        rng <- c(as.Date(paste0("0000-01-01")),
                 as.Date(paste0("2100-01-01")))
        xlim <- rng
    } else {
        if ( !inherits(xlim, "POSIXct") )
            xlim <- as.POSIXct(xlim)
        rng <- c(as.Date(sprintf("%04d-01-01",year(min(xlim)))),
                 as.Date(sprintf("%04d-01-01",year(max(xlim))+1)))
        xlim <- range(xlim)
    }
    
    ranges <- c("day","week","month","year")
    formats <- c("%m/%d","%m/%d","%m/%d","%Y")

    dt <- list()
    dt[[1]] <- seq(from=rng[1], to=rng[2],by="day")
    dt[[2]] <- seq(from=rng[1], to=rng[2],by="week")
    dt[[3]] <- seq(from=rng[1], to=rng[2],by="month")
    dt[[4]] <- seq(from=rng[1], to=rng[2],by="year")

    ## select scale based on diff
    ## 1: day/week/month
    ## 2: week/month/year
    dw <- as.numeric(difftime(xlim[2],xlim[1],unit="weeks"))
    ## select label format
    ##idx <- tail((1:5)[dw>c(0,1.5*5*10^c(0:2))],3)
    idx <- 2:4
    if ( dw<53) idx <- 1:3
    if ( dw<5) idx <- 1:2
    if ( dw<1) idx <- 1
    tcl <- par("tcl")
    for ( i in idx ) {
        frmt <- ifelse(i==tail(idx,1), formats[i], NA)
        axis.POSIXct(side, at=dt[[i]], format=frmt, labels=!is.na(frmt),
                     tcl=which(idx==i)*tcl/length(idx))
    }
}

## TODO: test this chatGPT-generated wrapper for base R's plot and
## matplot that can use ashaxis and logaxis


#' A wrapper for base R's plot.default that catches the log and adds an
#' ash argument, transforms data by log10 or arsinh, but plots the
#' original scale as logarithmic axes.
#' @inheritParams graphics::plot.default
#' @param ash used as the log argument in base R's plot.default, but
#'     using arcsinh transform instead of log.
#'@export
lplot <- function(x, y = NULL, ..., log = "",  ash = "") {

    axes_log <- unique(strsplit(log, "")[[1]])
    axes_ash <- unique(strsplit(ash, "")[[1]])
    if (length(intersect(axes_log, axes_ash))) {
        stop("Axis specified in both 'log' and 'ash'", call. = FALSE)
    }
    
    dots <- list(...)
    if (!is.null(dots$log)) log <- dots$log 
    dots$log <- NULL  # remove it so plot() never sees it

    ## Extract x/y if supplied via formula or list
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y

    ## Transform data
    if (grepl("x", log)) x <- log10(x)
    if (grepl("y", log)) y <- log10(y)
    if (grepl("x", ash)) x <- ash(x)
    if (grepl("y", ash)) y <- ash(y)

    ## Draw plot without log scaling
    plot.default(
        x, y,
        log = "",
        axes = FALSE,
        ...
    )

    ## Custom axes
    if (grepl("x", log)) logaxis(1)
    else if (grepl("x", ash)) ashaxis(1) 
    else axis(1)
    if (grepl("y", log)) logaxis(2)
    else if (grepl("y", ash)) ashaxis(2) 
    else axis(2)

}

#' Draw axis for arcsinh transformed data
#' @inheritParams graphics::axis
#' @param ... parameters passed on to \link[graphics:axis]{axis}
##xs <- c(-1e4, -10:10, 1e3)
##plot(xs, ash(xs))
##ashaxis(4)
#' @export
ashaxis <- function(side, at = c(-10^(7:1), 0, 10^(1:7)), ...) {

    for ( ax in side ) {

        ## axis with tick labels
        axis(ax, at = ash(at), labels = at, ...)
        ## select small ticks w/o labels
        ## TODO: select appropriate range for small ticks,
        ## as in logaxis
        if ( FALSE ) {
            n <- 20
            rng <- ash(ceiling(par('usr')[3:4]))
            lat <- seq(rng[1], rng[2], length.out=n)
            axis(ax, at = lat, labels=FALSE, tcl=par('tcl')/2)
        }
    }
}

## TODO: fix labels for base=exp(1)!!
##
#' Logarithmic axis ticks
#' @inheritParams graphics::axis
#' @param lat as argument \code{at} of \link[graphics:axis]{axis} but
#'     as the exponents of \code{log(y, base=base)}
#' @param base base of the logarithm
#' @param verb verbosity level
#' @export
logaxis <- function(side, lat=-100:100, base=10, labels, verb=0, ...) {

    ## handle lables argument
    if ( missing(labels) ) labels <- base^lat
    ## FALSE should suppress labels, TRUE: overwrite with labels
    if ( is.logical(labels) ) if ( labels ) labels <- base^lat

    ## plot all sides; NOTE: vectorized
    for ( ax in side ) {
        
        axis(ax, at=lat, labels=labels, ...)

        ## log-distance ticks at half tick length
        ## TODO: fix minor tick marks for log bases other than 10
        if ( base==10 ) {
          ##lticks<- log10(rep(1:10,  length(lat)) *   10^rep(lat-1, each=10))
            lticks <- log(rep(1:base, length(lat)) * base^rep(lat-1, each=base),
                          base=10)
            axis(ax, at=lticks,
                 tcl=par('tcl')/2, labels=FALSE, ...)
        } else if ( verb>0) {
            cat(paste('minor ticks marks only implement for log-base 10\n'))
        }
    }
}



## chatGPT; TODO: test and fix
logaxis2 <- function(side, lat = -10:10, base = 10, labels, ...) {
  if (missing(labels)) 
    labels <- base^lat
  if (is.logical(labels)) 
    if (labels) 
      labels <- base^lat
  
  for (ax in side) {
    # Major ticks
    axis(ax, at = lat, labels = labels, ...)
    
    # Minor ticks
    if (base > 1) {
      # Create minor ticks between each pair of major ticks
      minor_per_decade <- base - 1
      minor_ticks <- c()
      for (i in lat) {
        minors <- log(base^i * (1:(base - 1)), base = base)
        minor_ticks <- c(minor_ticks, minors)
      }
      axis(ax, at = minor_ticks, tcl = par("tcl") / 2, labels = FALSE, ...)
    } else {
      warning("Base must be > 1 for log scale")
    }
  }
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
#' @param x.circular treat x as circular coordinates in radian.
#' @param y.circular treat y as circular coordinates in radian.
#' @param colf color map function used to create a color gradient,
#'     eg. \code{\link[grDevices:colorRampPalette]{colorRampPalette}}
#'     or \code{viridis}.
#' @param ... arguments to plot (hint:cex can be useful).
#' @return the plotted data.frame, including local densities.
#' @export
dense2d <- function(x, y, pch=20, nbin=c(128,128), circular=FALSE,
                    x.circular=FALSE, y.circular=FALSE,
                    xlim, ylim,
                    colf=viridis::viridis,
                    ...) {

    ## NOTE, previous default color function was
    ## grDevices::colorRampPalette(c("#000099","#00FEFF",
    ##                               "#45FE4F", "#FCFF00",
    ##                               "#FF9400", "#FF3100"))

    
    xa <- x
    ya <- y

    if ( circular ) x.circular <- y.circular <- TRUE
    if ( x.circular & !y.circular ) {
        ## add data +/- 2*pi to get circular density
        xa <- c(x-2*pi,x,x+2*pi)
        ya <- c(y,y,y)
        if ( missing(xlim) ) xlim <- c(-pi,pi)
    } else if ( y.circular & !x.circular ) {
        ya <- c(y,
                y-2*pi,
                y+2*pi)
        xa <- c(x,x,x)        
        if ( missing(ylim) ) ylim <- c(-pi,pi)
    } else if ( circular ) {
        xa <- c(x-2*pi,x,x+2*pi,
                x-2*pi,x,x+2*pi,
                x-2*pi,x,x+2*pi)
        if ( missing(xlim) ) xlim <- c(-pi,pi)
        ya <- c(y,y,y,
                y-2*pi,y-2*pi,y-2*pi,
                y+2*pi,y+2*pi,y+2*pi)
        if ( missing(ylim) ) ylim <- c(-pi,pi)
    }

    df <- data.frame(x=xa, y=ya)

    if ( missing(xlim) ) xlim <- range(x[is.finite(x)], na.rm=TRUE)
    if ( missing(ylim) ) ylim <- range(y[is.finite(y)], na.rm=TRUE)
    
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


#' Color Ramp Palette `arno`
#'
#' generates colors similar to viridis::inferno but with a more visible
#' yellow, taken from viridis::viridis. Named after Arno because it
#' was his birthday.
#' @param n number of colors to return
#' @export
arno <- function(n) {
    
    mcol <- viridis::inferno(5)
    vcol <- viridis::viridis(5)
    mcol[5] <- vcol[5]
    colorRampPalette(mcol)(n)
}


#' R colors with alpha value.
#'
#' @param col and R color specifications such as numeric string colors.
#' @param alpha the desired alpha (opaqueness) level from 0
#'     (transparent) to 1 (opaque).
#' @export
col2alpha <- function(col='red', alpha=1)  {
    
    apply(grDevices::col2rgb(col, alpha=FALSE), 2, function(x)
        grDevices::rgb(x[1]/255, x[2]/255, x[3]/255, alpha))
}
