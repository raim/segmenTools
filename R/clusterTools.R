
### TOOLS for handling feature-based clusterings

## TODO:
## calculate average phases via phaseDist or directly
## auto-select oscillatory clusters from p-values
## clusterCluster/image_matrix : make higher level interface, see
## plot.hypergeoTables in analyzeSegments_2017.R
## integrate GO wrapper for clusterCluster

### COMPARE CLUSTERINGS 

#' Wilcoxon test wrapper that returns a normalized U-statistic
#'
#' This wrapper just calls \link[stats:wilcox.test]{wilcox.test} but
#' normalizes the test's U-statistic by the product of the lengths
#' (counts of non-NA values) of x and y such that values U>0 indicate
#' that x is generally larger than y and opposite for U<0.
#' @param x,y as x and y passed to
#'     \link[stats:wilcox.test]{wilcox.test}
#' @param ... further parameters passed on to
#'     \link[stats:wilcox.test]{wilcox.test}
#' @export
w.test <- function(x, y, ...) {
    res <- stats::wilcox.test(x, y, ...)
    ## normalized U-statistic
    tt <- res$statistic/(sum(!is.na(x))*sum(!is.na(y))) -0.5
    rt <- list()
    rt$statistic <- unlist(tt)
    rt$p.value <- unlist(res$p.value)
    rt
}

#' Calculate biases of numerical data between clusterings.
#'
#' Calculate t-tests or wilcoxon (rank sum) tests for each cluster in
#' a clustering against the total distributions. The resulting tables
#' can be plotted with \code{\link{plotOverlaps}} and sorted along one
#' axis by signficance with \code{\link{plotOverlaps}}.
#' @param x matrix of numeric values where columns are different data
#'     sets and rows must correspond to the clustering in argument
#'     \code{cls}
#' @param cls a clustering of the rows in argument \code{x} or a
#'     logical TRUE/FALSE table (matrix) with cluster labels as column
#'     names and rows corresponding to the rows in argument \code{x}
#' @param test test function to be applied, default is
#'     \code{\link[stats]{t.test}}. This can be any function that
#'     takes the cluster subset of \code{x[cls=<cl>,]} as first and
#'     the total distribution \code{x} as second argument, and returns
#'     an object with items \code{statistic} and
#'     \code{p.value}. Negative values in \code{statistic} will be
#'     differently colored in \code{\link{plotOverlaps}} and the sign
#'     copied to \code{p.value}. This package provides
#'     \code{\link{w.test}}, a wrapper for
#'     \code{\link[stats:wilcox.test]{wilcox.test}} that
#'     normalizes the tests' U-statistics such that it can be used in
#'     the same way as the t-statistic.
#' @param min.obs minimal number of non-NA observations
#' @param replace test with replacement, i.e., each cluster is tested
#'     against the whole data set, including the cluster items ##
#'     TODO: instead of passing function, pass a type and handle
#'     numbers ## betterer, eg. normalized U-statistic for wilcox - OR
#'     handle this ## in plotOverlaps by statistic type!
#' @export
clusterProfile <- function(x, cls, test=stats::t.test, min.obs=5,
                           replace=FALSE) {

    if ( !inherits(cls, "matrix") & !inherits(cls,"factor") )
        cls <- factor(cls, levels=unique(cls))
    logic <- FALSE
    if ( inherits(cls, "factor") )
        cls.srt <- levels(cls)
    else if ( inherits(cls, "matrix") ) {
        cls.srt <- colnames(cls)
        logic <- TRUE
    }


    if ( inherits(test, "character") )
        test <- get(test, mode="function")
    
    ## t-test statistic matrix
    tt <- matrix(0, ncol=ncol(x), nrow=length(cls.srt)) 
    colnames(tt) <- colnames(x)
    rownames(tt) <- cls.srt #levels(cls)
    
    tp <- sg <- tt
    tp[] <- 1 # p-value matrix
    sg[] <- 1 # sign matrix
    for ( i in 1:ncol(x) ) {
        for ( cl in cls.srt ) {

            if ( logic ) { # logic table
                y <- x[cls[,cl],i]
                if ( replace )
                    X <- x[,i]
                else X <- x[!cls[,cl],i] 
            } else { # clustering/factors
                y <- x[which(cls==cl),i]
                if ( replace )
                    X <- x[,i]
                else X <- x[which(cls!=cl),i]
            }
            if ( sum(!is.na(y))<min.obs ) next
            ttmp <- test(y, X)
            tt[cl,i] <- ttmp$statistic
            ## set p-values negative for two-sided plot!
            sgn <- sign(ttmp$statistic)
            sgn[sgn==0] <- 1
            sg[cl,i] <- sgn # store sign
            tp[cl,i] <- ttmp$p.value # store pvalue
        }
    }
    ## construct overlap object
    ova <- list()
    ova$p.value <- tp
    ova$statistic <- tt
    ova$sign <- sg

    ## counts
    ova$num.target <- t(as.matrix(apply(x,2,function(x) sum(!is.na(x)))))
    if ( logic )
        ova$num.query <- as.matrix(apply(cls,2,sum))
    else 
        ova$num.query <- as.matrix(table(cls)[levels(cls)])

    class(ova) <- "clusterOverlaps"
    ova
}


#' calculates overlaps between two clusterings
#' 
#' Calculates mutual overlaps between two clusterings of the same data set
#' using hypergeometric distribution statistics for significantly
#' enriched or deprived mutual overlaps. The resulting tables can
#' be plotted with \code{\link{plotOverlaps}} and sorted along one
#' axis by signficance with \code{\link{plotOverlaps}}.
#' TODO: specify wich cl will be rows/columns and to which
#' percent refers to
#' @param cl1 clustering 1; a vector of cluster associations or an
#' object of class "clustering" by segmenTier's
#' \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param cl2 clustering 2; see argument \code{cl1}
#' @param na.string replace NA or empty strings by `<na.string>'
#' @param cl1.srt optional cluster sorting of clustering 1
#' @param cl2.srt optional cluster sorting of clustering 2
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of `"greater"' to calculate enrichment (default),
#' `"less"' to calculate deprivation, or `"two.sided"' to report the more
#' more signficant p-value of both
#'@export
clusterCluster <- function(cl1, cl2, na.string="na", cl1.srt, cl2.srt,
                           alternative=c("greater")) {

   # read and return cluster variable names
   # can be important for plot Overlaps
   vn1=deparse(substitute(cl1))
   vn2=deparse(substitute(cl2))	 	
	 
    if ( inherits(cl1, "clustering") ) {
        K <- selected(cl1)
        if ( missing(cl1.srt) )
            cl1.srt <- cl1$sorting[[K]]
        cl1 <- cl1$clusters[,K]
    } else if ( inherits(cl1, "factor") ) {
        if ( missing(cl1.srt) ) cl1.srt <- levels(cl1)
        cl1 <- as.character(cl1)
    }
    if ( inherits(cl2, "clustering") ) {
        K <- selected(cl2)
        if ( missing(cl2.srt) )
            cl2.srt <- cl2$sorting[[K]]
        cl2 <- cl2$clusters[,K]
    } else if ( inherits(cl2, "factor") ) {
        if ( missing(cl2.srt) ) cl2.srt <- levels(cl2)
        cl2 <- as.character(cl2)
    }
   
  ## check cluster length
  if ( length(cl1) != length(cl2) ) 
      stop("ERROR cluster vectors of different size:", length(cl1),length(cl2))
  
  ## add NA cluster
  if ( sum(is.na(cl1)) )
    cl1[is.na(cl1)] <- na.string
  if ( sum(is.na(cl2)) )
    cl2[is.na(cl2)] <- na.string
  ## also rename empty strings!
  if ( sum(cl1=="") )
    cl1[cl1==""] <- na.string
  if ( sum(cl2=="") )
    cl2[cl2==""] <- na.string
  
  
  ## get requested values
  ## TODO : name p.value result matrices if both are requested!
  do.prich <- sum(alternative=="greater")>0
  do.ppoor <- sum(alternative=="less")>0
  if ( sum(alternative=="two.sided")>0 )
    do.ppoor <- do.prich <- TRUE
  
  ## get clusters
  f1 <- levels(as.factor(cl1))
  f2 <- levels(as.factor(cl2))
    ## TODO: remember if any was not sorted;
    ## and sort those by sortClusters at the end; sort smaller first
    if ( !missing(cl1.srt) ) f1 <-  as.character(cl1.srt)
    else cl1.srt <- as.character(f1)
    if ( !missing(cl2.srt) ) f2 <-  as.character(cl2.srt)
    else cl2.srt <- as.character(f2)
      

  ## result matrices
  overlap <- matrix(NA, nrow=length(f1), ncol=length(f2))
  rownames(overlap) <- f1
  colnames(overlap) <- f2
    jaccard <- percent <- frequency <- ratio <- overlap
  if ( do.prich ) prich <- overlap
  if ( do.ppoor ) ppoor <- overlap
  
  
  for ( i in 1:length(f1) ) { # white and black balls
      
      m <- sum(cl1==f1[i]); # number of white balls
      n <- sum(cl1!=f1[i]); # number of black balls
      N <- m+n # total number of balls
      
      for ( j in 1:length(f2) ) { # balls drawn
          
          q <- sum(cl1==f1[i] & cl2==f2[j]) # white balls drawn
          k <- sum(cl2==f2[j]) # number of balls drawn
          
          overlap[i,j] <- q
          frequency[i,j] <- q/k  # frequency 
          ratio[i,j] <- frequency[i,j] * N/m # frequency ratio
          percent[i,j] <- round(100*frequency[i,j], digits = 2) # TODO: rm this?
          
          ## Calculate cumulative HYPERGEOMETRIC distributions 
          
          ## enrichment:  p-value of finding q or more white balls,
          ## P[X >= x]
          if ( do.prich ) {
              prich[i,j] <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)
              if ( is.na(prich[i,j]) ) 
                prich[i,j]=1
          }
          
          
          ## deprivation: p-value for finding q or less white balls,
          ## P[X <= x]
          if ( do.ppoor ) {
              ppoor[i,j] <- phyper(q=q, m=m, n=n, k=k, lower.tail=TRUE)
              if ( is.na(ppoor[i,j]) )
                ppoor[i,j]=1
          }

          ## Jaccard:
          intersect <- q # overlap
          union <- k+m-q
          jaccard[i,j] <- intersect/union
      }
  }
  
    result <- list(overlap=overlap,   # TODO: rename to counts
                   frequency=frequency,
                   ratio=ratio,
                   percent=percent,
                   jaccard=jaccard)

  ## append p-values
  #p.value <- prich
  if ( do.prich & !do.ppoor ) p.value <- prich
  if ( do.ppoor & !do.prich ) p.value <- ppoor

    ## two.sided!
    ## take smaller p-value, if both are requested!
    if ( do.ppoor &  do.prich ) {
        p.value <- prich
        p.value[ppoor<prich] <- ppoor[ppoor<prich]
        result$sign <- ifelse(ppoor<prich,-1,1)

    }

    result <- append(result, list(p.value=p.value))
    result$alternative <- alternative
    result$varnames=c(vn1,vn2)

    ## TODO: add test number for later bonferroni correction
    ## TODO: add total numbers as result$num.query/total
    ## TODO: align with nomenclature in segmentOverlaps and plotOverlaps
    
    result$num.target <- t(as.matrix(table(cl2)[cl2.srt]))
    result$num.query <- as.matrix(table(cl1)[cl1.srt])

    class(result) <- "clusterOverlaps"
  return(result)
}

#' plot a legend for \code{\link{plotOverlaps}} plots
#' @param p.min significance cutoff, p-values equal or smaller to
#' this cutoff will appear black (one-sided tests) or red/blue
#' (two-sided tests) 
#' @param p.txt p-value cutoff for showing overlap numbers as white instead
#' of black text
#' @param type 1-sided or 2-sided tests
#' @param round round parameter for numeric text
#' @param dir horizontal (1) or vertical (2) orientation
#' @param labels add axis labels
#' @param show.text show log10(p) as text fields (to indicated text p-value
#' cutoff)
#' @param text optional single character to be shown as black/white text
#' in the legend (e.g. \code{text="n"} to indicate the counts used
#' in plots of \code{\link{clusterCluster}} results
## @param side plot side to draw label, use NA to supress
#' @param l number of fields to show in plot
#' @param n number of color shades between \code{p=1} (white)
#' and \code{p >= p.min} (black)
#' @param col color ramp, default are grey values (one-sided) or
#' red (above) and blue (below) for two-sided tests, length of
#' this vector overrules parameter \code{n}
#' @param ... further arguments to \code{\link{plotOverlaps}} 
#' @export
plotOverlapsLegend <- function(p.min=1e-10, p.txt=1e-5, type=1, round=0,
                               show.text=TRUE, text,
                               l=5, n=100, col, dir=1, labels=TRUE, ...) {

    leg <- list()
    pn <- -log10(p.min)
    leg$p.value <- t(t(10^-seq(pn,0,length.out=l)))
    ## use p.values as legend text
    leg$overlap <- t(t(-seq(pn,0,length.out=l)))
    rmz <- FALSE
    if ( !show.text ) { # dirty hack to supress text, TODO: clean up!
        leg$overlap[] <- 0
        rmz <- TRUE
    }
    rownames(leg$overlap) <- rownames(leg$p.value) <- leg$overlap

    ## "negative" p-values indicate two directions, eg. from t-tests
    if ( type==2 ) { # 2-sided

        leg$sign <- cbind(rep(1,l), rep(-1,l))
        leg$p.value <- cbind(leg$p.value, leg$p.value)
        leg$overlap <- cbind(leg$overlap, leg$overlap)

        if ( missing(col) ) {
            docols <- colorRampPalette(c("#FFFFFF","#FF0000"))(n/2)
            upcols <- colorRampPalette(c("#FFFFFF","#0000FF"))(n/2)
            col <- c(rev(docols), upcols)
        }
     } else { # 1-sided
        if ( missing(col) )
            col <- grDevices::gray(seq(1,0,length.out=n))
     }

    if ( !missing(text) ) {
        leg$text <- leg$overlap
        leg$text[] <- text
        leg$overlap <- NULL
    }

    if ( dir==2 )  
        leg <- lapply(leg, t)

    plotOverlaps(leg, p.min=p.min, p.txt=p.txt, round=round,#values="",
                 col=col, rmz=rmz,
                 ylab=NA, axis1.las=1, axis=NULL,xlab=NA,...)
    box()
    if ( labels ) {
        if ( dir==2 ) {
            side <- 1
            at <- 1:2
        } else {
            side <- 2
            at <- 2:1
        }
        side <- ifelse(dir==2, 1, 2)
        
        mtext(expression(log[10](p)), side, par("mgp")[1])
        if ( type==2 )
            axis(dir, at=at, labels=c("lower","higher"),
                 las=ifelse(dir==2, 1, 2))
    }
    x<-leg # silent return of used overlap object
}

#' dotplot
#'
#' plots dotplot profiles of `clusterOverlaps' objects.
#' @export
## TODO: implement different automatic coloring by type, and integrate with
## other functions
dotplot <- function(x,
                    p.min=0.01,
                    dot.sze=c(.3,2), 
                    value=c("overlap","count","statistic","jaccard","text"),
                    n=100, col=viridis::viridis(n), breaks,
                    xpd=FALSE, show.total=FALSE, tot.cex=.8,
                    file, ...) {

    ## values for coloring
    vals <- x[[value]]
    ##if ( lg2 ) {
    ##    vals <- log2(vals)
    ##    if ( missing(mxr) ) mxr <- max(abs(vals))
    ##    vals[vals >  mxr] <-  mxr
    ##    vals[vals < -mxr] <- -mxr
    ##}
    
    ## breaks
    if ( missing(breaks) )
        breaks <- seq(min(vals,na.rm=TRUE), max(vals, na.rm=TRUE),
                     length.out=length(col)+1)

    ## plot empty image
    navals <- vals
    navals[] <- NA
    image_matrix(navals, breaks=breaks, col=col, ...)

    ## dot size scaling by p-value
    p <- x$p.value
    p <- -log10(p)
    p[p>-log10(p.min)] <- -log10(p.min)
    z <- p/-log10(p.min)
    d.sze <- dot.sze[1]+dot.sze[2]*c(t(z))

    ## intersect data to colors
    cols <- col[findInterval(t(vals),
                               seq(min(breaks), max(breaks),
                                   length.out=length(breaks)),
                               all.inside = TRUE)]
    points(x = rep(1:ncol(z), nrow(z)),
           y = rep(nrow(z):1, each= ncol(z)),
           cex=d.sze, pch=19,
           col=cols, xpd=xpd)
    
    toty <- totx <- FALSE
    if (is.logical(show.total)) {
        if (show.total) 
            toty <- totx <- TRUE
    }
    else if (is.character(show.total)) {
        if (show.total == "x") 
            totx <- TRUE
        if (show.total == "y") 
            toty <- TRUE
        if (show.total %in% c("xy", "yx")) 
            toty <- totx <- TRUE
    }
    if (toty) 
        if ("num.query" %in% names(x)) 
            axis(4, at = length(x$num.query):1, mgp=c(0,0,0),
                 labels = format(x$num.query, big.mark=",", trim=TRUE), 
                las = 2, lwd = 0, lwd.ticks = 0, cex.axis=tot.cex)
    if (totx) 
        if ("num.target" %in% names(x)) 
            axis(3, at = 1:length(x$num.target), mgp=c(0,0,0),
                 labels = format(x$num.target, big.mark=",", trim=TRUE), 
                 las = 2, lwd = 0, lwd.ticks = 0, cex.axis=tot.cex)
    
}


#' plot cluster-cluster or segment-segment overlaps
#'
#' Plots the significance distribution of cluster-cluster or segment-segment
#' overlap statistics provided by \code{\link{clusterCluster}},
#' \code{\link{clusterProfile}} or \code{\link{segmentOverlaps}}, where
#' a color gradient is
#' calculated from \code{-log(p)}, and the text shows the overlap numbers,
#' e.g., the number of overlapping features for \code{\link{clusterCluster}},
#' or the Jaccard index or relative intersect values for
#' \code{\link{segmentOverlaps}}. Option \code{text} allows to select
#' which values to plot as text. Only "overlap" is available for a
#' \code{\link{clusterCluster}} results, while for \code{\link{segmentOverlaps}}
#' the Jaccard index ("jaccard"), or the relative intersect size can
#' be shown: "intersect.target" and "intersect.query" are the intersect
#' divided by total target or query length, respectively.
#' Note that two
#' distinct p-value cutoffs can be visualyized: p-values \code{<=p.min}
#' are shown in black, and p-values \code{p.txt} are shown in white instead
#' of black text (thus becoming visible on the black of significant
#' overlaps).
#' @param x a `clusterOverlaps' object returned by
#' \code{\link{clusterCluster}}
#' @param p.min significance cutoff, p-values equal or smaller to
#' this cutoff will appear black (one-sided tests) or red/blue
#' (two-sided tests) 
#' @param p.txt p-value cutoff for showing overlap numbers as white instead
#' of black text
#' @param p.max p-value cutoff for starting the color scale, e.g. 0.05; useful
#' for generally high p-values.
#' @param n number of color shades between \code{p=1} (white)
#' and \code{p >= p.min} (black)
#' @param col color ramp, default are grey values (one-sided) or
#' red (above) and blue (below) for two-sided tests, length of
#' this vector overrules parameter \code{n}
#' @param values selection of text (numeric values) to plot, depends
#' on available data in \code{x}, see "Description"
#' @param type 1 for one-sided or 2 for two-sided tests color scheme,
#' negative values of two-sided tests are reflected in negative p-values
#' @param txt.col two colors used for the plot text, ie., the
#' overlap counts; the second is used if `p<p.txt` as a discrete
#' signficance cutoff
#' @param rmz remove 0 from text values
#' @param short logical, indicating whether to cut higher overlap
#' numbers; currently: division by 1000 and replacement by \code{k}
#' @param scale factor to divide overlap numbers with, useful for
#' low numbers in Jaccard index
#' @param round number of digits to round overlap numbers to (useful
#' for Jaccard index)
#' @param axis integer vector, sets whether x-axis (1,3) and/or
#' y-axis (2,4) are drawn; the column and row names of \code{dat} will
#' be used as tick labels
#' @param show.sig only for overlap lists sorted by
#' \code{\link[segmenTier:sortClusters]{sortClusters}} from package
#' \code{segmenTier}:
#' draws a red line where unsorted non-significant hits start
#' @param show.total show total numbers (counts) of overlapping features
#' on top and right axes
#' @inheritParams image_matrix
#' @param ... arguments to \code{\link{image_matrix}}
## TODO: define as plot.clusterOverlaps method?
## TODO: sort by significance?
## TODO: handle jaccard vs. hypergeo better (see comments)
#' @export
plotOverlaps <- function(x, p.min=0.01, p.txt=p.min*5, p.max, n=100, col,
                         values=c("overlap","count","statistic",
                                  "jaccard","text"),
                         type=1, txt.col = c("black","white"), text.cex=1,
                         rmz=TRUE,
                         short=TRUE, scale=1, round, axis=1:2,
                         show.sig=TRUE, show.total=FALSE, ...) {

    ## set up p-value and colors
    pval <- x$p.value

    ## set p-values above p.max to 1
    if ( !missing(p.max) ) pval[pval>p.max] <- 1

    ## "negative" sign indicated two-sided test
    if ( "sign"%in%names(x) | type==2 ) {

        pval <- pval * x$sign
        sgn <- x$sign[abs(pval)<=p.min] #sign(pval[abs(pval)<=p.min])
        sgn[sgn==0] <- 1
        pval[abs(pval)<=p.min] <- p.min * sgn
        pval <- -log2(abs(pval))*sign(pval)
        
        if ( !missing(col) ) 
            n <- length(col)
        else {
                docols <- colorRampPalette(c("#FFFFFF","#0000FF"))(n/2)
                upcols <- colorRampPalette(c("#FFFFFF","#FF0000"))(n/2)
                col <- unique(c(rev(docols), upcols))
                n <- length(col)
        }
        breaks <- seq(log2(p.min),-log2(p.min),length.out=n+1)
    } else { # NO NEGATIVE P-VALS

        pval[pval<=p.min] <- p.min
        pval <- -log2(pval) # TODO: argument for log-type!
        if ( !missing(col) )
            n <- length(col)
        else
            col <- grDevices::gray(seq(1,0,length.out=n))
        breaks <- seq(0,-log2(p.min),length.out=n+1)
    }
    
    ## set up text (overlap numbers) and text colors
    type <- "" # TODO, add info to input results on the name of the statistics
    for ( i in 1:length(values) ) {
        ## take first available
        if ( values[i]%in%names(x) ) {
            type <- values[i]
            break
        } 
    }

    ## overrule some arguments
    ## TODO: better handle incompatible options round vs. short
    if ( !missing(round) )
        short <- FALSE


    ## NO TEXT:
    if ( type=="" | any(is.na(txt.col)) ) {
        txt <- round <- NA
        image_matrix(pval, breaks=breaks, col=col, axis=axis, ...)
    } else {
        
        ## cat(paste("text values:", type, "\n"))
        
        ## TODO: handle jaccard vs. hypergeo vs. t-test vs. wilcox.test better
        ## and align scale/round/short options
        
        ## shorten large overlap numbers?
        if ( !type%in%c("overlap","count","statistic") )
            short <- FALSE

        txt <- x[[type]]

        if ( type!="text" ) {
            ## parse and process text values
            if ( scale>1 ) txt <- txt*scale
            if ( !missing(round) ) txt <- round(txt,digits=round)
            else round <- NA
        }
             
        ## select color based on p.txt
        if ( rmz) txt[txt=="0"] <- ""
        tcol <- txt
        tcol[] <- txt.col[1]
        tcol[abs(pval) >= -log2(p.txt)] <- txt.col[2]

        ## cut high numbers by 1000: 1321 -> 1.3k
        if ( short ) {
            txt <- x[[type]]
            hg <-txt>1e3
            txt[hg] <- signif(txt[hg]/1e3,2)
            if ( rmz) txt[txt=="0"] <- ""
            txt[hg]  <- paste(txt[hg],"k",sep="")
        }

        ## TODO: allow the following, but "main" needs to be
        ## removed from ... !?
        ##args <- list(...)
        ##if ( "main" %in% names(args) )
        ##    main <- args[["main"]]
        ##else
        ##    main <- paste0("p.min=", p.min, ", p.txt=",p.txt)
        
        image_matrix(pval, breaks=breaks, col=col, axis=axis,
                     text=txt, text.col=tcol, text.cex=text.cex, ...)
    }
    
    ## overlaps sorted by sortClusters indicate where
    ## unsorted non-significant hits start - draw a line: 
    if ( show.sig & "nsig"%in%names(x) )
        if ( x$nsigdir== 2 ) {
            if ( nrow(pval)>x$nsig )
                abline(h=nrow(pval)-x$nsig+.5, col=2, lwd=2)
        } else {
            if ( ncol(pval)>x$nsig )
                abline(v=x$nsig+.5, col=2, lwd=2)
        }

    toty <- totx <- FALSE
    if ( is.logical(show.total) ) {
        if ( show.total )
            toty <- totx <- TRUE
    } else if ( is.character(show.total) ) {
        if ( show.total=="x" ) totx <- TRUE
        if ( show.total=="y" ) toty <- TRUE
        if ( show.total%in%c("xy","yx") ) toty <- totx <- TRUE
    }
    
        
    if ( toty ) 
        if ( "num.query"%in%names(x) )
            axis(4, at=length(x$num.query):1, labels=x$num.query,
                 las=2, lwd=0, lwd.ticks=1, cex.axis=text.cex)
    if ( totx )
        if ( "num.target"%in%names(x) )
            axis(3, at=1:length(x$num.target), labels=x$num.target,
                 las=2, lwd=0, lwd.ticks=1, cex.axis=text.cex)
    ## TODOL place this smarter
    ##if ( toty|totx )
    ##    figlabel("total", region="figure", pos="topright",cex=par("cex"))
    

    ## return plot settings silently
    invisible(list(type=type, short=short, round=round, scale=scale, text=txt,
                   col=col, breaks=breaks, pval=pval))
}

#' transpose cluster overlap object
#'
#' TODO: solve the problem that target/query info is lost, by
#' indicating whether target/query is in row or columns.
#' 
#' 
#' @param x object of class `clusterOverlaps`
#' @seealso \code{\link{clusterCluster}},
#'     \code{\link{clusterProfile}}, \code{\link{clusterAnnotation}}
#' @export
t.clusterOverlaps <- function(x) {
    for ( i in 1:length(x) )
        if ( inherits(x[[i]], "matrix") ) 
            x[[i]] <- t(x[[i]])
    
    ## switch names
    ## TODO: this is a bad solution, since the information
    ## is lost, what the actual query was
    tmpid <- paste(sample(LETTERS, 5, TRUE),collapse="")
    names(x) <- sub("target", tmpid, names(x))
    names(x) <- sub("query", "target", names(x))
    names(x) <- sub(tmpid, "query", names(x))
    
    ## increase transposed counter
    if ( !"parameters"%in%names(x) ) 
        x$parameters <- list()
    if ( !"transposed"%in%names(x$parameters) )
        x$parameters$transposed <- 0
    x$parameters$transposed <- x$parameters$transposed  + 1

    ## significane cutoff direction
    if ( "nsigdir"%in%names(x) ) # direction of sign. cut from sortOverlaps
        x$nsigdir <- ifelse(x$nsigdir==1, 2, 1)
    x
}


#' sorts cluster overlap structure by p-values
#' 
#' Sorts one dimension of a cluster overlap structure by their
#' p-values along the sorting of the other axis.  For each cluster
#' along the pre-sorted (non-selected) axis, the most significant
#' overlaps (\code{p<p.min}) are chosen and moved to the top of the
#' matrix.

## TODO:
## * align axis selection with nomenclature in clusterCluster

#' @param ovl a `clusterOverlaps' object returned by
#'     \code{\link{clusterCluster}}
#' @param p.min significance cutoff during sorting
#' @param cut remove all overlaps without any \code{p<p.min}
#' @param axis axis to sort (2 for y-axis/rows, 1 for x-axis/columns)
#' @param srt sorting vector for rows, if this is passed the
#'     significance sorting is skipped
#' @param sign sort only by p-values of this sign (\code{-1,+1}) for
#'     two-sided tests.
#' @param symmetric indicate whether the overlap matrix is symmetric
#'     with \code{symmetric="upper"} or \code{symmetric="lower"}.
#' @export
sortOverlaps <- function(ovl, axis=2, p.min=.05, cut=FALSE, srt,
                         sign=0, symmetric="no") {

    
    ## handle triangle matrix
    ## NOTE: currently only produced by segmentOverlaps, where
    ## p.values=1 and counts=0 in the lower triangle
    if ( symmetric!="no" ) {

        if ( symmetric=="upper" )
            symm.tri <- lower.tri ## NOTE: assume upper part is filled!
        else if ( symmetric=="lower" )
            symm.tri <- upper.tri
        
        ## copy upper to lower
        pvl <- abs(ovl$p.value)

        ## only consider p-values of the correct sign for two-sided tests
        if ( sign!=0 & "sign"%in%names(ovl) )
            pvl[sign(ovl$sign) != sign(sign)] <- 1
        
        n <- nrow(pvl)
        m <- ncol(pvl)
        if ( n!=m )
            stop("symmetric handling requested for non-symmetric matrix")
        for ( i in 1:length(ovl) ) {
            x <- ovl[[i]]
            if ( inherits(x, "matrix") ) 
                if ( nrow(x)==n & ncol(x)==m )
                    x[symm.tri(x)] <- t(x)[symm.tri(x)]
            ovl[[i]] <- x
        }
    }

    
    ## transpose all, if sorting of x-axis (1) is requested
    if ( axis==1 ) 
        ovl <- t.clusterOverlaps(ovl)

    pvl <- abs(ovl$p.value)

    ## only consider p-values of the correct sign for two-sided tests
    if ( sign!=0 & "sign"%in%names(ovl) )
            pvl[sign(ovl$sign) != sign(sign)] <- 1
        
    
    ## sort by significance
    if ( missing(srt) ) {
        
        cls.srt <- colnames(pvl)

        if ( is.null(colnames(pvl)) )
            stop("missing column names in p-value matrix")
        
        sig.srt <- NULL
        ## first, get highly significant
        for ( cl in cls.srt ) {
            tmp.srt <- order(pvl[,cl], decreasing=FALSE)
            ## cut by p value
            sig.srt <- c(sig.srt,
                         tmp.srt[tmp.srt %in% which(pvl[,cl] < p.min)])
        }
        ## second, sort rest by increasing pval
        rest.srt <- which(!(1:nrow(pvl)) %in% sig.srt)
        rest.srt <- rest.srt[order(apply(pvl[rest.srt,,drop=FALSE],1,max),
                                   decreasing=FALSE)]
        new.srt <- sig.srt[!duplicated(sig.srt)]
        if ( !cut ) new.srt <- c(new.srt, rest.srt)

        ## remember row split between sig and non-sig
        nsig <- sum(!duplicated(sig.srt))
    } else {
        ## used passed sorting!
        new.srt <- srt
        nsig <- NULL

        ## 202307 - tested well in clusterGo.R
        ## warning("custom sorting via `srt` is untested!")
    }
    
    ## resort all matrices in overlap structure (overlap, pvalue, jaccard, ...)
    ## TODO: do this safer, check if everything got sorted?
    n <- nrow(pvl)
    m <- ncol(pvl)
    for ( i in 1:length(ovl) )
        if ( inherits(ovl[[i]], "matrix") ) {
            ## check if matrix is of same dim and if rows are >1
            ## to avoid clash of num.target when new.srt is only of length 1
            if ( nrow(ovl[[i]])==n & nrow(ovl[[i]])>1 )  {
                ovl[[i]] <- ovl[[i]][new.srt,,drop=FALSE]
            }
            ## symmetric case!
            if ( symmetric!="no" & (ncol(ovl[[i]])==m & ncol(ovl[[i]])>1) ) 
                ovl[[i]] <- ovl[[i]][,new.srt,drop=FALSE]
        }
    ## transpose back
    if ( axis==1 )
        ovl <- t.clusterOverlaps(ovl)

    ## symmetric case: set other to NA
    if ( symmetric!="no" ) {

        ## copy upper to lower
        pvl <- abs(ovl$p.value)
        n <- nrow(pvl)
        m <- ncol(pvl)
        if ( n!=m )
            stop("symmetric handling requested for non-symmetric matrix")
        for ( i in 1:length(ovl) ) {
            x <- ovl[[i]]
            replace <- ifelse(names(ovl)[i]=="p.value", 1, 0)
            if ( inherits(x, "matrix") ) 
                if ( nrow(x)==n & ncol(x)==m )
                    x[symm.tri(x)] <- replace
            ovl[[i]] <- x
        }
    }


    ## add cutoff and number of sorted sig
    if ( missing(srt) ) {
        if ( cut )
            ovl$p.min <- p.min
        ovl$p.srt <- p.min
    } else ovl$srt <- srt
    
    ovl$nsig <- nsig
    ovl$nsigdir <- axis # remember direction


    ovl
}



#' Alluvial plot for cluster matrices.
#'
#' Adapted from the code of \cite{alluvial::alluvial}.
#' @param clusters a matrix of different clusterings (columns).
#' @param srt a list of cluster orders for each column in \cite{clusters}.
#' @param cls.col a list of cluster colors used for boxes at start and end
#' of flow; must have the same order as clusters in \code{srt}.
#'@export
clusterFlow <- function (clusters, srt, cls.col, 
                         col = "gray", border = 0, layer,
                         hide = FALSE, alpha = 0.5, gap.width = 0.05,
                         xw = 0.1, cw = 0.1, blocks = TRUE, 
                         axis_labels = NULL, cex = par("cex"),
                         font = par("font"),
                         cex.axis = par("cex.axis")) 
{

    ## TODO: make sure all is correct: srt, cls.col

    ## renumber clusters to reflect sorting!
    ## pad 0 since, sorting is alphabetical
    if ( !missing(srt) ) {
        for ( j in seq_along(srt) ) {
            clusters[,j] <- sprintf("%02d",
                                    match(as.character(clusters[,j]),
                                          rev(srt[[j]])))
         }
    }

    ## flows along clustering columns
    flows <- apply(clusters, 1, paste, collapse=" ")
    ## occurence of each flow pattern
    freq <- table(flows)
    tab <- do.call(rbind, strsplit(names(freq), " "))
    colnames(tab) <- colnames(clusters)
    if ( length(col)>1 )
        cols <- rev(col)[as.numeric(tab[,1])] # todo: colorings for all
    else cols <- col
    
    p <- data.frame(tab, freq = unlist(c(freq)), 
                    col=cols, alpha, border, hide, 
                    stringsAsFactors = FALSE)
    np <- ncol(p) - 5


    n <- nrow(p)
    if (missing(layer)) {
        layer <- 1:n
    }
    p$layer <- layer
    d <- p[, 1:np, drop = FALSE]
    p <- p[, -c(1:np), drop = FALSE]
    p$freq <- with(p, freq/sum(freq))
    col <- col2rgb(p$col, alpha = TRUE)
    if (!identical(alpha, FALSE)) {
        col["alpha", ] <- p$alpha * 256
    }
    p$col <- apply(col, 2, function(x) do.call(rgb, c(as.list(x), 
                                                      maxColorValue = 256)))

    ## really suppress colors!
    p$col[is.na(cols)] <- NA
    
    isch <- sapply(d, is.character)
    d[isch] <- lapply(d[isch], as.factor)
    if (length(blocks) == 1) {
        blocks <- if (!is.na(as.logical(blocks))) {
            rep(blocks, np)
        }
        else if (blocks == "bookends") {
            c(TRUE, rep(FALSE, np - 2), TRUE)
        }
    }
    if (is.null(axis_labels)) {
        axis_labels <- names(d)
    }
    else {
        if (length(axis_labels) != ncol(d)) 
            stop("`axis_labels` should have length ", names(d), 
                ", has ", length(axis_labels))
    }
    getp <- function(i, d, f, w = gap.width) {
        a <- c(i, (1:ncol(d))[-i])

        ## top-down order generated here!
        o <- do.call(order, d[a])

        x <- c(0, cumsum(f[o])) * (1 - w)
        x <- cbind(x[-length(x)], x[-1])
        gap <- cumsum(c(0L, diff(as.numeric(d[o, i])) != 0))
        mx <- max(gap)
        if (mx == 0) 
            mx <- 1
        gap <- gap/mx * w
        (x + gap)[order(o), ]
    }
    dd <- lapply(seq_along(d), getp, d = d, f = p$freq)


    
    rval <- list(endpoints = dd)
    ##op <- par(mar = c(2, 1, 1, 1))
    plot(NULL, type = "n", xlim = c(1 - cw, np + cw), ylim = c(0, 
        1), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", 
        ylab = "", frame = FALSE)
    ind <- which(!p$hide)[rev(order(p[!p$hide, ]$layer))]
    for (i in ind) {
        for (j in 1:(np - 1)) {
            
            graphics::xspline(c(j, j, j + xw, j + 1 - xw, j + 1, j + 1, 
                                j + 1 - xw, j + xw, j) +
                              rep(c(cw, -cw, cw), c(3, 4, 2)),
                              c(dd[[j]][i, c(1, 2, 2)],
                                rev(dd[[j + 1]][i, c(1, 1, 2, 2)]),
                                dd[[j]][i, c(1, 1)]), 
                              shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
                              open = FALSE, 
                              col = p$col[i], border = p$border[i])
        }
    }
    for (j in seq_along(dd)) {
        ax <- lapply(split(dd[[j]], d[, j]), range)
        if (blocks[j]) {
            for (k in seq_along(ax)) {
                colk <- border <- NA ## TODO make nicer
                if ( !missing(cls.col) )
                    colk <- rev(cls.col[[j]])[k]
                else border <- 1
                graphics::rect(j - cw, ax[[k]][1], j + cw, ax[[k]][2],
                               col=colk, border=border)
            }
        }
        else {
            for (i in ind) {
                x <- j + c(-1, 1) * cw
                y <- t(dd[[j]][c(i, i), ])
                w <- xw * (x[2] - x[1])
                
                graphics::xspline(x = c(x[1], x[1], x[1] + w, x[2] - w, 
                                        x[2], x[2], x[2] - w, x[1] + w, x[1]),
                                  y = c(y[c(1, 2, 2), 1],
                                        y[c(2, 2, 1, 1), 2],
                                        y[c(1, 1),  1]),
                                  shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0), 
                  open = FALSE, col = p$col[i], border = p$border[i])
            }
        }
        for (k in seq_along(ax)) {

            nme <- names(ax)[k]
            if ( !missing(srt) )
                nme <- rev(srt[[j]])[k]
            bg <- "white"
            if ( !missing(cls.col) )
                    bg <- rev(cls.col[[j]])[k]
            ## TODO: best color selection?
            shadowtext(j, mean(ax[[k]]), labels = nme, cex = cex,
                       font=font, col=1)
        }
    }
    axis(1, at = rep(c(-cw, cw), ncol(d)) + rep(seq_along(d), 
        each = 2), line = 0.5, col = "white", col.ticks = "black", 
        labels = FALSE)
    axis(1, at = seq_along(d), tick = FALSE, labels = axis_labels, 
        cex.axis = cex.axis)
###par(op)
     return(p)
    invisible(rval)
}


#' Parse a matrix of ID/annotation mappings
#' 
#' Parses a 2D ID/annotation mapping where annotations are a list of
#' terms with a separator. The function returns a TRUE/FALSE table that
#' can be used as input to \code{\link{clusterAnnotation}}.
#' @param got a 2D matrix with IDs in the first column and a list of
#'     annotation terms in the second
#' @param sep separator used for the list of terms in the second
#'     column
#' @export
parseAnnotationList <- function(got, sep=";") {

    terms <- unique(unlist(strsplit(got[,2],sep)))
    terms <- terms[!is.na(terms)]
    mat <- matrix(FALSE, nrow=nrow(got), ncol=length(terms))
    rownames(mat) <- got[,1]
    colnames(mat) <- terms
    
    for ( term in terms )
        mat[grep(term,got[,2]),term] <- TRUE
    mat
}

#' Parse an annotation file (a bidirectional map)
#' 
#' Parses a bidirectional map of feature IDs vs. annotation terms, e.g.
#' the GO annotation file at \url{ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt.gz} and returns a TRUE/FALSE table that
#' can be used as input to \code{\link{clusterAnnotation}}.
#' @param got input table, e.g. a GO annotation 
#' @param idcol column where feature IDs can be found 
#' @param keycol column where annotation terms are found 
#' @param termcol optional column where a readable description of annotation
#' terms is found; will only be used if keycol and termcol are a constant map,
#' ie. each key is always associated with the same term; TODO: terminology?
#' @param rm.empty rm all empty annotations; TODO: why does this occur?
#' @export
parseAnnotation <- function(got, idcol=1, keycol=6, termcol, rm.empty=TRUE) {
    ## EXAMPLE DATA are file
    ## ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt.gz
    ## TODO: allow use of evidence filter
    ## TODO: handle specifier in column 4 ("is downregulated by")
    ## TODO: handle GO category in column 8 
    ## map terms
    gokeys <- unique(got[,keycol]) # no terms; for other columns; TODO: cleaner
    ## readable description - only works for GO terms
    goterms <- unique(got[,termcol])
    if ( length(gokeys)!=length(goterms) )
        goterms <- NULL
    else names(goterms) <- gokeys
    
    ## reduce to unique IDs
    got <- got[!is.na(got[,idcol]),]
    ## generate T/F table
    genes <- unique(got[,idcol])
    gotable <- matrix(FALSE,nrow=length(genes),ncol=length(gokeys))
    rownames(gotable) <- genes
    colnames(gotable) <- gokeys
    for ( gok in gokeys ) 
        gotable[as.character(got[got[,keycol]==gok,idcol]),gok] <- TRUE
    if ( rm.empty ) {
        keep <- apply(gotable,2,any)
        gotable <- gotable[,keep]
        gokeys <- gokeys[keep]
        if ( !is.null(goterms) )
            goterms <- goterms[gokeys]
    }
    list(table=gotable, terms=goterms)
}


#' Cluster annotation enrichment scan
#' 
#' Scans for overlap enrichments of a clustering in a matrix of
#' clusterings, potentially simply a TRUE/FALSE table, eg. indicating
#' annotation with a specific Gene Ontology or other term. This input
#' table can be generated eg. with \code{\link{parseAnnotation}} or
#' \code{\link{parseAnnotationList}}. The function reports
#' tables of all overlap sizes and their enrichment p-values. The
#' function is a wrapper around \code{\link{clusterCluster}}, which
#' performs cumulative hypergeometric distribution tests between
#' two clusterings. The result object can be used with
#' \code{\link{sortOverlaps}} and \code{\link{plotOverlaps}}.
## TODO : re-activate and align usage with contStatTable for continuous data
## TODO 2017: move away from T/F table, use simple ;-sep string of
## annotation terms instead;
## * generalize T/F table usage
## provide pval table for image_matrix
#' @param cls a character or numeric vector with the main clustering
#' @param cls.srt a numeric ar string vector indicating a cluster sorting,
#' allows also to analyze only a subset of clusters in argument \code{cls}
#' @param data a string, numeric or logical matrix with other categorizations
#' of genes. Data rows must correspond to clustering in argument \code{cls}
#' @param p p-value threshold reporting overlaps
#' @param terms optional map of annotation key to descriptions
#' @param replace.terms replace annotation terms by the descriptions
#' in argument \code{terms}
#' @param bin.filter string, indicating bins (categories in \code{data})
#' to be globally omitted; useful eg. for logical data to omit
#' all enrichments with category "FALSE" (indicating deprivement)
#' @param verbose print progress messages
#' @export
clusterAnnotation <- function(cls, data, p=1,
                              cls.srt, terms=NULL, replace.terms=FALSE,
                              bin.filter, verbose=TRUE) {

    ## sorted list of clusters!
    if ( missing(cls.srt) ) {
        cls.srt <- unique(cls)

        ## sort numeric if clusters can be converted!
        if ( all(!is.na(as.numeric(cls.srt[!is.na(cls.srt)]))) ) {
            numeric.sorting <- as.numeric(cls.srt[!is.na(cls.srt)])
            cls.srt <- as.character(sort(numeric.sorting))
        }
      }

    ## filter non-available clusters
    cls.srt <- cls.srt[cls.srt %in% unique(cls)]
        
    cls.num <- length(cls.srt)
    if (verbose) cat(paste("CATEGORICAL DATA ANALYSIS:\n",
                           cls.num, "clusters\t:",
                           paste(cls.srt,collapse=","),"\n"))

    ## TODO : make work with logic table
    

    ## PERFORM HYPERGEOMETRIC DISTRIBUTION TESTS
    ## TODO: generate matrix and dont use rbind; or generate as list
    ## sum(cls.num*length(expected))
    overlap <- NULL ## CONTIGENCY TABLE!
    pvalues <- NULL ## hypergeometric distribution test p.value
    hyp.names <- NULL

    ## ANNOTATION DATA

    ## convert data.frame to matrix
    if ( inherits(data, 'data.frame') ) data <- as.matrix(data)

    ## is it a TRUE/FALSE table?
    logic <- typeof(data)=="logical"
    if ( logic ) { # for TRUE/FALSE table, rm FALSE and rm TRUE from names
        rm.pat <- "_FALSE"
        rp.pat <- "_TRUE"
    }
     
    if (verbose) cat(paste("HYPERGEOMETRIC STATISTICS FOR:\n",
                           ncol(data), "categories\t:\n"))
    ## TODO : parallelize this
    for ( j in 1:ncol(data) ) {

        bins <- data[,j]
        name <- colnames(data)[j]
        expected <- unique(bins)

        if (verbose) cat(paste(j, "-",name, ", ", sep=""))
        
        ## cumulative hypergeometric distribution
        ## TODO: take sum of bonferroni correction factors
        ## correct below in bin.filter
        if ( logic )
            tmp <- clusterCluster(bins, cls, cl1.srt=c("TRUE"),
                                  alternative=c("greater"),cl2.srt=cls.srt)
        else 
            tmp <- clusterCluster(bins, cls,
                                  alternative=c("greater"),cl2.srt=cls.srt)
        
        # get overlap and p.values in original bin order
        if ( !is.null(tmp) ) {
            ovl <- tmp$overlap # contingency table
            pvl <- tmp$p.value
            tnm <- name
            if ( !logic )
                tnm <- paste0(name,': ',rownames(pvl))
            rownames(ovl) <- rownames(pvl) <- tnm
            overlap <- rbind(overlap, ovl)
            pvalues <- rbind(pvalues, pvl)
            # copy name in same size to allow easy cbind below
            hyp.names <- c(hyp.names, rep(name, nrow(ovl)))
          } else {
              warning(paste("WARNING: hypergeo test failed for bin:", name))
          }
    }
    if (verbose) {
        cat(paste("... done\nFINISHED STATISTIC ANALYSIS:\n"))
        cat(paste("\toverlaps  \t", nrow(pvalues), "\n"))
    }
    
    ## CONSTRUCT RESULT TABLES
    ## with descriptions from description and bin list
            
    ## write results to HTML file for each cluster
    total.in.bin <- rowSums(overlap)
    bin.per.genome <- 100*total.in.bin/length(cls)

    ## CLUSTER EVALUATION
    ##  COLLECT FRACTION OF SIGNIFICANT P-VALUES PER CLUSTER
    # fraction of significant (by option p) p-values

    # column for each cluster and one column for total fraction
    psig <- matrix(NA, nrow=(cls.num + 2), ncol=2)
    rownames(psig) <- c(cls.srt, "total", "avg") 
    colnames(psig) <- c("size", "h")


    # total clustered genes 
    total.clustered <- sum(!is.na(cls))
    
    # TOTAL FRACTION by clustering,
    # count fraction of p-values below threshold in ANY of the clusters
    # hypergeometric tests
    hp <- NULL
    for ( row in 1:nrow(pvalues) )
      hp<- c(hp, min(pvalues[row,], na.rm=T))
    hp[is.na(hp)] <- 1 # NA shouldnt happen!?!


    # store fraction of significant
    psig["total",] <- c(total.clustered, sum(hp<p)/length(hp))

    ## EVALUATE EACH CLUSTER
    if (verbose) cat(paste("EVALUATING EACH CLUSTER:\n"))

    hyp.tables <- list()
    for ( j in cls.srt ) {
        if (verbose) cat(paste("cluster", j))

        
        # number of genes in cluster - na.rm required?
        total.in.cluster <- sum(j==cls,na.rm=TRUE)
        psig[j, "size"] <- total.in.cluster
        
        # TODO : get pvalues and overlap matrices from
        # result list instead of storing them??

        if ( j %in% colnames(overlap) ) {
            cluster.per.bin <- 100*overlap[,j]/total.in.bin
            bin.per.cluster <- 100*overlap[,j]/total.in.cluster
            p.values <- pvalues[,j]
            
            # rank by sorting p-value
            rank <- order(order(pvalues[,j],decreasing=FALSE))

            # store fraction of hypergeo p-values below threshold
            #filter <- pvalues[, j] < p
            psig[j, "h"] <- sum(pvalues[, j] < p)/nrow(pvalues)

            if (verbose) cat(paste(" table "))

            ## REPLACE NaN where bin.per.genome is 0
            ## by neutral value 1
            enrichment <- round(bin.per.cluster/bin.per.genome,digits=3)
            enrichment[bin.per.genome==0] <- 1
            
            # generate results table for categorical (hypergeo.)
            # statistics tests of clusters
            hyp.table <- cbind.data.frame(rank, hyp.names, names(p.values),
                                          total.in.bin, overlap[,j],
                                          round(bin.per.genome,digits=1),
                                          round(bin.per.cluster,digits=1),
                                          enrichment,
                                          round(cluster.per.bin,digits=1),
                                          signif(p.values),
                                          stringsAsFactors=FALSE)
        
            colnames(hyp.table) <- c("rank", "category", "bin", 
                                     "bin in genome",  "bin in cluster",
                                     "% of genome", "% of cluster",
                                     "enrichment",
                                     "cluster % of bin", "p-value")
          }
        else { hyp.table<-NULL } # TODO: does this cause downstream problems?


        hyp.tables <- append(hyp.tables, list(hyp.table=hyp.table))
        
        if (verbose) cat(paste(" ... done;\n"))
    }
    names(hyp.tables) <- cls.srt

    if (verbose) cat(paste(" ... done;\n\tSUMMARIZING\n"))

    ## cluster averages in summary table
    only.clusters <- rownames(psig)!="total" 
    psig["avg",] <- rowMeans(t( psig[only.clusters, ] ),na.rm=T)

    if ( !missing(bin.filter) ) 
        rm.pat <- paste0("_",bin.filter,"$")
    
   
    ## FILTER BY P-VALUE, BINS and add DESCRIPTIONS
    ## TODO: add GO classes (MF,CC,BP)
    ## filter significant and order by p-value
    sig <- NULL
    if ( p<1 ) {
        sig <- lapply(hyp.tables, function(x) {
            x <- x[x[,"p-value"] <= p,] #smaller then passed p-value threshold?
            x[order(x[,"p-value"]),]})

        ## filter by bins
        if ( !missing(bin.filter) ) 
            sig <- lapply(sig, function(x)
                x[grep(rm.pat,x[,"bin"], invert=TRUE),] )
        
    
        ## add description of terms
        if ( !is.null(terms) ) # not required for other columns (no keys)
            sig <- lapply(sig, function(x)
                cbind.data.frame(description=terms[x[,"category"]],x))
    }
        
    ##cat(paste("USED FILTER p<", p, "\n"))

    ## process full tables
    ## remove filtered (usually bin.filter==FALSE for a T/F table input)
    if ( !missing(bin.filter) ) 
        pvalues <- pvalues[grep(rm.pat, rownames(pvalues), invert=TRUE),]


    ## remove all where p-value is below threshold
    rm.pvl <- apply(pvalues,1,function(x) !any(x<=p))
    pvalues <- pvalues[!rm.pvl,]
    overlap <- overlap[rownames(pvalues),]

    ## add total counts
    num.target <- t(as.matrix(table(cls)[cls.srt]))
    ## NOTE: only works for logical
    ## TODO: find good alternative
    if ( logic )
        num.query <- as.matrix(apply(data,2,function(x) sum(x)))
    else num.query <- as.matrix(setNames(rep(NA, nrow(pvalues)),
                                         rownames(pvalues)))

    ## filter those reported in p.values/overlap tables
    num.query <- num.query[rownames(pvalues),,drop=FALSE]
    num.target <- num.target[,colnames(pvalues),drop=FALSE]

    ## replace terms by description
    if ( !is.null(terms) & replace.terms ) {
        rownames(pvalues) <- terms[rownames(pvalues)]
        rownames(overlap) <- terms[rownames(overlap)]
        rownames(num.query) <- terms[rownames(num.query)]
    }


    ovl <- list(tables=sig, psig=psig, p.value=pvalues, overlap=overlap,
                num.query=num.query, num.target=num.target)
    class(ovl) <- "clusterOverlaps"
    return(ovl)
           
}

#' analyze internal overlaps in a TRUE/FALSE feature annotation table,
#' where rows are features (eg. genes) and columns are annotation
#' terms. Columns must be named with the annotation terms/IDs.
#' @param x a logical (TRUE/FALSE) feature annotation table
#' @export
annotationOverlap <- function(x) {

    if ( typeof(x)!="logical" )
        stop("input must be a logical (TRUE/FALSE) table")

    cids <- colnames(x)
    cnt <- matrix(0,nrow=length(cids), ncol=length(cids))
    num.query <- num.target <- rep(0, length(cids))
    colnames(cnt) <- rownames(cnt) <- cids
    pvl <- cnt; pvl[] <- 1
    for ( i in 1:length(cids) )
        for ( j in i:length(cids) ) {
            ## TODO: call hypergeo directly here instead of re-routing
            ## via clusterCluster (check if this would be more efficient)
            tmp <- clusterCluster(x[,i], x[,j],
                                  cl1.srt=c("TRUE"), cl2.srt="TRUE",
                                  alternative=c("greater"))
            cnt[i,j] <- cnt[j,i] <-tmp$overlap[1,1]
        pvl[i,j] <- pvl[j,i] <-tmp$p.value[1,1]
        }
    num.query <- as.matrix(apply(x,2,sum))
    num.target <- t(as.matrix(apply(x,2,sum)))
    ovl <- list(overlap=cnt, p.value=pvl,
                num.query=num.query, num.target=num.target)
    ovl
}

#' a wrapper for  \code{gprofiler2}'s \code{gost} function
#'
#' takes a gene clustering, a named vector of gene categories, where
#' the names are gene IDs, and calls the \code{gost} function of the
#' \code{gprofiler2} to do enrichment analysis for all gene
#' annotations available. It then parses the output of the \code{gost}
#' function of the \code{gprofiler2} R package into overlap tables for
#' use with \code{segmenTools}' cluster analysis pipeline.
#' 
#' @param cls a name vector of a gene classification, where the names
#' are gene IDs recognized by \code{gprofiler2} in the respective
#' \code{organism},
#' @param organism the organism ID to which gene IDs of the \code{cls}
#' argument refer to; see \code{?gprofiler2::gost},
#' @param cls.srt an optional sorting of gene classes in \code{cls},
#' @param categories annotation source categories to report,
#' @param terms list of annotation terms to report (over all categories),
#' a name list of terms to report,
#' @param significant only report significant hits, option to \code{gost}.
#' @param verb level of verbosity, 0: no output, 1: progress messages,
#' @param ... further arguments to \code{gprofiler2}'s \code{gost} function.
#' @export
runGost <- function(cls, organism="hsapiens",
                    cls.srt, categories, terms, verb=1,
                    significant=FALSE, ...) {


    ## prepare clustering
    cls.lst <- split(names(cls), f=as.character(cls))

    if ( missing(cls.srt) )
        cls.srt <- names(cls.lst)

    ## sorting and size
    cls.lst <- cls.lst[cls.srt]
    cls.sze <- table(cls)[cls.srt]



    ## get categories
    if ( missing(categories)  ) {
        if ( verb>0 ) cat(paste("loading categories\n"))
        categories <- gprofiler2::get_version_info(organism=organism)$sources
        categories <- names(categories)
    }    

    if ( verb>0 ) cat(paste("calculating enrichments in organism",
                            organism, "\n"))
    gores <- gprofiler2::gost(query=cls.lst, organism = organism,
                              source=categories, significant=FALSE, ...)


    ## parse results into overlap enrichment tables (class "clusterOverlaps")

    if ( verb>0 ) cat(paste("parsing calculated enrichments\n"))
    
    ovll <- list()


    ## USE REQUESTED TERMS!
    ## NOTE: major alternative use case, using requested terms
    ## as categories
    if ( !missing(terms) ) {
        newres <- list()
        if ( verb>0 )
            cat(paste("filtering and restructuring results",
                      "by requested terms\n"))
        for ( i in seq_along(terms) ) {
            trms <- terms[[i]]
            idx <- numeric()
            for ( j in seq_along(trms) ) {
                idx <- c(idx,
                         grep(trms[j], gores$result$term_name))
            }
            idx <- unique(idx)
            trm.results <- gores$result[idx,]
            ## add source to tern name
            trm.results$term_name <- paste0(trm.results$source, ": ",
                                            trm.results$term_name)
            ## replace source by function category
            trm.results$source <- names(terms)[i]
            newres[[i]] <- trm.results
        }
        newres <- do.call(rbind, newres)
        gores$result <- newres
        categories <- names(terms)
    }

    for ( ctgy in categories ) { 
        
        go <- gores$result
        go <- go[go$source==ctgy,]

        if ( verb>0 )
            cat(paste("analyzing category", ctgy, "with",
                      length(unique(go$term_id)), "entries\n"))
    
        ## collect terms
        terms <- go[,c("term_id","term_name","term_size")]
        terms <- terms[!duplicated(terms$term_id),]
        rownames(terms) <- terms$term_id
        

        ## construct overlap matrix
        govl <- matrix(NA, nrow=nrow(terms), ncol=length(cls.lst))
        colnames(govl) <- names(cls.lst)
        rownames(govl) <- terms$term_name
        
        gopvl <- gocnt <- govl
        
        ## generate overlap structure:
        ## fill matrices with overlap p-values and counts, num.query, num.target
        ## TODO: do this more efficiently, one loop and vectors.
        for ( i in 1:nrow(govl) ) 
            for ( j in 1:ncol(govl) ) {
                idx <- which(go$query==colnames(govl)[j] &
                             go$term_id==terms$term_id[i])
                if ( length(idx)>1 ) {
                    stop("too many fields; shouldn't happen")
                    break
                } else if ( length(idx)==0 ) {
                    gopvl[i,j] <- 1
                    gocnt[i,j] <- 0
                } else {
                    gopvl[i,j] <- go$p_value[idx]
                    gocnt[i,j] <- go$intersection_size[idx]
                }
            }
        
        ## construct overlap class
        ovl <- list()
        ovl$p.value <- gopvl
        ovl$count <- gocnt

        ## add annotation term sizes
        ovl$num.query <- as.matrix(terms[,"term_size",drop=FALSE])
        rownames(ovl$num.query) <- terms$term_name

        ## add cluster sizes
        ovl$num.target <- t(as.matrix(cls.sze))
        
        class(ovl) <- "clusterOverlaps"

        ## append to list
        ovll[[ctgy]] <- ovl
    }
    return(ovll)
}


### PLOT cluster-cluster overlaps
#' wrapper for \code{\link[graphics]{image}} plotting a data matrix
#' in the orientation of text display
#' 
#' Wrapper around \code{\link[graphics]{image}} to plot a matrix as
#' it is displayed in R console, i.e. the field \code{dat[1,1]}
#' is at the top left corner. It further allows to plot text
#' into individual fields and have colored axis tick labels.
#' @param z the numeric data matrix to be plotted
#' @param x optional x coordinates, corresponding to columns of z
#' @param y optional y coordinates, corresponding to rows of z
#' @param text a matrix of characters corresponding to \code{dat}
#' which will be plotted on the image
#' @param text.col individual colors for text fields
#' @param text.cex relative font size for text fields
#' @param axis integer vector, sets whether x-axis (1,3) and/or
#' y-axis (2,4) are drawn; the column and row names of \code{dat} will
#' be used as tick labels
#' @param axis1.col invididual colors for x-axis tick labels, length must
#' equal the number of columns of \code{dat}
#' @param axis2.col invididual colors for y-axis tick labels, length must
#' equal the number of rows of \code{dat}
#' @param axis.cex if axis[1|2].col is provided, this sets the tick label size
#' @param axis1.las parameter \code{las} (tick label orientation) for x axis
#' @param axis2.las parameter \code{las} (tick label orientation) for y axis
#' @param breaks breakst to use for color selection.
#' @param cut if TRUE data is cut at min and max of breaks.
#' @param ... further arguments to \code{\link[graphics]{image}}, e.g., col
#' to select colors
#' @export
image_matrix <- function(z, x, y, text, text.col, text.cex=1,
                         axis=1:2, axis.cex=1, cut=FALSE, breaks,
                         axis1.col, axis1.las=2,
                         axis2.col, axis2.las=2, ...) {


    ## TODO: 20240206, clarify strange behaviour for SAAP/LSCC data
    ## ## these behave differently, first is wrong, WHY?
 
    if ( !missing(breaks) )
        cut <- TRUE
    if ( cut ) {
            z[z<min(breaks)] <- min(breaks)
            z[z>max(breaks)] <- max(breaks)
    }

    ## reverse columns and transpose
    if ( nrow(z)>1 )
        imgdat <- t(apply(z, 2, rev))
    else
        imgdat <- t(z)
    axis1.numeric <- !missing(x)
    axis2.numeric <- !missing(y)
    if ( missing(x) ) x <- 1:ncol(z)
    if ( missing(y) ) y <- 1:nrow(z)

    ## call image with transformed data, pass breaks and all arguments
    image(x=x, y=y, z=imgdat, axes=FALSE, breaks=breaks, ...)

    ## add text
    if ( !missing(text) ) {
        if ( missing(text.col) )
            text.col <- rep(1, length(c(text)))
        text(x=rep(1:ncol(z),nrow(z)), y=rep(nrow(z):1,each=ncol(z)),
             paste(t(text)),col=t(text.col),cex=text.cex)
    }

    
    ## add axes
    ## TODO : handle axes=FALSE
    if ( !missing(axis) ) {
        if ( 1 %in% axis ) 
            if ( !missing(axis1.col) ) # colored ticks
                for ( i in 1:ncol(z) )
                    axis(1, at=i,colnames(z)[i],
                         col.axis=axis1.col[i], col=axis1.col[i],
                         las=axis1.las, cex.axis=axis.cex, lwd=2)
            else if ( !axis1.numeric )
                axis(1, at=1:ncol(z), labels=colnames(z), las=axis1.las,
                     cex.axis=axis.cex)
            else axis(1)               
        if ( 3 %in% axis ) 
            if ( !missing(axis1.col) ) # colored ticks
                for ( i in 1:ncol(z) )
                    axis(3, at=i,colnames(z)[i],
                         col.axis=axis1.col[i], col=axis1.col[i],
                         las=axis1.las, cex.axis=axis.cex, lwd=2)
            else if ( !axis1.numeric )
                axis(3, at=1:ncol(z), labels=colnames(z), las=axis1.las,
                     cex.axis=axis.cex)
            else axis(3)               
        if ( 2 %in% axis )
            if ( !missing(axis2.col) ) # colored ticks
                for ( i in 1:nrow(z) )
                        axis(2, at=nrow(z)-i+1, rownames(z)[i],
                             col.axis=axis2.col[i], col=axis2.col[i],
                             las=axis2.las, cex.axis=axis.cex, lwd=2)
            else if ( !axis2.numeric )
                axis(2, at=nrow(z):1, rownames(z),las=axis2.las,
                     cex.axis=axis.cex)        
            else axis(2)               
        if ( 4 %in% axis )
            if ( !missing(axis2.col) ) # colored ticks
                for ( i in 1:nrow(z) )
                        axis(4, at=nrow(z)-i+1, rownames(z)[i],
                             col.axis=axis2.col[i], col=axis2.col[i],
                             las=axis2.las, cex.axis=axis.cex, lwd=2)
            else if ( !axis2.numeric )
                axis(4, at=nrow(z):1, rownames(z),las=axis2.las,
                     cex.axis=axis.cex)
            else axis(4)
    }
} 

### CLUSTERING OBJECT UTILS

#' get clustering ID
#'
#' get the column name or index of the selected
#' clustering from a `clustering' object from segmenTier's
#' clustering wrappers
#' \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}) and
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}).
#' The retuned integer or string allows to interface data for this clustering
#' from the `clustering' object.
#' @param cset  a structure of class 'clustering' as returned by
#' segmenTier's \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param K cluster number; if missing the `selected' clustering will
#' be taken
#' @param name logical, if TRUE the name of of the selected clustering
#' will be returned, if FALSE its index number
#' @export
selected <- function(cset, K, name=TRUE) {
    if ( missing(K) )
      K <- cset$selected
    if ( is.numeric(K) )
      kCol <- paste0("K:", K)
    else kCol <- K
    ## add K
    if ( length(grep("^K:",K))==0)
        kCol <- paste0("K:", K)

    # return
    if ( name )
        return(kCol)
    else
        return(which(colnames(cset$clusters)==kCol))
}

#' get color vector for clustered features
#'
#' Retrieves a color vector for clustered features from
#' segmenTier's clustering object. Default is to return colors for the
#' pre-selected clustering via option \code{K} to function
#' \code{\link{selected}}. Not to mix up with segmenTier's
#' function \code{\link[segmenTier:colorClusters]{colorClusters}} which
#' assigns a cluster coloring scheme.
#' @param cset a structure of class 'clustering' as returned by
#' segmenTier's \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param expand logical indicating whether colors should be expanded
#' to all data
#' @param ... arguments to cluster selection with function
#' \code{\link{selected}}, e.g., cluster number K 
#' @export
clusterColors <- function(cset, expand=TRUE, ...) {
    K <- selected(cset, ...)
    if ( expand )
        return(cset$colors[[K]][as.character(cset$clusters[,K])])
    else return(cset$colors[[K]][as.character(cset$sorting[[K]])])
}


#' Cluster a processed time-series with k-means or flowClust
#' 
#' A wrapper for clustering a time-series object \code{tset} provided by
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}},
#' where specifically the DFT of a time-series and requested data
#' transformation were calculated. The clustering is performed
#' on the \code{tset$dat} matrix by \code{\link[stats:kmeans]{kmeans}} or
#' model-based by packages \pkg{flowClust} & \pkg{flowMerge}.
#' 
#' This function attempts to combine the previous separate clustering
#' wrappers from package \code{segmenTier}, \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}} and k-mean's based \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}, the latter of which is used for
#' the segmenTier algorithm. Please see the corresponding help files for
#' details on the clustering parameters in argument \code{parameters}.
#' @param tset a timeseries processed by \code{\link[segmenTier:processTimeseries]{processTimeseries}}
#' @param K selected cluster numbers, the argument \code{centers}
#' of \code{\link[stats:kmeans]{kmeans}}
#' @param method string specifying the clustering algorithm to use, currently
#' "kmeans" and "flowClust" are implemented
#' @param parameters named vector of parameters for clustering algorithms,
#' currently ALL required parameters MUST be specified
#' @param selected a pre-selected cluster number  which is then
#' used as a start clustering for \code{flowMerge} (if option
#' \code{merge==TRUE})
#' @param nui.thresh threshold correlation of a data point to a cluster
#' center; if below the data point will be added to nuissance cluster 0
#' @param ncpu number of cores available for parallel mode of
#' \pkg{flowClust}. NOTE: parallel mode of
#' \code{\link[flowClust:flowClust]{flowClust}} is often non-functional.
#' Alternatively, you can set \code{options(mc.cores=ncpu)} directly.
#' @param verb level of verbosity, 0: no output, 1: progress messages
#' @param ... further parameters to \code{flowClust} or \code{\link[stats:kmeans]{kmeans}}
#'@export
clusterTimeseries2 <- function(tset, K=16, method="flowClust", selected, 
                               parameters=c(iter.max=100000, nstart=100,
                                            B=500, tol=1e-5, lambda=1,
                                            nu=4, nu.est=0, trans=1,
                                            merge=FALSE, randomStart=0),
                               nui.thresh=-Inf,  ncpu=1, verb=1,...) {

    if ( method=="flowClust" ) {
        # dont ## suppress parallel mode!
        # ncpu=1
        oldcpu <- unlist(options("cores"))
        oldcpu2 <- unlist(options("mc.cores"))
        options(cores=ncpu)
        options(mc.cores=ncpu)
    }
    
    ## get time series data
    id <- tset$id
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    N <- nrow(dat)

    ## enought distinct values?
    ## TODO: issue segment based on low-filter
    ## OOR: postprocessing - extend segments into low levels?
    warn <- NULL
    if ( sum(!rm.vals)<10 ) 
        warn <- "not enough data"
    else if ( sum(!duplicated(dat[!rm.vals,]))<2 ) 
        warn <- "not enough data diversity"
    if ( !is.null(warn) ) {
        warning(warn)
        return(NULL)
    }
    
    ## CLUSTERING
    ## stored data
    clusters <- matrix(NA, nrow=nrow(dat), ncol=length(K))
    rownames(clusters) <- rownames(dat)
    results <- centers <- Pci <- Ccc <- rep(list(NA), length(K))

    ## BIC/AIC 
    bic <- rep(NA, length(K))
    names(bic) <- as.character(K)
    icl <- aic <- bic
    
    if ( verb>0 ) {
        cat(paste("Timeseries N\t",N,"\n",sep=""))
        cat(paste("Used datapoints\t",sum(!rm.vals),"\n",sep=""))
    }

    ## run flowClust over all K (multithreaded)
    if ( method=="flowClust" ) {
        ## TODO: defaults if missing
        B <- parameters["B"]
        tol <- parameters["tol"]
        lambda <- parameters["lambda"]
        nu <- parameters["nu"]
        nu.est <- parameters["nu.est"]
        trans <- parameters["trans"]
        randomStart <- parameters["randomStart"]
        
        fcls <- flowClust::flowClust(dat[!rm.vals,], K=K, B=B, tol=tol,
                                    lambda=lambda, randomStart=randomStart,
                                    nu=nu, nu.est=nu.est, trans=trans,...)
    }

    ## loop over K: run k-means or retrieve values from
    ## pre-calculated flowClust
    
    usedk <- K
    for ( k in 1:length(K) ) {
        
        ## get cluster number K
        Kused <- min(c(K[k],sum(!duplicated(dat[!rm.vals,]))))
        
        if ( verb>0 )
            cat(paste("Clusters K\t", Kused, "\n",sep=""))
        
        ## cluster
        if ( method=="kmeans" ) {

            iter.max <- parameters["iter.max"]
            nstart <- parameters["nstart"]
            
            km <- stats::kmeans(dat[!rm.vals,], Kused, iter.max=iter.max,
                                nstart=nstart, algorithm="Hartigan-Wong", ...)
            ## use alternative algo if this error occured
            if (km$ifault==4) {
                km <- stats::kmeans(dat[!rm.vals,], Kused,
                                    iter.max=iter.max,nstart=nstart,
                                    algorithm="MacQueen")
                warn <- "quick-transfer error in kmeans algorithm Hartigan-Wong, taking MacQueen"
                warning(warn)
            }
        
            ## prepare cluster sequence
            seq <- rep(0, N) ## init. to nuissance cluster 0
            seq[!rm.vals] <- km$cluster
            
            ## store which K was used, the clustering and cluster centers
            usedk[k] <- Kused
            clusters[,k] <- seq
            centers[[k]] <- km$centers
            
            ## calculate BIC/AIC
            bic[k] <- stats::BIC(km)
            aic[k] <- stats::AIC(km)

            ## store full results from kmeans
            results[[k]] <- km

        } else if ( method=="flowClust" ) {

            ## get flowClust object
            if ( length(fcls) > 1 ) fc <- fcls[[k]]
            else fc <- fcls
      

            ## prepare cluster sequence
            seq <- rep(0, N) ## init. to nuissance cluster 0
            seq[!rm.vals] <- flowClust::Map(fc,rm.outliers=FALSE)
            
            ## store which K was used, the clustering and cluster centers
            usedk[k] <- fc@K
            clusters[,k] <- seq

            ## get cluster centers!
            x <- fc@mu
            rownames(x) <- 1:nrow(x)
            colnames(x) <- colnames(dat)
            centers[[k]] <- x

            ## BIC/ICL
            bic[k] <- fc@BIC
            icl[k] <- fc@ICL

            ## store full result object from flowCluster
            results[[k]] <- fc
        }

        ## C(c,c) - cluster X cluster cross-correlation matrix
        cr <- stats::cor(t(centers[[k]]))

        Ccc[[k]] <- cr
        
        ## P(c,i) - position X cluster correlation
        P <- matrix(NA,nrow=N,ncol=Kused)
        P[!rm.vals,] <- segmenTier::clusterCor_c(dat[!rm.vals,], centers[[k]])

        Pci[[k]] <- P

    }

    ## re-assign by correlation threshold
    ## NOTE: this only affects scoring function ccor
    for ( k in 1:ncol(clusters) ) {
        cls <- clusters[,k]
        for ( p in 1:nrow(Pci[[k]]) )
            if ( !any(Pci[[k]][p,] > nui.thresh, na.rm=TRUE) )
                cls[p] <- 0
        clusters[,k] <- cls
    }           
    
    ## count duplicate K
    if ( any(duplicated(K)) ) {
        sel <- paste(K,".1",sep="")
        cnt <- 2
        while( sum(duplicated(sel)) ) {
            sel[duplicated(sel)] <- sub("\\..*",paste0(".",cnt),
                                        sel[duplicated(sel)])
            cnt <- cnt+1
        }
        K <- sub("\\.1$","",sel)
    }
    ## name all results by K, will be used!
    colnames(clusters) <- names(centers) <- names(results) <- 
        names(Pci) <- names(Ccc) <- names(usedk) <- paste0("K:",K) 

    ## max BIC and ICL
    max.bic <- max(bic, na.rm=T)
    max.clb <- K[which(bic==max.bic)[1]]
    if ( any(!is.na(aic)) )
        max.aic <- max(aic, na.rm=T)
    else max.aic <- NA
    max.cla <- K[which(aic==max.aic)[1]]
    if ( any(!is.na(icl)) ) {
        max.icl <- max(icl, na.rm=T)
        max.cli <- K[which(icl==max.icl)[1]]
    } else max.icl <- max.cli <- NA
    ## best K selection
    ## use K with max BIC
    selected <- max.clb
    
    ## clustering data set for use in segmentCluster.batch 
    cset <- list(clusters=clusters, 
                 N=sum(!rm.vals), # number of clustered data
                 M=ncol(dat), # dimension of clustered data
                 centers=centers, Pci=Pci, Ccc=Ccc,
                 K=K, usedk=usedk, selected=selected,
                 bic=bic, aic=aic, icl=icl,
                 max.clb=max.clb, max.cla=max.cla, max.cli=max.cli,
                 warn=warn, ids=colnames(clusters),
                 tsid=rep(id,ncol(clusters)),
                 method=method, results=results)
    class(cset) <- "clustering"

    ## add cluster sorting & colors
    cset <- segmenTier::colorClusters(cset)


    ##
    if ( method=="flowClust" ) {
        options(cores=oldcpu)
        options(mc.cores=oldcpu2)
    }

    
    ## silent return
    tmp <- cset
}

#' re-cluster clustering by \code{\link[stats:kmeans]{kmeans}}
#'
#' Use cluster centers from an initial clustering to initialize
#' \code{\link[stats:kmeans]{kmeans}}. This is still experimental,
#' and used to re-associated data rows to cluster centers from
#' a best clustering found by
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}.
#' While the latter clustering works best to extract specific time-courses
#' from the data set, it often comes with a high fraction of badly
#' associated individual data sets. Re-clustering with
#' \code{\link[stats:kmeans]{kmeans}} seems to clean this up, e.g., the
#' phase distributions of re-clustered clusterings are often tighter.
#' TODO: allow to generate cluster centers from novel data, to
#' account for different/more data then during clustering!
#' @param tset the `timeseries' object from segmenTier's
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}} used
#' for initial clustering
#' @param cset the `clustering' object from segmenTier's 
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}
#' @param k colum name or index of the clustering that should be
#' re-clustered; defaults to the pre-selected clustering if missing
#' @param select use the re-clustered clustering as the new pre-selected
#' clustering
#' @param ... parameters to \code{\link[stats:kmeans]{kmeans}}
#' @export
reCluster <- function(tset, cset, k, select=TRUE, ...) {

    if ( missing(k) )
      k <- selected(cset, name=TRUE)

    recls <- tryCatch(stats::kmeans(tset$dat[!tset$rm.vals, ],
                                    centers = cset$centers[[k]], 
                                    algorithm ="Hartigan-Wong", ...),
                      error = function(e) return(list(ifault=4))
    )    
	## use alternative algo if this error occured
    warn <- NULL
    if (recls$ifault==4) {
        recls <- stats::kmeans(tset$dat[!tset$rm.vals,],
                               centers=cset$centers[[k]],
                               algorithm="MacQueen", ...)
        warn <- "quick-transfer error in kmeans algorithm Hartigan-Wong, taking MacQueen"
        warning(warn)
    }
    
    cls <- rep(0, nrow(cset$clusters))
    cls[!tset$rm.vals] <- recls$cluster
    cset$clusters <- cbind(cset$clusters,cls)
    ## copy existing sorting and coloring
    cset$sorting <- append(cset$sorting, cset$sorting[k])
    cset$colors <- append(cset$colors, cset$colors[k])

    K <- nrow(cset$centers[[k]])
    N <- length(cls)
    ## add cluster data
    ## cluster centers
    cset$centers <- append(cset$centers, list(recls$centers))
    ## C(c,c) - cluster X cluster cross-correlation matrix
    cset$Ccc <- append(cset$Ccc, list(stats::cor(t(recls$centers))))
    ## P(c,i) - position X cluster correlation
    P <- matrix(NA,nrow=N, ncol=K)
    P[!tset$rm.vals,] <- segmenTier::clusterCor_c(tset$dat[!tset$rm.vals,],
                                                  recls$centers)
    cset$Pci <- append(cset$Pci, list(P))
    ## warning message from kmeans
    if ( !is.null(warn) )
        if ( "warn" %in% names(cset) )
            cset$warn <- c(cset$warn, warn)
        else cset$warn <- warn

    ## TODO: BIC, ICL? calculate or add NA

    cset$K <- c(cset$K, paste0(K,"_re")) # character: pre-selected K.duplicates
    cset$usedk <- c(cset$usedk, K) # numeric: actual K used

    ## add name
    newKcol <- paste0(k,"_re")
    idx <- ncol(cset$clusters)
    colnames(cset$clusters)[idx] <- names(cset$sorting)[idx] <-
      names(cset$colors)[idx] <- names(cset$centers)[idx] <-
        names(cset$Ccc)[idx] <- names(cset$Pci)[idx] <-
        names(cset$usedk)[idx] <- newKcol



    cset$reclustered <- newKcol
    if ( select )
        cset$selected <- newKcol
    cset

}

## TODO: use flowMerge to merge selected or best BIC clustering
## from cset (clusterTimeseries2)
#' using \pkg{flowMerge} to merge clusterings 
#' @param tset the `timeseries' object from segmenTier's
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}} used
#' for initial clustering
#' @param cset the `clustering' object from segmenTier's 
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}
#' @param selected optional pre-selection of clustering to merge; if missing
#' the pre-selected clustering (usually max. BIC)  in \code{cset} will be used
#'@export
mergeCluster <- function(tset, cset, selected) {

    if ( missing(selected) )
        selected <- selected(cset)

    ## get original data
    dat <- tset$dat
    rm.vals <- tset$rm.vals
    clsDat <- dat[!rm.vals,]

    ## MERGE CLUSTERS, starting from best BIC by flowMerge
    mrg.orig <- mrg.cl <- mrg.id <-  obj <- NULL

    ## get start clustering
    K <- cset$K
    fc <- cset$results[[selected]]


    ## prepare results
    mcls <- rep(0, nrow(dat))
    mrg.id <- mrg.cl <- "NA"

    ## initiate flowMerge object
    obj <- try(flowMerge::flowObj(fc, flowCore::flowFrame(clsDat)))

    ## start flowMerge
    if ( !inherits(obj, "try-error") ) {
        mrg <- try(flowMerge::merge(obj))
        if ( !inherits(mrg, "try-error") ) {
            mrg.cl <- flowMerge::fitPiecewiseLinreg(mrg)
            obj <- mrg[[mrg.cl]]
            mcls[!rm.vals] <- flowClust::Map(obj, rm.outliers=FALSE)
            mrg.id <- paste0(selected,"m",mrg.cl) # merged K, column name

            ## NOTE/TODO: rownames of obj@mu contain merge info
            ## where else? use to get consistent names?
            ## perhaps in mtree structure?
            ## see plot(obj@mtree)
            centers <- obj@mu
            rownames(centers) <- 1:nrow(centers) # new cluster labels=row order
            colnames(centers) <- colnames(clsDat)
            
            ## add clusters to matrix
            cset$clusters <- cbind(cset$clusters, mcls)

            ## add cluster sorting, colors, centers, Ccc, Pci
            cset$centers <- append(cset$centers, list(centers))
            ## C(c,c) - cluster X cluster cross-correlation matrix
            cset$Ccc <- append(cset$Ccc, list(stats::cor(t(centers))))
            ## P(c,i) - position X cluster correlation
            P <- matrix(NA,nrow=nrow(dat), ncol=mrg.cl)
            P[!tset$rm.vals,] <-
                segmenTier::clusterCor_c(tset$dat[!tset$rm.vals,],
                                         centers)
            cset$Pci <- append(cset$Pci, list(P))

            ## add K and usedk
            cset$usedk <- c(cset$usedk, mrg.cl)
            cset$K <- c(cset$K, mrg.id)
            
            cset$merged <- mrg.id
            cset$merged.K <- mrg.cl
            cset$merged.origK <- cset$usedk[selected(cset)]

            ## add name
            idx <- ncol(cset$clusters)
            colnames(cset$clusters)[idx] <- names(cset$centers)[idx] <-
              names(cset$Ccc)[idx] <- names(cset$Pci)[idx] <-
                names(cset$usedk)[idx] <- mrg.id

            ## add sorting and colors!
            ## tmp: remove old sorting and colors and redo all
            ## TODO: just do this for added merged cluster!
            cset$sorting <- NULL
            cset$colors <- NULL
            cset <- segmenTier::colorClusters(cset)

            
            cset$results <- append(cset$results,obj)
        } else cat(paste("error: flowObj failed\n")) 
    } else cat(paste("error: merge failed\n")) 
    return(cset)
}


## TODO; finish implementation
## sort clusters by time series phase
##
## takes a time-series and a clustering object,
## calculates phases for each cluster and sorts clusters
## by phase
phasesortClusters <- function(ts, cls, phase, cycles) {

    
    orig <- NULL
    if ( inherits(ts, "timeseries") )
        ts <- ts$ts
    if ( is.vector(cls) )
        cls <- matrix(cls,ncol=1)
    if ( inherits(cls, "clustering") ) {
        ## store if cls is a clustering set
        orig <- cls
        cls <- cls$clusters
    }

    ## calculate segment phase
    if ( missing(phase) )
        phase <- calculatePhase(ts, cycles=cycles)

    for ( k in 1:ncol(cls) ) {
        cl <- cls[,k]
        
        cl.srt <- unique(cl)
        names(cl.srt) <- cl.srt
        
        phase <- sapply(cl.srt, function(x) phaseDist(phase[cl==x,]))
        #names(sort(phase[grep("mean",rownames(phase)),]))
        colnames(phase) <- cl.srt
        cl.srt <- cl.srt[order(phase[1,])]
        
    }
    
    ## sort by phase!
    phases <- calculatePhase(ts$ts, cycles=3)
    nset <- cls
    for ( i in 1:length(cls$sorting) ) {
        cls <- as.character(cls$clusters[,i])
        cls.srt <- cls$sorting[[i]]
        phs <- sapply(cls.srt, function(x) phaseDist(phases[cls==x]))
        colnames(phs) <- cls.srt
        nset$sorting[[i]] <- cls.srt[order(phs[1,])]
    }
    nset <- relabelClusters(nset) # re-label by sorting
    cls <- nset
}

#' relabels cluster labels by their sorting
#'
#' relabels all cluster labels in a `clustering' object (as returned
#' by \code{\link{clusterTimeseries2}},
#' \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}} and
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}})
#' by the cluster sorting (in \code{cls$sorting}), such that the new
#' numeric cluster labels reflect this sorting.
#' @param cls the `clustering' object from segmenTier's 
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}
## TODO: alternatively just add a list of cluster labels
## to the clustering object and use this in plot functions!
#' @export
relabelClusters <- function(cls) {
    if ( !inherits(cls, "clustering") )
        stop("function requires class 'clustering',",
             " as returned by clusterTimeseries2")
    for ( i in 1:length(cls$sorting) ) {
        k <- names(cls$sorting)[i]
        srt <- 1:length(cls$sorting[[i]])
        names(srt) <- cls$sorting[[i]]

        ## re-order and re-name matrices
        cls$centers[[k]] <- cls$centers[[k]][names(srt),]
        rownames(cls$centers[[k]]) <- srt
        ## TODO: upstream - provide colnames for Pci!!!
        cls$Pci[[k]] <- cls$Pci[[k]][,as.numeric(names(srt))]
        colnames(cls$Pci[[k]]) <- srt
        cls$Ccc[[k]] <- cls$Ccc[[k]][names(srt),names(srt)]
        rownames(cls$Ccc[[k]]) <- colnames(cls$Ccc[[k]]) <- srt

        ## add 0-cluster and re-write clusters, sorting and colors
        if ( !0 %in% srt )
            srt <- c(srt,"0"=0)

        cls$clusters[,k] <- srt[as.character(cls$clusters[,k])]

        cls$sorting[[k]] <- as.character(srt[cls$sorting[[k]]])
        names(cls$colors[[k]]) <- srt[names(cls$colors[[k]])]
    }
    cls
}

### PLOT CLUSTERED TIME-SERIES

#' plot polar coordinates
#'
#' Plots the components of a Discrete Fourier Transform (DFT)
#' as polar coordinates (Re and Im of the complex numbers in the DFT).
#' Arguments \code{dft} and \code{col} can be segmenTier timeseries
#' and clustering objects.
#' @param dft the Fourier transform of a time series as returned
#' by \code{t(mvfft(t(timeseries)))}, or alternatively, a `timeseries' object
#' from segmenTier's
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}} when
#' run with (\code{use.fft=TRUE})
#' @param cycles the number of cycles (index of non-DC DFT component)
#' to be plotted
#' @param radius radius of the polar plot circle as a fraction
#' of data to be contained within the radius (smaller amplitude)
#' @param col a color vector for the rows in argument \code{dft} or
#' alternatively,  `clustering' object as returned by
#' segmenTier's \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' with coloring information
#' @param lambda parameter for Box-Cox transformation of DFT data; has no
#' effect for \code{lambda==1}
#' @param bc type of Box-Cox transformation (\code{if lambda!=1});
#' "component": separate transformation of real and imaginary parts of
#' the DFT; "amplitude": Box-Cox transformation of the amplitude
#' @param ... arguments to the base \code{\link[graphics:plot]{plot}} 
#' and/or \code{\link[graphics:points]{points}} functions
#' @export
plotDFT <- function(dft, col, cycles=3, radius=.9, lambda=1, bc="component", ...) {

    ## dft
    ## can be a segmenTier timeseries object
    if ( inherits(dft, "timeseries") )
        dft <- dft$dft
 
    ## colors
    if ( missing(col) )
        col <- rep("#00000077",nrow(dft))
    else if ( inherits(col, "clustering") )
        col <- clusterColors(col, expand=TRUE)
    else if ( length(col)==1 )
        col <- rep(col,nrow(dft))

    ## split into Re/Im parts
    re <- Re(dft)
    colnames(re) <- paste("Re_",colnames(re),sep="")
    im <- Im(dft)
    colnames(im) <- paste("Im_",colnames(im),sep="")
    ## filter 0 imaginary components: DC and Nyquist!
    im <- im[,apply(im,2,function(x) any(x!=0,na.rm=TRUE))]
    dat <- cbind(re,im)

    
    ## temporary for exploration of Box-Cox
    ## box-cox trafo for negative values (Bickel and Doksum 1981)
    ## as used in flowclust
    bccmp <- function(x,lambda) (sign(x)*abs(x)^lambda-1)/lambda
    ## amplitude box-cox trafo for complex polar coordinates 
    bcdft <- function(x, lambda) {
        if ( inherits(x, "matrix") )
            return(apply(x,2, bcdft, lambda))
        ## Box-Cox transform amplitude
        y <- bccmp(abs(x), lambda)
        ## amplitude scaling factor
        sf <- (y-min(y,na.rm=T))/abs(x)
        x*sf
    }

    ## angles for drawing circle
    theta <- seq(0, 2 * pi, length = 200)
    ## radii - before box cox
    circles <- rep(list(NA), length(cycles))
    names(circles) <- cycles
    for ( cycle in cycles ) {
        rd <- quantile(abs(dft[!is.na(col),cycle+1]),
                       probs=radius, na.rm=TRUE)
        circle <- xy.coords(x = rd * cos(theta),
                            y = rd * sin(theta))
        circles[[as.character(cycle)]] <- circle
    }
    
    ## Box-Cox Transformation
    ori.line <- 0
    if ( lambda!=1 ) {
        if ( bc == "component" ) {
            ## Box-Cox for components
            for ( i in 1:ncol(dat) ) 
                dat[,i] <- bccmp(dat[,i], lambda=lambda)
            ## re-construct DFT
            for ( cycle in cycles ) {
                dft[,cycle+1] <-
                    complex(real=dat[,paste0("Re_",cycle)],
                            imaginary=dat[,paste0("Im_",cycle)])
                ## box-cox for circles
                circ <- circles[[as.character(cycle)]]
                circ$x <- bccmp(circ$x, lambda=lambda)
                circ$y <- bccmp(circ$y, lambda=lambda)
                circles[[as.character(cycle)]] <- circ                
            }
            ori.line <- bccmp(0,lambda)
        }
        if ( bc == "amplitude" ) # note: box-cox calculated below
            dft <- bcdft(dft,lambda)
    }
    
    for ( cycle in cycles ) {
        plot(dft[,cycle+1],col=NA,
             xlab=bquote("Real(X"[.(cycle)]~")"),
             ylab=bquote("Imaginary(X"[.(cycle)]~")"),axes=FALSE,
             main=paste0(cycle, ifelse(cycle==1," cycle"," cycles")), ...)
        axis(1);axis(2)
        abline(v=ori.line,col=1,lwd=1)
        abline(h=ori.line,col=1,lwd=1)
        points(dft[,cycle+1], col=col, ...)

        ## draw circle
        ## calculated above for component Box-Cox
        if ( bc=="component" | lambda==1)
            lines(circles[[as.character(cycle)]])
        else if ( bc=="amplitude" ) {
            rd <- quantile(abs(dft[!is.na(col),cycle+1]),
                           probs=radius, na.rm=TRUE)
            lines(x = rd * cos(theta) + ori.line,
                  y = rd * sin(theta) + ori.line)
        }       
     }
    ## silent return
    invisible(list(bc=bc, lambda=lambda))
}


#' calculates cluster averages
#' 
#' calculates average values and distributions for each cluster and
#' time point of a time series
#' @param ts a matrix of time series, with time points in columns
#' @param cls a clustering of the time series, \code{length(cls)} must
#'     equal \code{nrow(cls)}
#' @param cls.srt optional sorting of the clusters
#' @param avg a function (or the name of a function as a string) for
#'     calculating an `average' value for each cluster; default is the
#'     \code{median} ## @param dev a function for calculating
#'     deviation from the the average, ## default is the standard
#'     deviation (‘sd’)
#' @param q either numeric (0.5-1, or 0), or a string. If numeric it
#'     indicates the fraction of data shown as transparent ranges,
#'     e.g. \code{q=0.9} indicates that 90% of data are within the
#'     ranges, i.e. the 5% and 95% quantiles are the lower and upper
#'     borders, and for q=0.5 the ranges correspond to the box limits
#'     in a boxplot. If q is a string it must be a function name for
#'     calculating variance (eg. "sd", "var"), which will be added and
#'     subtracted from the average (argument \code{avg}) for high and
#'     low data cut-offs.
#' @param rm.inf remove infinite values (e.g. from log
#'     transformations)
#' @export
clusterAverages <- function(ts, cls, cls.srt, avg="median", q=.9, rm.inf=TRUE) {

    if ( missing(cls.srt) )
      cls.srt <- sort(unique(cls))
    ## ensure present
    cls.srt <- cls.srt[cls.srt%in%cls]

    ## ensure use of rownames
    cls.srt <- as.character(cls.srt)

    ## set up result matrices
    clavg <- matrix(NA, ncol=ncol(ts), nrow=length(cls.srt))
    rownames(clavg) <- cls.srt
    cllow <- clhig <- clavg

    ## get requested upper/lower range function
    if ( is.numeric(q) ) {
        if ( (q<0.5 | q>1) & q!=0 )
            stop("only q within 0.5 and 1 (or q=0) is allowed; tip: if data is normally distributed you can use 'var', 'sd' or 'stderr' and combine with avg='mean'.")
        ## q is the fraction of shown values, (1-q)/2 are the calculated
        ## quantiles
        qf <- (1-q)/2
    } else {
        if ( q%in%c("sd","var","stderr") & avg!="mean" ) {
            avg <- "mean"
            warning("average avg='mean' enforced, since sd, var and stderr subtraction only make sense for mean.")
        }
        qf <- get(q, mode="function")
    }

    ## get requested average functions
    if ( mode(avg)!="function" )
        avg <- get(avg, mode="function")
    
    
    ## rm Inf, eg.
    ts[is.infinite(ts)] <- NA
    
    for ( cl in cls.srt ) {
        ## average and deviation
        clavg[cl,] <- apply(ts[cls==cl,,drop=FALSE],2,
                            function(x) avg(x,na.rm=T))
        if ( is.numeric(q) ) {
            ## upper/lower quantiles
            cllow[cl,]<- apply(ts[cls==cl,,drop=FALSE],2,
                               function(x) quantile(x,  qf,na.rm=T))
            clhig[cl,]<- apply(ts[cls==cl,,drop=FALSE],2,
                               function(x) quantile(x,1-qf,na.rm=T))
        } else  {
            df <- apply(ts[cls==cl,,drop=FALSE],2, function(x) qf(x,na.rm=T))
            cllow[cl,] <- clavg[cl,] - df
            clhig[cl,] <- clavg[cl,] + df
        }
        #cat(paste(cl,"\n"))
    }
    res <- list(avg=clavg, low=cllow, high=clhig,
                functions=list(average=avg, q=q))
    class(res) <- "clusteraverages"
    res
}

#' plot BIC from flowclusterTimeseries
#'
#' @param cls "clustering" object from function
#' @param norm logical, indicating whether the plotted BIC/ICL values
#' should be normalized by number of (clustered, \code{cluster!="0"})
#' data points
#' @param legend position of the legend in the plot
#' @param ... arguments to plot
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}
#' @export
plotBIC <- function(cls, norm=FALSE, legend="right", ...) {
    
    ## BIC/ICL
    bic <- cls$bic
    icl <- cls$icl
    K <- as.numeric(names(bic))
    krng <- range(K)
    max.clb <- K[which(cls$K==cls$max.clb)]
    max.cli <- K[which(cls$K==cls$max.cli)]
    sel.bic <- bic[which(cls$K==cls$selected)]
    sel <- K[which(cls$K==cls$selected)]
    
    if ( norm ) {
        
        N <- cls$N  # number of clustered data points
        M <- cls$M  # data dimension

        ## normalize BIC by N and M
        bic <- bic/N/M
        icl <- icl/N/M
        sel.bic <- sel.bic/N/M
    }
    max.bic <- max(bic, na.rm=TRUE) # max BIC
    max.icl <- max(icl, na.rm=TRUE) # max ICL

    ## legend
    leg <- c("BIC","ICL","failed","max. BIC","max. ICL", "selected")
    lpch <- c(1,1,4,4,4,19)
    llty <- c(1,NA,NA,NA,NA,NA)
    lcol <- c(1,4,2,1,4,2)
    lcex <- c(1,.5,1,1.5,1,1.5)
    
    ## flowMerge clusters present?
    if ( !is.null(cls$merged.K) ) {
        mrg.cl <- cls$merged.K
        mrg.origK <- sub("K:","",names(cls$merged.origK))
        mrg.orig <- cls$merged.origK
        mrg.bic <- bic[which(cls$K==mrg.origK)]
        krng <- range(c(krng,mrg.orig,mrg.cl))
    }
    
    plot(K, bic, ylim=range(c(bic,icl),na.rm=T),xlab="K",xlim=krng,
         ylab=paste0(ifelse(norm,"norm. ",""),"BIC/ICL"), axes=FALSE, ...)
    axis(2)
    axis(1,at=1:max(krng))
    lines(K,bic)
    points(K[is.na(bic)], rep(min(bic,na.rm=T),sum(is.na(bic))), pch=4, col=2)
    points(max.clb, max.bic,lty=2,pch=4,cex=1.5)
    abline(v=max.clb,lty=2)
    if ( any(!is.na(icl)) ) {
        points(K, icl, col=4,cex=.5)
        points(max.cli, max.icl,lty=2,pch=4,cex=1,col=4)
    }
    points(sel, sel.bic, col=2, pch=19, cex=1.5)
    legend(legend,legend=leg,pch=lpch,lty=llty,col=lcol,pt.cex=lcex,
           bg="#FFFFFF77")
    if ( !is.null(cls$merged.K) ) {
        arrows(x0=mrg.orig,y0=mrg.bic, x1=mrg.cl,y1=mrg.bic)
        abline(v=c(mrg.orig,mrg.cl),lty=2)
        text(mrg.cl, max.bic,"merge", pos=2,cex=.9)
    }
}

#' plot indivividual time series in cluster context
#'
#' plot indivividual time series (GOI: genes of interest) from
#' timeseries and clustering objects; a wrapper for \code{\link{plotClusters}},
#' allowing to only plot individual time series in their cluster context
#' (colors and individual panels for clusters) and using the same style and
#' providing all functionality from \code{\link{plotClusters}}.
#' @param x either a simple data matrix with time points (x-axis) in columns,
#' or a processed time-series as provided by
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}}
#' @param cls "clustering" object from function
#' \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}} or
#' \code{\link[segmenTier:flowclusterTimeseries]{flowclusterTimeseries}}
#' @param goi list of feature ids (rownames in cls$clusters !) to plot
#' @param grep logical, if TRUE \code{goi} are searched by pattern
#' match (as argument \code{pattern} in
#' \code{\link[base:grep]{grep}}) instead of perfect matches
#' @param each plot separate panels for each cluster
#' @param plot.legend logical indicating whether to plot feature names in a
#' legend
#' @param lwd line width of single time-series
#' @param leg.xy position of the legend, see
#' \code{\link[graphics:legend]{legend}}
#' @param y.intersp legend line interspacing, see
#' \code{\link[graphics:legend]{legend}}
#' @param ... arguments to \code{\link{plotClusters}}
#' @export
plotSingles <- function(x, cls, goi, grep=FALSE,
                        each=TRUE, plot.legend=each, lwd=2,
                        leg.xy="topleft", y.intersp=1, ...) {
    ## TODO: set all genes
    ## in cls$cluster to "-1"

    ## rm none-present goi
    rm <- !goi%in%rownames(cls$clusters)
    if ( grep )
      rm <- sapply(goi, function(x) length(grep(x, rownames(cls$clusters)))==0)
    if ( sum(rm) )
      cat(paste(paste(goi[rm],collapse=";"),"not found\n"))
    goi <- goi[!rm]
    
    if ( length(goi)==0 ) {
        cat(paste("no GOI found\n"))
        return(NULL)
    }

    ## get goi and set all other clusterings to -1
    kp <- which(rownames(cls$clusters)%in%goi)
    if ( grep )
      kp <- unlist(sapply(goi, grep, rownames(cls$clusters)))
    cls$clusters[-kp,] <- -1

    ## remap goi names to use for legend ids
    if ( is.null(names(goi)) )
        names(goi) <- goi
    leg.ids <- names(goi)
    names(leg.ids) <- goi
    
    avg <- plotClusters(x, cls, avg.col=NA, lwd=lwd, avg.lwd=0, each=each,
                        alpha=1, use.lty=TRUE, type=c("all"),
                        plot.legend=each&plot.legend,
                        leg.xy=leg.xy, leg.ids=leg.ids, ...)
    leg <- do.call(rbind,avg$legend)
    if ( !is.null(names(goi)) ) {
        if ( grep ) {
            nms <- names(goi)[unlist(sapply(goi, grep, leg[,"id"]))]
            leg[!is.na(nms),"id"]  <- nms[!is.na(nms)]
        } else
            leg[,"id"] <- names(goi)[match(leg[,"id"], goi)]
    }
    ## TODO: auto-select y-intersp if too many goi
    if ( !each & plot.legend )
      legend(leg.xy, legend=leg[,"id"],
             lty=leg[,"lty"], col=leg[,"col"], lwd=lwd,
             bg="#FFFFFFAA",bty="o", y.intersp=y.intersp)
    avg$legend <- leg
    invisible(avg)
}

#' plots cluster averages 
#'
#' plots average time series of clusters as calculated by
#' \code{\link{clusterAverages}}, including variations around the
#' average as transparent areas, or all individual data.  The function
#' is quite flexible and allows to normalize the data and set
#' automatic y-limit selections.
#' @param x either a simple data matrix with time points (x-axis) in
#'     columns, or a processed time-series as provided by
#'     \code{\link[segmenTier:processTimeseries]{processTimeseries}}
#' @param cls eiter a vector (\code{length(cls)==nrow(x)}) or a
#'     structure of class 'clustering' as returned by segmenTier's
#'     \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param k integer or string specifiying the clustering (k: cluster
#'     numbers) to be used if cls is of class 'clustering'; if missing
#'     (default) the `selected' clustering from \code{cls} is chosen
#'     (see \code{\link{selected}}).
#' @param type string specifying the type of plot: "rng" for plotting
#'     only data ranges (see argument \code{q}) or "all" to plot each
#'     individual time-course (as thin lines)
#' @param fill fill data ranges with transparent version of cluster
#'     color
#' @param border line width of border line at data ranges
#' @param each logical value indicating whether to plot all cluster
#'     averages on one panel (\code{FALSE}) or each cluster on a
#'     separate panel (\code{TRUE})
#' @param time optional numeric vector specifiying x-axis time
#'     coordinates
#' @param time.at argument \code{at} for the x-axis (\code{axis(1,
#'     at=time.at)})
#' @param avg a function (or the name of a function as a string) for
#'     calculating an `average' value for each cluster; default is the
#'     \code{median}
#' @param norm normalization of the time-series data, must be a
#'     function that transforms the data, available via
#'     \link{segmenTools} are \code{lg2r}, \code{ash}, \code{log_1},
#'     \code{meanzero} normalizations
#' @param q either numeric (0.5-1, or 0), or a string. If numeric it
#'     indicates the fraction of data shown as transparent ranges,
#'     e.g. \code{q=0.9} indicates that 90% of data are within the
#'     ranges, i.e. the 5% and 95% quantiles are the lower and upper
#'     borders, and for q=0.5 the ranges correspond to the box limits
#'     in a boxplot. If q is a string it must be a function name for
#'     calculating variance (eg. "sd", "var"), which will be added and
#'     subtracted from the average (argument \code{avg}) for high and
#'     low data cut-offs.  Note that this parameter can also influence
#'     \code{ylim}, the limits of the y axis.
#' @param cls.srt optional cluster sorting, can be used for selection
#'     of subsets of clusters; if cls is of class 'clustering' it is
#'     taken from there
#' @param cls.col optional named cluster color vector, where names
#'     indicate the clusters (as.caracter(cls)); if cls is of class
#'     'clustering' it is taken from there
#' @param axes add axes (bottom and left)
#' @param xlab x-axis label (auto-selected if missing)
#' @param xlim \code{xlim} parameter for plot
#' @param ylab y-axis label (only used if \code{each==FALSE})
#' @param ylab.size add cluster size to cluster-wise ylab
#' @param ylab.cex font size of ylab
#' @param ylim either conventional range of the y-axis, or a string
#'     specifying whether ylim should be calculated from the average
#'     (\code{ylim="avg"}), for all data (\code{ylim="all"}), or from
#'     the lower/upper ranges (\code{ylim="rng"}); in the latter case
#'     the y-axis depends on argument \code{q}
#' @param ylim.scale if \code{ylim=="avg"}, the calculated ylim will
#'     be extended by this fraction of the total range on both sides
#' @param avg.col color for average line; used only if
#'     \code{type="all"}
#' @param avg.lty line type for average plots
#' @param avg.lwd line width for average plots
#' @param avg.pch point symbol for average plots
#' @param avg.cex point size for average plots
#' @param lwd line width for indidiual time series plots (if
#'     \code{type=="all"})
#' @param use.lty use individual line types (if \code{type=="all"});
#'     this is only useful for very small clusters and is mainly used
#'     in the \code{\link{plotSingles}} wrapper
#' @param alpha set alpha value of range or individual time series
#'     colors (color opaqueness)
#' @param plot.legend add a legend, useful for very small clusters and
#'     mainly used in the \code{\link{plotSingles}} wrapper
#' @param embed logical, if TRUE and argument \code{each=TRUE} (one
#'     plot for each cluster), the automatic is suppressed allowing to
#'     embed multiple plots (for each cluster) into external
#'     \code{\link[graphics:layout]{layout}} or
#'     \code{par(mfcol/mfrow)} setups
#' @param leg.xy position of the legend, see
#'     \code{\link[graphics:legend]{legend}}
#' @param sample.ticks draw ticks on top axis for each actual sample
#' @param vline x coordinates for vertical lines, useful e.g. to mark
#'     replicates
#' @param vl_col color for vertical line (default to cluster colour)
#' @param vl_lty vertical line type
#' @param vl_lwd vertical line width
#' @param ref.xy reference data (either x or xy coordinates, see
#'     \code{\link[grDevices:xy.coords]{xy.coords}}) for reference
#'     data to be plotted in the background, e.g. light/dark phases or
#'     dissolved oxygen traces
#' @param ref.col color for reference data plot
#' @param ref.ylab y-axis label (right y-axis) for reference data
#' @param ref.log logarithmic axis for reference data
#' @param leg.ids a named vector providing alternative IDs for
#'     legends; the names should correspond to the rownames of
#'     clusterings in \code{cls}
#' @param ... further arguments to the basic setup call to
#'     \code{\link[graphics:plot]{plot}} ## TODO: clean up mess
#'     between plot.clustering, plot.clusteraverages and this ##
#'     plot.clusteraverages should become a function of
#'     plot.clustering, ## and clusterAverages can be a private
#'     function!  ## make function `timeseriesPlot' or `clusterPlot',
#'     that takes ## either tset/cset or matrix/vector
#' @export
plotClusters <- function(x, cls, k, each=TRUE, type="rng",
                         fill=TRUE, border=0,
                         time, time.at,
                         avg="median",  q=.9, norm, 
                         cls.col, cls.srt,  
                         axes=TRUE, xlab, xlim,
                         ylab, ylab.size=TRUE, ylab.cex=1,
                         ylim=ifelse(each,"avg","rng"), ylim.scale=.1,
                         avg.col="#000000",avg.lty=1,avg.lwd=3,
                         avg.cex=1,avg.pch=1,
                         lwd=.5, use.lty=FALSE, alpha=.2,
                         embed=FALSE,
                         plot.legend=FALSE, leg.xy="topleft", leg.ids,
                         sample.ticks=TRUE,
			 vline='',vl_col = 0,vl_lwd=3,vl_lty = 1,
                         ref.xy, ref.col="#C0C0C080", ref.ylab="",
                         ref.log=FALSE, ...) 
{

    
    if ( inherits(cls, "clustering") ) {

        if ( missing(k) )
            k <- selected(cls, name=TRUE)
        if ( is.numeric(k) )
            k <- selected(cls, K=k, name=TRUE)

        ## TODO: why import problem in R CMD check?
        ##if ( !"colors" %in% names(cls) )
        ##    cls <- segmenTier::colorClusters(cls)
        ##if ( !"sorting" %in% names(cls) )
        ##    cls <- segmenTier::sortClusters(cls)
        cls.col <- cls$colors[[k]]
        if( missing(cls.srt) )
            cls.srt <- cls$sorting[[k]]

        cls <- cls$clusters[,k]
    } else {

        ## cluster sorting
        if ( missing(cls.srt) )
            cls.srt <- sort(unique(cls))
        ## cluster colors
        if ( missing(cls.col) ) {
            cls.col <- rainbow(length(cls.srt))
            names(cls.col) <- cls.srt
        }
    }
    ## convert to character for use as rownames
    cls.srt <- as.character(cls.srt)
    cls <- as.character(cls)

    ## filter for present clusters (allows plotSingles)
    cls.srt <- cls.srt[cls.srt%in%unique(cls)]
    cls.sze <- table(cls)
    
    ## average and polygon colors
    pol.col <- add_alphas(cls.col,rep(alpha,length(cls.col)))
    if ( type[1]=="rng" ) 
        avg.col <- add_alphas(cls.col,rep(1,length(cls.col)))
    if ( type[1]=="all" ) {
        avg.col <- rep(avg.col, length(cls.col))
        names(avg.col) <- names(cls.col)
    }

    ## TIME SERIES
    if ( inherits(x, "timeseries") )
        ts <- x$ts
    else ts <- data.matrix(x)

    ## add rownames if missing
    if ( any(is.null(rownames(ts))) )
      rownames(ts) <- 1:nrow(ts)
    
    ## normalize
    if ( !missing(norm) )
        ts <- get(norm,mode="function")(ts)
    else norm <- "raw"

    ## y-axis
    ##if ( missing(ylab) )
    ##  ylab <- norm
    
    ## x-axis
    if ( missing(time) ) {
        time <- 1:ncol(ts)
        if ( missing(xlab) )
            xlab <- "index"
    } 
    if ( missing(xlim) )
        xlim <- range(time)
    if ( missing(xlab) )
        xlab <- "time"
    if ( missing(time.at) )
        time.at <- pretty(time)

    ## set all but goi to 
    #if ( !missing(goi) ) {
    #    #cls.col[rownames(ts)%in%goi
    #}

    #if ( 
    all.col <- pol.col[cls]
    #all.lty <- rep(1, length(cls))
    
    ## average values
    if ( type[1]=="all" ) q <- 1
    avg <- clusterAverages(ts, cls=cls, cls.srt=cls.srt, avg=avg, q=q)
    #cat(paste(cls.srt))

    ## calculate ylim from full ranges? 
    if ( typeof(ylim)=="character" ) {
        if ( ylim == "rng" )
            ylim <- c(min(avg$low[cls.srt,],na.rm=TRUE),
                      max(avg$high[cls.srt,],na.rm=TRUE))
        else if ( ylim == "avg" ) {
            ylim <- range(avg$avg[cls.srt,],na.rm=TRUE)
            ylim <- c(ylim[1]-diff(ylim)*ylim.scale,
                      ylim[2]+diff(ylim)*ylim.scale)
        }
        else if ( ylim=="all" ) {
            ylim <- range(ts, na.rm=TRUE)
        }
    }

    ## PLOT
    ## plot each cluster on separate panel?
    if ( each ) {
        mai <- par("mai")
        mfc <- par("mfcol")
        newmai <- c(0,mai[2],0,mai[4])
        if ( !embed ) # don't set mfcol, to embed into externally set layouts
            par(mfcol=c(length(cls.srt),1),mai=newmai)
     } else {
         if ( !missing(ref.xy) ) {
             ylog <- ifelse(ref.log,"y", "")
             rlm <- range(ref.xy[,2])
             if ( ylog!="y" ) rlm <- round(rlm)
             plot(ref.xy,type="l",lty=1,lwd=2,col=ref.col,
                  axes=FALSE,xlab=NA,ylab=NA,
                  ylim=rlm,xlim=xlim, log=ylog)
             polygon(x=c(ref.xy[1,1],ref.xy[,1],ref.xy[nrow(ref.xy),1]),
                    y=c(min(ref.xy[,2]),ref.xy[,2],min(ref.xy[,2])),
                    col=ref.col,border=NA)
             if ( axes ) {
                 if ( ylog=="y" )
                     axis(4,labels=FALSE, las=2,
                          col=ref.col, col.ticks=ref.col, col.axis=ref.col)
                 axis(4,at=rlm,las=2,
                      col=ref.col, col.ticks=ref.col, col.axis=ref.col)
                 mtext(ref.ylab, 4, .35, col=ref.col)
             }
             par(new=TRUE)
         }
        ## use normalization as ylab
        if ( missing(ylab) ) ylab <- norm
        plot(1,col=NA,axes=FALSE,
             xlab=xlab,xlim=xlim,
             ylab=NA,ylim=ylim, ...)
         mtext(ylab, 2, par("mgp")[1], cex=ylab.cex)
         if ( axes ) {
             axis(1, at=time.at);axis(2)
             if ( sample.ticks ) {
                 axis(3, at=time, labels=FALSE, tcl=-par("tcl"))
                 mtext("samples", 3, 1.2)
             }
        }
    }
    ## plot each cluster in cls.srt
    used.pars <- rep(list(NA), length(cls.srt))
    names(used.pars) <- cls.srt
    for ( cl in cls.srt ) {
        if ( each ) {
            if ( !missing(ref.xy) ) {
                ylog <- ifelse(ref.log,"y", "")
                rlm <- range(ref.xy[,2])
                if ( ylog!="y" ) rlm <- round(rlm)
                plot(ref.xy,type="l",lty=1,lwd=2,col=ref.col,
                     axes=FALSE,xlab=NA,ylab=NA,
                     ylim=rlm,xlim=xlim, log=ylog)
                polygon(x=c(ref.xy[1,1],ref.xy[,1],ref.xy[nrow(ref.xy),1]),
                        y=c(min(ref.xy[,2]),ref.xy[,2],min(ref.xy[,2])),
                        col=ref.col,border=NA)
                if ( axes ) {
                    if ( ylog=="y" ) 
                        axis(4,labels=FALSE, las=2,
                             col=ref.col, col.ticks=ref.col, col.axis=ref.col)
                    axis(4,at=rlm,las=2,
                         col=ref.col, col.ticks=ref.col)
                    mtext(ref.ylab, 4, .35, col=ref.col, col.axis=ref.col,
                          cex=ylab.cex)
                }
                par(new=TRUE)
            }
            if ( missing(ylab) )
                if ( ylab.size )
                    ylb <- paste(cl," (",cls.sze[cl],")",sep="")
                else
                    ylb <- cl
            else ylb <- ylab
            plot(1,col=NA,axes=FALSE, xlab=NA, xlim=xlim,
                 ylab=NA,ylim=ylim, ...)
            mtext(ylb, 2, par("mgp")[1], cex=ylab.cex)
            if ( axes ) {
                axis(1, at=time.at);axis(2)
            }
        }
        if ( "rng"%in%type ) {## polygon
            pcol <- ifelse(fill, pol.col[cl], NA)
            x <- c(time,rev(time))
            y <- c(avg$low[cl,],rev(avg$high[cl,]))
            if ( any(is.na(y)) ) warning("some data in plot ranges were NA")
            polygon(x=x[!is.na(y)], y=y[!is.na(y)],
                    col=pcol,border=NA)
            if ( border>0 ) {
                lines(time, avg$low[cl,], col=cls.col[cl], lwd=border)
                lines(time, avg$high[cl,], col=cls.col[cl], lwd=border)
            }
        }
        if ( "all"%in%type ) {
            idx <- cls==cl
            if ( use.lty )
                lty <- rep(1:6, length.out=sum(idx,na.rm=TRUE))
            else  lty <- rep(1, sum(idx,na.rm=TRUE)) #lty <- all.lty[idx]
            matplot(time, t(ts[idx,,drop=FALSE]), add=TRUE,
                    type="l", lty=lty, col=all.col[idx], lwd=lwd)
            
            ## store for external legend (eg. plotSingles with each=FALSE)
            used.pars[[cl]] <- data.frame(id=rownames(ts)[idx],
                                          lty=lty,col=all.col[idx],
                                          stringsAsFactors=FALSE)
            if ( plot.legend ) {
                nms <- rownames(ts)[idx]
                if ( !missing(leg.ids) ) {
                    nms[nms%in%names(leg.ids)] <-
                        leg.ids[nms[nms%in%names(leg.ids)]]
                }
                legend(leg.xy, nms, bg="#FFFFFFAA",bty="o",
                       col=all.col[idx], lty=lty, lwd=lwd)
            }
        }
        ## average last
        lines(time, avg$avg[cl,], lwd=avg.lwd, lty=avg.lty, col=avg.col[cl])
        points(time, avg$avg[cl,], pch=avg.pch, cex=avg.cex, col=avg.col[cl])
        ## plot decoration
        if ( each ) {
            if(vl_col==0) abline(v=vline,col=cls.col[cl],lwd=vl_lwd,lty=vl_lty)
            else abline(v=vline,col=vl_col,lwd=vl_lwd,lty=vl_lty)
        }
    }
    ## plot decoration
    ## TODO: allow reference data (DO traces, LD phases)
    if ( !each ) {
 	if(vl_col==0) vl_col <- 1
	abline(v=vline,col=vl_col,lwd=vl_lwd,lty=vl_lty)
    }
    
    ## reset plot pars
    if ( each & !embed) {
        ## this reset prohibits plotting of legends etc.
        par(mai=mai,mfcol=mfc)
    }
        
    avg$normalization <- norm
    avg$ylim <- ylim
    avg$legend <- used.pars
    invisible(avg) # silent return of averages
    #avg.col
}

#' plots cluster averages 
#' 
#' plots average time series of clusters as calculated by
#' \code{\link{clusterAverages}}, including the variations around the mean
#' as polygons
#' @param x cluster time series average object as calculated by
#' \code{\link{clusterAverages}}
#' @param cls.srt optional sorting of clusters; clusters will be plotted
#' in this order, i.e. the last in \code{cls.srt} is plotted last;
#' can be used for selecting a subset of clusters
#' @param cls.col optional coloring of clusters
#' @param each logical value indicating whether to plot all cluster
#' averages on one panel (\code{FALSE}) or each cluster on a separate panel
#' (\code{TRUE})
#' @param ranges logical indicating whether to plot the ranges `high' and
#' `low' in the cluster average object \code{avg}
#' @param border line width of border line at data ranges 
#' @param xlab x-axis label
#' @param time optional time-points of the x-axis
#' @param ylab y-axis label
#' @param ylim either conventional range of the y-axis, or a string
#' specifying whether ylim should be calculated from the average
#' (\code{ylim="avg"}) or from the lower/upper ranges (\code{ylim="rng"})
#' @param ylim.scale if \code{ylim=="avg"}, the calculated ylim will be
#' extended by this fraction of the total range on both sides
#' @param ... arguments to plot
#' @export
plot.clusteraverages <- function(x, cls.srt, cls.col,
                                each=FALSE, ranges=TRUE, border=0,
                                xlab, time, 
                                ylab="average", ylim=ifelse(each,"avg","rng"),
                                ylim.scale=.1,...) {

    ##TODO: this should replace plotClusters?
    
    avg <- x
    ## x-axis
    if ( missing(time) ) {
        time <- 1:ncol(avg$avg)
        if ( missing(xlab) )
          xlab <- "index"
    } 
    if ( missing(xlab) )
      xlab <- "time"

    ## cluster sorting
    if ( missing(cls.srt) )
      cls.srt <- rownames(avg$avg)

    ## cluster colors
    if ( missing(cls.col) ) {
        cls.col <- rainbow(length(cls.srt))
        names(cls.col) <- cls.srt
    }

    ## average and polygon colors
    pol.col <- add_alphas(cls.col,rep(.2,length(cls.col)))
    avg.col <- add_alphas(cls.col,rep(1,length(cls.col)))

    ## calculate ylim from full ranges? 
    if ( typeof(ylim)=="character" ) {
        if ( ylim == "rng" )
            ylim <- c(min(avg$low[cls.srt,],na.rm=TRUE),
                      max(avg$high[cls.srt,],na.rm=TRUE))
        else if ( ylim == "avg" ) {
            ylim <- range(avg$avg[cls.srt,])
            ylim <- c(ylim[1]-diff(ylim)*ylim.scale,
                      ylim[2]+diff(ylim)*ylim.scale)
        }
    }
    
    ## plot each cluster on separate panel?
    if ( each ) {
        mai <- par("mai")
        mfc <- par("mfcol")
        mai[c(1,3)] <- 0
        par(mfcol=c(length(cls.srt),1),mai=mai)
    } else {
        plot(1,col=NA,axes=FALSE,
             xlab=xlab,xlim=range(time),
             ylab=ylab,ylim=ylim, ...)
        axis(1);axis(2)
        axis(3, at=time, labels=FALSE, tcl=-par("tcl"))
        mtext("samples", 3, 1.2)
    }
    ## plot each cluster in cls.srt
    for ( cl in cls.srt ) {
        if ( each ) {
            plot(1,col=NA,axes=FALSE, xlab=NA,xlim=range(time),
                 ylab=paste(ylab,cl,sep=" - "),ylim=ylim, ...)
            axis(1);axis(2)
        }
        if ( ranges ) {
            x <- c(time,rev(time))
            y <- c(avg$low[cl,],rev(avg$high[cl,]))
            if ( any(is.na(y)) ) warning("some data in plot ranges were NA")

            polygon(x=x[!is.na(y)],y=y[!is.na(y)], col=pol.col[cl],border=NA)
            if ( border>0 ) {
                lines(time, avg$low[cl,], col=cls.col[cl], lwd=border)
                lines(time, avg$high[cl,], col=cls.col[cl], lwd=border)
            }
        }
        lines(time, avg$avg[cl,], col=avg.col[cl], lwd=3) ## average last
    }
    ## reset plot pars
    if ( each ) {
        par(mai=mai,mfcol=mfc)
    }
}




#' plot sorted clustering as color table
#' @param cset a structure of class 'clustering' as returned by
#' segmenTier's \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param k integer or string specifiying to clustering (K: cluster numbers)
#' to be used if cls is of class 'clustering'; if missing (default) the
#' `selected' clustering from \code{cls} is chosen
#' @param dir direction, "h" for horizontal (default), "v" for vertical
#' @param ... arguments to \code{\link{image_matrix}}
#' @export
plotClusterLegend <- function(cset, k=selected(cset), dir="h", ...) {

    ## sort colors by numeric cluster index
    cols <- lapply(cset$colors,
                   function(x) x[as.character(sort(as.numeric(names(x))))])
    ## get clusters as matrix and subtract .5
    mat <- matrix(cset$clusters[,k,drop=FALSE]-.5, nrow=1)
    ## calculate breaks
    breaks <- as.numeric(names(cols[[k]]))
    breaks <- c(min(breaks)-1,breaks)
    ## TODO: allow multiple k
    ## (for loop with add=TRUE, and NA in all non-k columns)
    if ( dir=="v" ) mat <- t(mat)
    image_matrix(mat, col=cols[[k]], breaks=breaks,...)
    ## TODO: use image_cls below
    ## image_cls(mat, col=cols[[k]], dir=dir,...)
}

image_cls <- function(cls, col, dir="h", ...) {

    mat <- matrix(cls-.5, nrow=1)
    breaks <- sort(as.numeric(unique(cls)))
    breaks <- c(min(breaks)-1,breaks)
    ## TODO: allow multiple k
    ## (for loop with add=TRUE, and NA in all non-k columns)
    if ( dir=="v") mat <- t(mat)
    image_matrix(mat, col=col, breaks=breaks,...)
}

### UTILS

## from lib TeachingDemos; used in plotFeatures
#' plot borders around text based on text colors
#' @param x x-coordinates of text labels
#' @param y y-coordinates of text labels
#' @param labels vector of text labels
#' @param col vector of text foreground colors
#' @param bg vector of text background colors
#' @param auto.bg auto-select background color black or white based
#' on a brightness color model
#' @param bg.thresh brightness threshold to switch from black to white
#' background color
#' @param r ratio of stringwidth and stringheight to use for background color
#' @param ... further arguments to \code{\link{text}} in both foreground
#' and background calls
#'@export
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       auto.bg=TRUE, bg.thresh=.75, r=0.05, ... ) {

    if ( auto.bg ) {
        bgn <- rep(bg, length(col))
        for ( i in 1:length(col) ) {
            crgb <- col2rgb(col[i])/255
            L <- 0.2126 * crgb[1,1] + 0.7152 * crgb[2,1] + 0.0722 * crgb[3,1]
            bgn[i] <- ifelse(L>bg.thresh, "#000000", "#FFFFFF")
        }
        bg <- bgn
    }
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  

    theta <- seq(pi/4, 2*pi, length.out=100)
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}



#' replace alpha values of an RGB string color vector
#'
#' adds (or replaces existing) alpha values (from 0 to 1) to
#' RGB string colors; extended to vectors from \url{http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html}
#' @param col a string vector of RGB color specifiers with ("#FF0000AA ")
#' or without ("#FF0000") existing alpha values
#' @param alpha a numeric vector of equal length as argument \code{col}
#' providing new alpha values from 0 to 1
#' @export
add_alphas <- function(col, alpha=rep(1,length(col))){
    nms <- names(col)
    nacol <- which(is.na(col))
    lcol <- lapply(col, function(x) c(col2rgb(x)/255))
    col <- sapply(1:length(col), function(x) {
        y <- lcol[[x]]
        rgb(y[1], y[2], y[3], alpha=alpha[x])})
    names(col) <- nms
    col[nacol] <- NA
    col
}

