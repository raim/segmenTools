
### TOOLS for handling feature-based clusterings

## TODO:
## calculate average phases via phaseDist or directly
## auto-select oscillatory clusters from p-values
## clusterCluster/image_matrix : make higher level interface, see
## plot.hypergeoTables in analyzeSegments_2017.R
## integrate GO wrapper for clusterCluster

### COMPARE CLUSTERINGS

#' calculates overlaps between two clusterings
#' calculates mutual overlaps between two clusterings of the same data set
#' using hypergeometric distribution statistics for significantly
#' enriched or deprived mutual overlaps
#' @param cl1 clustering 1
#' @param cl2 clustering 2
#' @param na.string replace NA or empty strings by `<na.string>'
#' @param cl1.srt optional cluster sorting of clustering 1
#' @param cl2.srt optional cluster sorting of clustering 2
#' @param req.vals requested statistics; one of `greater', `less'
#' or `two-sided' to calculate enrichment, depletion or the more
#' more signficant of both, respectively
#'@export
clusterCluster <- function(cl1, cl2, na.string="na", cl1.srt, cl2.srt,
                           req.vals=c("greater")) {
  ## check cluster length
  if ( length(cl1) != length(cl2) ) {
      print(paste("ERROR cluster vectors of different size:",
                  length(cl1),length(cl2)))
      return(NULL)
  }
  
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
  do.prich <- sum(req.vals=="greater")>0
  do.ppoor <- sum(req.vals=="less")>0
  if ( sum(req.vals=="two.sided")>0 )
    do.ppoor <- do.prich <- TRUE
  
  ## get clusters
  f1 <- levels(as.factor(cl1))
  f2 <- levels(as.factor(cl2))
  if ( !missing(cl1.srt) ) f1 <-  cl1.srt
  if ( !missing(cl2.srt) ) f2 <-  cl2.srt

      

  ## result matrices
  overlap <- matrix(NA, nrow=length(f1), ncol=length(f2))
  rownames(overlap) <- f1
  colnames(overlap) <- f2
  percent <- overlap
  if ( do.prich ) prich <- overlap
  if ( do.ppoor ) ppoor <- overlap
  
  
  for ( i in 1:length(f1) ) { # white and black balls
      
      m=sum(cl1==f1[i]); # number of white balls
      n=sum(cl1!=f1[i]); # number of black balls
      
      for ( j in 1:length(f2) ) { # balls drawn
          
          q <- sum(cl1==f1[i] & cl2==f2[j]); # white balls drawn
          k <- sum(cl2==f2[j]); # number of balls drawn
          
          overlap[i,j] <- q
          percent[i,j] <- round(100*q/k, digits = 2)
          
          ## Calculate cumulative HYPERGEOMETRIC distributions 
          
          ## enrichment:  p-value of finding q or more white balls,
          ## P[X >= x]
          if ( do.prich ) {
              prich[i,j] <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=F)
              if ( is.na(prich[i,j]) ) 
                prich[i,j]=1
          }
          
          
          ## deprivation: p-value for finding q or less white balls,
          ## P[X <= x]
          if ( do.ppoor ) {
              ppoor[i,j] <- phyper(q=q, m=m, n=n, k=k, lower.tail=T);
              if ( is.na(ppoor[i,j]) )
                ppoor[i,j]=1;
          }
      }
  }
  
  result <- list(overlap=overlap,percent=percent)

  ## append p-values
  p.value <- prich
  if ( do.prich & !do.ppoor ) p.value <- prich
  if ( do.ppoor & !do.prich ) p.value <- ppoor

  ## take smaller p-value, if both are requested!
  if ( do.ppoor &  do.prich )
    for ( prow in 1:nrow(prich) )
      for ( pcol in 1:ncol(prich) )
        p.value[prow,pcol] <- min(prich[prow,pcol],ppoor[prow,pcol])
  
  result <- append(result, list(p.value=p.value))
  
  
  class(result) <- "clusterOverlaps"
  return(result)
}

## Does hypergeometric distribution tests between a given clustering
## and a table of diverse other clusterings (character/factor table),
## e.g. single GO annotations, or other classifications.
## TODO : fuse with SCI and externalStats.
## TODO : fuse with contStatTable and switch via req.vals
## ARGUMENTS
## cls: a character or numeric vector with the main clustering
## data: a matrix with other clusterings (nrow(data)==length(cluster)
## to which the main clustering; or e.g. a TRUE/FALSE table
## TODO: generalize the latter case
## should be compared, hypergeometric distribution will be tested
## for all combinations of clusters
## summary.table: construct summary tables of 
## summary.filter: a regex pattern for classifications NOT to be included
## in the results summary, e.g. "FALSE" in simple logical tables
## add.expected: for simple binary tables with values "TRUE"/"FALSE" or
## colnames(data)/"not_".colnames(data), add both values to
## results, even if a given clustering in data contains only one of those!
## ...: arguments to clusterCluster, e.g. plot.type
## TODO 2017: move away from T/F table, use simple ;-sep string of
## annotation terms instead; make nicer report as for lehmann13_geneClusters.R
## and provide pval table for image_matrix
clusterAnnotation <- function(cls, data, cls.srt,
                              summary.table=FALSE, summary.filter="^not_",
                              summary.val="enrichment",
                              add.expected=FALSE,
                              p.thresh=0.01, 
                              verbose=TRUE,
                              ...) {

    ## sorted list of clusters!
    if ( missing(cls.srt) ) {
        cls.srt <- unique(cls)

        ## sort numeric if clusters can be converted!
        ## TODO: find better way to find whether clusters can be cast to numeric
        ## TODO: repair for cases where clusters contain NA
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
    overlap <- NULL ## CONTIGENCY TABLE!
    pvalues <- NULL ## hypergeometric distribution test p.value
    hyp.names <- NULL
    hyp.img <- NULL # store images to link from tables
    ## TODO : parallelize this
    if (verbose) cat(paste("HYPERGEOMETRIC STATISTICS FOR:\n",
                           ncol(data), "categories\t:\n"))
    for ( j in 1:ncol(data) ) {
        bins <- data[,j]
        name <- colnames(data)[j]
        expected <- unique(bins)

        if ( add.expected ) {
            ## TRUE/FALSE
            if ( sum(expected %in% c("TRUE","FALSE" ))>0 )
              expected <- unique(c(expected, "TRUE","FALSE"))
            ## name and "not_".name
            else
              expected <- unique(c(expected, name, paste("not",name,sep="_")))
        }

        if (verbose) cat(paste(j, "-",name, ", ", sep=""))
        
        ## cumulative hypergeometric distribution
        tmp <- clusterCluster(bins, cls, req.vals=c("greater"), ...)
       
        # get overlap and p.values in original bin order
        if ( !is.null(tmp) ) {
            ovl <- tmp$overlap # contingency table
            pvl <- tmp$p.value
            ## add expected values!
            if ( add.expected )
              for ( exp in expected )
                if ( ! exp %in% rownames(ovl) ) {
                    nam <- rownames(ovl)
                    ovl <- rbind(ovl, rep(0, ncol(ovl)))
                    pvl <- rbind(pvl, rep(1, ncol(ovl)))
                    rownames(ovl) <- c(nam,exp)
                    rownames(pvl) <- c(nam,exp)
                }
              
            overlap <- rbind(overlap, ovl)
            pvalues <- rbind(pvalues, pvl)
            # copy name in same size to allow easy cbind below
            hyp.names<-c(hyp.names, rep(name, nrow(ovl)))
          } else {
            if (verbose)
              cat(paste("WARNING: hypergeo test failed for bin:", name, "\n"))
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
    # fraction of significant (by option p.thresh) p-values

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
    psig["total",] <- c(total.clustered, sum(hp<p.thresh)/length(hp))

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
            #filter <- pvalues[, j] < p.thresh
            psig[j, "h"] <- sum(pvalues[, j] < p.thresh)/nrow(pvalues)

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
        
            colnames(hyp.table) <- c("rank", "data name", "bin", 
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

    ## collect results into one table
    ## for amounts of each cluster/category overlap
    ## and one table for a p-value of this overlap
    if ( summary.table ) {
        for ( k in 1:length(hyp.tables) ) {            
            tab <- hyp.tables[[k]]

            if ( k==1 ) {
                bins <- tab[,"bin"]
                if ( summary.filter != "" )
                  bins <- grep(summary.filter, bins, invert=TRUE, value=TRUE)
                
                bin.cls <- matrix(NA, nrow=length(hyp.tables),
                                   ncol=length(bins))
                rownames(bin.cls) <- names(hyp.tables)
                colnames(bin.cls) <- bins
                # convert to logic filter
                bins <- tab[,"bin"] %in% bins
                bin.pvl <- bin.cls
                bin.tot <- bin.cls
                
            }
            
            bin.cls[k,] <- as.numeric(tab[bins,summary.val]) # eg. "enrichment"
            bin.tot[k,] <- as.numeric(tab[bins,"bin in cluster"]) 
            bin.pvl[k,] <- as.numeric(tab[bins,"p-value"])
        }
        summary.table <- list(type=summary.val,
                              statistic=bin.cls, p.value=bin.pvl, total=bin.tot)
    }

    if (verbose) cat(paste(" ... done;\n"))
    
    
    # cluster averages
    only.clusters <- rownames(psig)!="total" 
    psig["avg",] <- rowMeans(t( psig[only.clusters, ] ),na.rm=T)
    
    return(list(psig=psig, summary=summary.table,
                tables=hyp.tables, images=hyp.img))
           
}

### PLOT cluster-cluster overlaps

#' Wrapper around \code{\link[graphics]{image}} to plot a matrix as
#' it is displayed in R console, i.e. the field \code{dat[1,1]}
#' is at the top left corner. It further allows to plot text
#' into individual fields and have colored axis tick labels.
#' @param dat the numeric data matrix to be plotted
#' @param text a matrix of characters corresponding to \code{dat}
#' which will be plotted on the image
#' @param text.col individual colors for text fields
#' @param axis integer vector, sets whether bottom (1) and/or left
#' (2) axis are draw; the column and row names of \code{dat} will
#' be used as tick labels
#' @param axis1.col invididual colors for x-axis tick labels, length must
#' equal the number of columns of \code{dat}
#' @param axis2.col invididual colors for y-axis tick labels, length must
#' equal the number of rows of \code{dat}
#' @param axis.cex if axis[1|2].col is provided, this sets the tick label size
#' @param ... further arguments to \code{\link[graphics]{image}}, e.g., col
#' to select colors
#' @export
image_matrix <- function(dat, text, text.col, axis=1:2, axis1.col, axis2.col, axis.cex=1.5, ...) {

    ## reverse columns and transpose
    if ( nrow(dat)>1 )
        imgdat <- t(apply(dat, 2, rev))
    else
        imgdat <- t(dat)
    image(x=1:ncol(dat), y=1:nrow(dat), z=imgdat, axes=FALSE, ...)

    ## add text
    if ( !missing(text) ) {
        if ( missing(text.col) )
            text.col <- rep(1, length(c(text)))
        text(x=rep(1:ncol(dat),nrow(dat)), y=rep(nrow(dat):1,each=ncol(dat)),
             paste(t(text)),col=t(text.col))
    }
    
    ## add axes
    ## TODO : handle axes=FALSE
    if ( !missing(axis) ) {
        if ( 1 %in% axis ) 
            if ( !missing(axis1.col) ) # colored ticks
                for ( i in 1:ncol(dat) )
                    axis(1, at=i,colnames(dat)[i],
                         col.axis=axis1.col[i], col=axis1.col[i],
                         las=2, cex.axis=axis.cex, lwd=2)
            else
                axis(1, at=1:ncol(dat), labels=colnames(dat), las=2)
        if ( 2 %in% axis )
            if ( !missing(axis2.col) ) # colored ticks
                for ( i in 1:nrow(dat) )
                        axis(2, at=nrow(dat)-i+1, rownames(dat)[i],
                             col.axis=axis2.col[i], col=axis2.col[i],
                             las=2, cex.axis=axis.cex, lwd=2)
            else
                axis(2, at=nrow(dat):1, rownames(dat),las=2)        
    }
} 


### PLOT CLUSTERED TIME-SERIES


#' calculates cluster averages
#' calculates average values and distributions for each cluster
#' and time point of a time series
#' @param ts a time series, with time points in columns
#' @param cls a clustering of the time series, \code{length(cls)}
#' must equal \code{nrow(cls)}
#' @param cls.srt optional sorting of the clusters
#' @param avg a function for calculating an `average' value for
#' each cluster; default is the \code{median}
#' @param dev a function for calculating deviation from the
#' the average, default is the standard deviation (\code{sd})
#' @param q the \code{\link{quantile}} cut-off, a real number between 0 and 1;
#' \code{1-q} of the time-series are smaller and \code{q}; the range
#' of reported quantile time-series `high' and `low' will cover
#' (1- 2q) of all values in the cluster
#' of the time series are larger then the reported values.
#' @export
getClsAvg <- function(ts, cls, cls.srt,
                      avg=get("median", mode="function"),
                      dev=get("sd", mode="function"), q=.1) {

    if ( missing(cls.srt) )
      cls.srt <- sort(unique(cls))
    
    clavg <- matrix(NA, ncol=ncol(ts), nrow=length(cls.srt))
    rownames(clavg) <- cls.srt
    clstd  <- cllow <- clhig <- clavg

    for ( cl in cls.srt ) {
        clavg[cl,] <- apply(ts[cls==cl,],2, function(x) avg(x,na.rm=T))
        clstd[cl,] <- apply(ts[cls==cl,],2, function(x) dev(x,na.rm=T))
        cllow[cl,]<- apply(ts[cls==cl,],2,
                           function(x) quantile(x,  q,na.rm=T))
        clhig[cl,]<- apply(ts[cls==cl,],2,
                           function(x) quantile(x,1-q,na.rm=T))
        
    }
    list(avg=clavg, std=clstd, low=cllow, high=clhig)
}

#' plots cluster averages
#' plots average time series of clusters as calculated by
#' \code{\link{getClsAvg}}, including the variations around the mean
#' as polygons
#' @param avg cluster time series average object as calculated by
#' \code{\link{getClsAvg}}
#' @param cls.srt optional sorting of clusters; clusters will be plotted
#' in this order, i.e. the last in \code{cls.srt} is plotted last
#' @param cls.col optional coloring of clusters
#' @param each logical value indicating whether to plot all cluster
#' averages on one panel (\code{FALSE}) or each cluster on a separate panel
#' (\code{TRUE})
#' @param polygon logical indicating whether to plot the ranges `high' and
#' `low' in the cluster average object \code{avg}
#' @param xlab x-axis label
#' @param time optional time-points of the x-axis
#' @param ylab y-axis label
#' @param ylim range of the y-axis
#' @export
plotClsAvg <- function(avg, cls.srt, cls.col,
                       each=FALSE, polygon=TRUE,
                       xlab, time, 
                       ylab="average",ylim) {

    if ( missing(time) ) {
        time <- 1:ncol(avg$avg)
        if ( missing(xlab) )
          xlab <- "index"
    } 
    if ( missing(xlab) )
      xlab <- "time"

    if ( missing(cls.srt) )
      cls.srt <- rownames(avg$avg)

    if ( missing(cls.col) ) {
        cls.col <- rainbow(length(cls.srt))
        names(cls.col) <- cls.srt
    }

    ## average and polygon colors
    pol.col <- add.alphas(cls.col,rep(.2,length(cls.col)))
    avg.col <- add.alphas(cls.col,rep(1,length(cls.col)))

    if ( missing(ylim) )
      ylim <- range(avg$avg[cls.srt,])*1

    if ( each ) {
        mai <- par("mai")
        mai[c(1,3)] <- 0
        old.par <- par(mfcol=c(length(cls.srt),1),mai=mai)
    } else {
        plot(1,col=NA,axes=FALSE,
             xlab=xlab,xlim=range(time),
             ylab=ylab,ylim=ylim)
        axis(1);axis(2)
        axis(3, at=time, labels=FALSE)
        mtext("samples", 3, 1.2)
    }
    for ( cl in cls.srt ) {
        if ( each ) {
            plot(1,col=NA,axes=FALSE, xlab=NA,xlim=range(time),
                 ylab=paste(ylab,cl,sep=" - "),ylim=ylim)
            axis(1);axis(2)
        }
        if ( polygon )
          polygon(c(time,rev(time)),c(avg$low[cl,],rev(avg$high[cl,])),
                  col=pol.col[cl],border=NA)
        lines(time, avg$avg[cl,], col=avg.col[cl], lwd=3)
    }
    ## reset plot pars
    if ( each ) par(old.par)
    ##pol.col
}

## util to add/replace alpha values of color vectors
## after http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add.alphas <- function(col, alpha=rep(1,length(col))){
    nms <- names(col)
    lcol <- lapply(col, function(x) c(col2rgb(x)/255))
    col <- sapply(1:length(col), function(x) {
        y <- lcol[[x]]
        rgb(y[1], y[2], y[3], alpha=alpha[x])})
    names(col) <- nms
    col
}
