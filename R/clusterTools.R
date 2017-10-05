
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
          
          q <- sum(cl1==f1[i] & cl2==f2[j]) # white balls drawn
          k <- sum(cl2==f2[j]) # number of balls drawn
          
          overlap[i,j] <- q
          percent[i,j] <- round(100*q/k, digits = 2)
          
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

    ## TODO: add test number for later bonferroni correction
  result <- append(result, list(p.value=p.value))
  class(result) <- "clusterOverlaps"
  return(result)
}

#' parse an annotation file (a bidirectional map)
#' parses a bidirectional map of feature IDs vs. annotation terms, e.g.
#' the GO annotation file at \url{ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt.gz}
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


#' Clustering enrichment scan
#' Scans for overlap enrichments of a clustering in a matrix of
#' clusterings, potentially simply a TRUE/FALSE table, eg. indicating
#' annotation with a specific Gene Ontology or other term. It reports
#' tables of all overlap sizes and their enrichment p-values. The
#' function is a wrapper around \code{\link{clusterCluster}}, which
#' performs cumulative hypergeometric distribution tests between
#' two clusterings.
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
#' @param bin.filter string, indicating bins (categories in \code{data})
#' to be globally omitted; useful eg. for logical data to omit
#' all enrichments with category "FALSE" (indicating deprivement)
#' @param verbose print progress messages
#' @export
clusterAnnotation <- function(cls, data, p=1,
                              cls.srt, terms=NULL, bin.filter, verbose=TRUE) {

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
    ## TODO : parallelize this
    if (verbose) cat(paste("HYPERGEOMETRIC STATISTICS FOR:\n",
                           ncol(data), "categories\t:\n"))
    for ( j in 1:ncol(data) ) {
        bins <- data[,j]
        name <- colnames(data)[j]
        expected <- unique(bins)


        if (verbose) cat(paste(j, "-",name, ", ", sep=""))
        
        ## cumulative hypergeometric distribution
        ## TODO: take sum of bonferroni correction factors
        ## correct below in bin.filter
        tmp <- clusterCluster(bins, cls, req.vals=c("greater"))
       
        # get overlap and p.values in original bin order
        if ( !is.null(tmp) ) {
            ovl <- tmp$overlap # contingency table
            pvl <- tmp$p.value
            overlap <- rbind(overlap, ovl)
            pvalues <- rbind(pvalues, pvl)
            # copy name in same size to allow easy cbind below
            hyp.names<-c(hyp.names, rep(name, nrow(ovl)))
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

    if (verbose) cat(paste(" ... done;\n"))

    ## cluster averages in summary table
    only.clusters <- rownames(psig)!="total" 
    psig["avg",] <- rowMeans(t( psig[only.clusters, ] ),na.rm=T)

    ## FILTER BY P-VALUE, BINS and add DESCRIPTIONS
    ## TODO: add GO classes (MF,CC,BP)
    ## filter significant and order by p-value
    if ( p<1 )
        sig <- lapply(hyp.tables, function(x) {
            x <- x[x[,"p-value"] <= p,] #smaller then passed p-value threshold?
            x[order(x[,"p-value"]),]})
    ## filter by bins
    if ( !missing(bin.filter) )
        sig <- lapply(sig, function(x) x[!x[,"bin"]%in%bin.filter,] )
        
    ## add description of terms
    if ( !is.null(terms) ) # not required for other columns (no keys)
        sig <- lapply(sig, function(x)
            cbind.data.frame(description=terms[x[,"category"]],x))
    
    return(list(tables=sig, psig=psig))
           
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
#' @param ts a matrix of time series, with time points in columns
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
clusterAverages <- function(ts, cls, cls.srt,
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
    res <- list(avg=clavg, std=clstd, low=cllow, high=clhig)
    class(res) <- "clusteraverages"
    res
}

#' plots cluster averages
#' plots average time series of clusters as calculated by
#' \code{\link{clusterAverages}}, including the variations around the mean
#' as polygons
#' @param x cluster time series average object as calculated by
#' \code{\link{clusterAverages}}
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
#' @param ylim range of the y-axis, will be calculated from the average
#' values, extended \code{ylim.scale}
#' @param ylim.scale if ylim is missing, the calculated ylim will be
#' extended by this fraction of the total range on both sides
#' @param ... arguments to plot
#' @export
plot.clusteraverages <- function(x, cls.srt, cls.col,
                                each=FALSE, polygon=TRUE,
                                xlab, time, 
                                ylab="average",ylim,ylim.scale=.1,...) {
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
    pol.col <- add.alphas(cls.col,rep(.2,length(cls.col)))
    avg.col <- add.alphas(cls.col,rep(1,length(cls.col)))

    ## calculate ylim
    ## TODO: use ranges?
    if ( missing(ylim) ) {
      ylim <- range(avg$avg[cls.srt,])
      ylim <- c(ylim[1]-diff(ylim)*ylim.scale,
                ylim[2]+diff(ylim)*ylim.scale)
    }

    ## plot each cluster on separate panel?
    if ( each ) {
        mai <- par("mai")
        mai[c(1,3)] <- 0
        old.par <- par(mfcol=c(length(cls.srt),1),mai=mai)
    } else {
        plot(1,col=NA,axes=FALSE,
             xlab=xlab,xlim=range(time),
             ylab=ylab,ylim=ylim, ...)
        axis(1);axis(2)
        axis(3, at=time, labels=FALSE)
        mtext("samples", 3, 1.2)
    }
    ## plot each cluster in cls.srt
    for ( cl in cls.srt ) {
        if ( each ) {
            plot(1,col=NA,axes=FALSE, xlab=NA,xlim=range(time),
                 ylab=paste(ylab,cl,sep=" - "),ylim=ylim, ...)
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
