
### TOOLS for handling feature-based clusterings

## TODO:
## calculate average phases via phaseDist or directly
## auto-select oscillatory clusters from p-values
## clusterCluster/image_matrix : make higher level interface, see
## plot.hypergeoTables in analyzeSegments_2017.R
## integrate GO wrapper for clusterCluster

### COMPARE CLUSTERINGS

#' calculates overlaps between two clusterings
#' 
#' calculates mutual overlaps between two clusterings of the same data set
#' using hypergeometric distribution statistics for significantly
#' enriched or deprived mutual overlaps;
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

    if ( class(cl1)=="clustering" ) {
        K <- paste("K:",cl1$selected,sep="")
        if ( missing(cl1.srt) )
            cl1.srt <- cl1$sorting[[K]]
        cl1 <- cl1$clusters[,K]
    }
    if ( class(cl2)=="clustering" ) {
        K <- paste("K:",cl2$selected,sep="")
        if ( missing(cl2.srt) )
            cl2.srt <- cl2$sorting[[K]]
        cl2 <- cl2$clusters[,K]
    }
   
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
  do.prich <- sum(alternative=="greater")>0
  do.ppoor <- sum(alternative=="less")>0
  if ( sum(alternative=="two.sided")>0 )
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
    jaccard <- percent <- overlap
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

          ## Jaccard:
          intersect <- q # overlap
          union <- k+m-q
          jaccard[i,j] <- intersect/union
      }
  }
  
  result <- list(overlap=overlap,percent=percent,jaccard=jaccard)

  ## append p-values
  #p.value <- prich
  if ( do.prich & !do.ppoor ) p.value <- prich
  if ( do.ppoor & !do.prich ) p.value <- ppoor

  ## take smaller p-value, if both are requested!
  if ( do.ppoor &  do.prich )
    for ( prow in 1:nrow(prich) )
      for ( pcol in 1:ncol(prich) )
        p.value[prow,pcol] <- min(prich[prow,pcol],ppoor[prow,pcol])

    ## TODO: add test number for later bonferroni correction
  result <- append(result, list(p.value=p.value))
    result$alternative <- alternative
  class(result) <- "clusterOverlaps"
  return(result)
}

#' plot cluster-cluster overlaps
#'
#' Plots the significance distribution of a cluster-cluster
#' overlap statistics provided by \code{\link{clusterCluster}},
#' where the a white-black gradient is calculated from \code{-log(p)}.
#' @param x a `clusterOverlaps' object returned by
#' \code{\link{clusterCluster}}
#' @param p.min significance cutoff, p-values equal or smaller to
#' this cutoff will appear black; TODO: plot legend
#' @param n number of gray shades between \code{p=1} (white)
#' and \code{p >= p.min} (black)
#' @param p.txt p-value cutoff for showing white text
#' @param ... arguments to \code{\link{image_matrix}}
## TODO: define as plot.clusterOverlaps method?
## TODO: select white text colors close to p.min !
## TODO: sort by significance?
#' @export
plotOverlaps <- function(x, p.min=0.001, n=100, p.txt=-log2(p.min)*.85, ...) {

    ## set up p-value and colors
    pval <- x$p.value
    pval[pval<=p.min] <- p.min
    pval <- -log2(pval)
    breaks <- seq(0,-log2(p.min),length.out=n+1)
    colors <- grDevices::gray(seq(1,0,length.out=n))
    ## set up text and text colors
    txt <- x$overlap
    txt[txt=="0"] <- ""
    txt.col <- txt
    txt.col[] <- "black"
    txt.col[pval >= -log2(p.txt)] <- "white"

    image_matrix(pval,breaks=breaks,col=colors,axis=1:2,text=txt, text.col=txt.col, ...)
}

#' parse an annotation file (a bidirectional map)
#' 
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


#' clustering enrichment scan
#' 
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
        tmp <- clusterCluster(bins, cls, alternative=c("greater"))
       
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
#' data.matrix `as-is' wrapper for \code{\link[graphics]{image}}
#' 
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

### CLUSTERING OBJECT UTILS

#' get clustering ID
#'
#' get the column name or index of the selected
#' clustering; allows to interface values for the selected clustering
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
    kCol <- paste0("K:", K)
    if ( name )
        return(kCol)
    else
        return(which(colnames(cset$clusters)==kCol))
}

## TODO; finish implementation
## sort clusters by time series phase
##
## takes a time-series and a clustering object,
## calculates phases for each cluster and sorts clusters
## by phase
phasesortClusters <- function(x, cls) {
    cls.srt <- unique(cls)
    names(cls.srt) <- cls.srt
    phase <- calculatePhase(x)
    phase <- sapply(cls.srt, function(x) phaseDist(phase[cls==x,]))
    names(sort(phase[grep("mean",rownames(phase)),]))
}

### PLOT CLUSTERED TIME-SERIES


#' calculates cluster averages
#' 
#' calculates average values and distributions for each cluster
#' and time point of a time series
#' @param ts a matrix of time series, with time points in columns
#' @param cls a clustering of the time series, \code{length(cls)}
#' must equal \code{nrow(cls)}
#' @param cls.srt optional sorting of the clusters
#' @param avg a function for calculating an `average' value for
#' each cluster; default is the \code{median}
## @param dev a function for calculating deviation from the the average,
## default is the standard deviation (‘sd’)
#' @param q either numeric 0-1, the fraction of data for which high
#' and low data cut-offs are calculated, or a function name for
#' calculating variance (eg. "sd", "var"), which will be added and
#' subtracted from the average (argument \code{avg}) for high
#' and low data cut-offs
#' @param rm.inf remove infinite values (e.g. from log transformations)
#' @export
clusterAverages <- function(ts, cls, cls.srt, avg="median", q=.9, rm.inf=TRUE) {

    ## get requested functions
    avg <- get(avg, mode="function")
    #dev <- get(dev, mode="function")
    
    if ( missing(cls.srt) )
      cls.srt <- sort(unique(cls))
    ## ensure present
    cls.srt <- cls.srt[cls.srt%in%cls]

    ## ensure use of rownames
    cls.srt <- as.character(cls.srt)
    
    clavg <- matrix(NA, ncol=ncol(ts), nrow=length(cls.srt))
    rownames(clavg) <- cls.srt
    cllow <- clhig <- clavg

    if ( is.numeric(q) )
        qf <- (1-q)/2
    else qf <- get(q, mode="function")

    ## rm Inf, eg.
    ts[is.infinite(ts)] <- NA
    
    for ( cl in cls.srt ) {
        ## average and deviation
        clavg[cl,] <- apply(ts[cls==cl,,drop=F],2, function(x) avg(x,na.rm=T))
        if ( is.numeric(q) ) {
            ## upper/lower quantiles
            cllow[cl,]<- apply(ts[cls==cl,,drop=F],2,
                               function(x) quantile(x,  qf,na.rm=T))
            clhig[cl,]<- apply(ts[cls==cl,,drop=F],2,
                               function(x) quantile(x,1-qf,na.rm=T))
        } else  {
            df <- apply(ts[cls==cl,,drop=F],2, function(x) qf(x,na.rm=T))
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

#' plot indivividual time-courses
#'
#' plot indivividual time-courses (GOI: genes of interest) from
#' timeseries and clustering object (simple list)
#' @param x "timeseries" object from function
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}}
#' @param cls "clustering" object from function
#' \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param goi list of feature ids (rownames in cls$clusters) to plot
#' @param grep logical, if TRUE \code{goi} are searched by pattern
#' match (as argument \code{pattern} in
#' \code{\link[base:grep]{grep}}) instead of perfect matches; useful
#' eg. if features in \code{cls}
#' @param each plot separate panels for each cluster
#' @param lwd line width of single time-series
#' @param leg.xy position of the legend, see
#' @param ... arguments to \code{\link{plotClusters}}
#' \code{\link[graphics:legend]{legend}}
#' @export
plotSingles <- function(x, cls, goi, grep=FALSE,
                        each, lwd=2, leg.xy="topleft", ...) {
    ## TODO: set all genes
    ## in cls$cluster to "-1"

    ## rm none-present goi
    rm <- !goi%in%rownames(cls$clusters)
    if ( grep )
      rm <- sapply(goi, function(x) length(grep(x, rownames(cls$clusters)))==0)
    if ( sum(rm) )
      cat(paste(paste(goi[rm],collapse=";"),"not found\n"))
    goi <- goi[!rm]
    
    if ( length(goi)==0 )
        stop()

    ## get goi and set all other clusterings to -1
    kp <- which(rownames(cls$clusters)%in%goi)
    if ( grep )
      kp <- sapply(goi, grep, rownames(cls$clusters))
    cls$clusters[-kp,] <- -1
    
    avg <- plotClusters(x, cls, avg.col=NA, lwd=lwd, lwd.avg=0, each=each, alpha=1, use.lty=TRUE, type=c("all"), ...)
    leg <- do.call(rbind,avg$legend)
    if ( !is.null(names(goi)) ) {
        if ( grep )
          leg[,"id"] <- names(goi)[unlist(sapply(goi, grep, leg[,"id"]))]
        else
          leg[,"id"] <- names(goi)[match(leg[,"id"], goi)]
    }
    legend(leg.xy, legend=leg[,"id"], lty=leg[,"lty"], col=leg[,"col"], lwd=lwd,
           bg="#FFFFFFAA",bty="o")
    avg$legend <- leg
    invisible(avg)
}

#' plots cluster averages 
#'
#' as \code{\link{plot.clusteraverages}} but using
#' package \link[segmenTier]{segmenTier}'s time series and clustering object
#' as input; calls, and plots and returns results from
#' \code{\link{clusterAverages}}.
#' @param x either a simple data matrix with time points (x-axis) in columns,
#' or a processed time-series as provided by
#' \code{\link[segmenTier:processTimeseries]{processTimeseries}}
#' @param cls eiter a vector (length(cls)==nrow(x) or a structure of
#' class 'clustering' as returned by segmenTier's
#' \code{\link[segmenTier:clusterTimeseries]{clusterTimeseries}}
#' @param k integer or string specifiying to clustering (K: cluster numbers)
#' to be used if cls is of class 'clustering'; if missing (default) the
#' `selected' clustering from \code{cls} is chosen
#' @param norm normalization of the time-series data, must be a function
#' that transforms the data, available via \link{segmenTools} are
#' \code{lg2r}, \code{ash}, \code{log_1}, \code{meanzero} normalizations
#' @param avg a function for calculating an `average' value for
#' each cluster; default is the \code{median}
#' @param q either numeric 0-1, the fraction of data to be shown,
#' or a function name for calculating variance (eg. "sd", "var")
#' @param cls.srt optional cluster sorting, can be used for selection of
#' subsets of clusters; if cls is of class 'clustering' it is taken
#' from there
#' @param type string specifying the type of plot: "rng" for plotting
#' only data ranges (see argument \code{q}) or "all" to plot
#' each individual time-course (as thin lines)
#' @param cls.col optional named cluster color vector, where names indicate
#' the clusters (as.caracter(cls)); if cls is of class 'clustering' it
#' is taken from there
#' @param each logical value indicating whether to plot all cluster
#' averages on one panel (\code{FALSE}) or each cluster on a separate panel
#' (\code{TRUE})
#' @param time optional numeric vector specifiying x-axis time coordinates
#' @param xlab x-axis label (auto-selected if missing)
#' @param ylab y-axis label (only used if \code{each==FALSE})
#' @param ylim either conventional range of the y-axis, or a string
#' specifying whether ylim should be calculated from the average
#' (\code{ylim="avg"}) or from the lower/upper ranges (\code{ylim="rng"})
#' @param ylim.scale if \code{ylim=="avg"}, the calculated ylim will be
#' extended by this fraction of the total range on both sides
#' @param avg.col color for average line; used only if \code{type="all"}
#' @param lwd.avg line width for average plots (if \code{type=="all"})
#' @param lwd line width for indidiual time series plots (if \code{type=="all"})
#' @param use.lty use individual line type expansion (if \code{type=="all"})
#' @param alpha set alpha value of range or individual time series
#' colors (color opaqueness)
#' @param ... further arguments to \code{\link{plot.clusteraverages}} 
## TODO: clean up mess between plot.clustering, plot.clusteraverages and this
## plot.clusteraverages should become a function of plot.clustering,
## and clusterAverages can be a private function!
## make function `timeseriesPlot' or `clusterPlot', that takes
## either tset/cset or matrix/vector
#' @export
plotClusters <- function(x, cls, k, cls.col, cls.srt, each=TRUE, type="rng", 
                         norm, avg="median",  q=.9, 
                         ylab, ylim=ifelse(each,"avg","rng"), ylim.scale=.1,
                         time, xlab, avg.col="#000000",
                         lwd=.5, lwd.avg=3, use.lty=FALSE, alpha=.2, ...) {

    
    if ( class(cls)=="clustering" ) {

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
    if ( class(x)=="timeseries" )
        ts <- x$ts
    else ts <- data.matrix(x)

    ## normalize
    if ( !missing(norm) )
        ts <- get(norm,mode="function")(ts)
    else norm <- "raw"
    
    ## x-axis
    if ( missing(time) ) {
        time <- 1:ncol(ts)
        if ( missing(xlab) )
            xlab <- "index"
    } 
    if ( missing(xlab) )
        xlab <- "time"

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
            ylim <- range(avg$avg[cls.srt,])
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
        mai[c(1,3)] <- 0
        par(mfcol=c(length(cls.srt),1),mai=mai)
    } else {
        plot(1,col=NA,axes=FALSE,
             xlab=xlab,xlim=range(time),
             ylab=norm,ylim=ylim, ...)
        axis(1);axis(2)
        axis(3, at=time, labels=FALSE)
        mtext("samples", 3, 1.2)
    }
    ## plot each cluster in cls.srt
    used.pars <- rep(list(NA), length(cls.srt))
    names(used.pars) <- cls.srt
    for ( cl in cls.srt ) {
        if ( each ) {
            plot(1,col=NA,axes=FALSE, xlab=NA, xlim=range(time),
                 ylab=paste(cl," (",cls.sze[cl],")",sep=""),ylim=ylim, ...)
            axis(1);axis(2)
        }
        if ( "rng"%in%type ) ## polygon
          polygon(c(time,rev(time)),c(avg$low[cl,],rev(avg$high[cl,])),
                  col=pol.col[cl],border=NA)
        if ( "all"%in%type ) {
            idx <- cls==cl
            if ( use.lty )
                lty <- rep(1:6, len=sum(idx,na.rm=TRUE))
            else  lty <- rep(1, sum(idx,na.rm=TRUE)) #lty <- all.lty[idx]
            matplot(time, t(ts[idx,,drop=F]), add=TRUE,
                    type="l", lty=lty, col=all.col[idx], lwd=lwd)
            ## store
            #if ( use.lty )
            used.pars[[cl]] <- data.frame(id=rownames(ts)[idx],
                                          lty=lty,col=all.col[idx],
                                          stringsAsFactors=FALSE)
        }
        lines(time, avg$avg[cl,], lwd=lwd.avg, col=avg.col[cl]) ## average last
        points(time, avg$avg[cl,], col=avg.col[cl])
    }
    ## reset plot pars
    if ( each ) {
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
#' @param polygon logical indicating whether to plot the ranges `high' and
#' `low' in the cluster average object \code{avg}
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
                                each=FALSE, polygon=TRUE,
                                xlab, time, 
                                ylab="average", ylim=ifelse(each,"avg","rng"),
                                ylim.scale=.1,...) {
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
        if ( polygon ) ## polygon
          polygon(c(time,rev(time)),c(avg$low[cl,],rev(avg$high[cl,])),
                  col=pol.col[cl],border=NA)
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
image_clustering <- function(cset, k=selected(cset), dir="h", ...) {

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
    lcol <- lapply(col, function(x) c(col2rgb(x)/255))
    col <- sapply(1:length(col), function(x) {
        y <- lcol[[x]]
        rgb(y[1], y[2], y[3], alpha=alpha[x])})
    names(col) <- nms
    col
}
