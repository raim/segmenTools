
### A VARIETY OF BIO FILE PARSERS

## GEO MICROARRY

#' parse GO Soft archives
#'
#' parses GEO Soft microarray experiment archives into a
#' data matrix and an ID-probe mapping table
#' @param file a GEO Soft archive file (eg. GSE18902_family.soft.gz)
#' @param idcol columns to retrieve from the platform table
#' @param valcol value column (currently only 1 is allowed) in the sample tables
#' @param title if TRUE a sample title will be retrieved
#' @param desc if TRUE, value descriptions in field "#<valcol>" will be retrieved
#' @export
parseGEOSoft <- function(file, idcol, valcol="VALUE", title=TRUE, desc=TRUE) {

    ## working with IDs, set this to false!
    old <- unlist(options("stringsAsFactors"))
    op <- options(stringsAsFactors=F)
    
    ## parse all lines
    cnx <- gzfile(file)
    lines <- readLines(cnx)
    ## find probe-ID map
    start <- grep("^!platform_table_begin",lines)
    end <- grep("^!platform_table_end",lines)
    ids <- read.delim(cnx,skip=start, nrow=end-start-2,row.names=1)
    if ( missing(idcol) )
        idcol <- colnames(ids)
    ids <- ids[,idcol,drop=FALSE]
    #ids[ids=="",idcol[1]] <- rownames(ids)[ids==""] 
    ## find samples
    idx <- grep("^\\^SAMPLE", lines)

    
    ## reduce by those with actual data
    cidx <- grep("^!Sample_data_row_count", lines) ## data rows present?
    samplecnt <- as.numeric(sub(".*= ","",lines[cidx]))
    idx <- idx[samplecnt>0]

    ## sample IDs
    sampleids <- sub(".*= ","",lines[idx])
    
    ## sample titles
    tit <- NULL
    if ( title ) {
        tcol <- "^!Sample_title = "
        tidx <- grep(tcol, lines)
        tidx <- tidx[samplecnt>0]
        tit <- sub(tcol,"",lines[tidx])
        names(tit) <- sampleids
    }

    ## data description - should be only present for samples with data
    description <- NULL
    if ( desc ) {
        descol <- paste("^#",valcol," = ",sep="")
        vidx <- grep(descol,lines)
        description <- sub(descol,"",lines[vidx])
        names(description) <- sampleids
    }

    ## data tables
    sidx <- grep("^!sample_table_begin", lines)
    eidx <- grep("^!sample_table_end", lines)-1
    dat <- matrix(NA, nrow=nrow(ids), ncol=length(idx))
    rownames(dat) <- rownames(ids)
    colnames(dat) <- sampleids
    for ( i in 1:length(idx) ) {
        tmp <- read.delim(gzfile(file),
                          skip=sidx[i], nrow=eidx[i]-sidx[i]-1,row.names=1)
        dat[rownames(tmp),i] <- tmp[,valcol]
    }
    res <- list(data=dat, ids=ids, title=tit, description=description)
    class(res) <- "geosoft"
    # reset option
    options(stringsAsFactors=old)
    res
}

#' summarize GEOSoft probes
#'
#' takes as input the list returned by \code{\link{parseGEOSoft}},
#' a GEO microarray data set, including a table of probe mappings,
#' as provided via GEO soft family files.
#' Argument \code{id} specifies a probe mapping column ID.
#' Multiples probes for the same features (same string in column ID) will
#' then be summarized by the function in argument \code{avg}
#' (default: median). If the package \pkg{farms} installed, the more
#' sophisticated probe summarization of this package can be used
#' \url{https://www.bioconductor.org/packages/release/bioc/html/farms.html}.
#' @param data a list as returned by function \code{\link{parseGEOSoft}}
#' @param id a column ID; probes with equal strings in this column will
#' be summarized
#' @param avg a function (or the name of a function as a string)
#' for calculating an `average' value for each cluster; default is
#' the \code{mean}
#' @param farms logical to indicate to use the bioconductor package
#' \pkg{farms} for probe summarization (function \code{\link[farms:generateExprVal.method.farms]{generateExprVal.method.farms}})
#' @param keep.empty keep columns with no entry in column ID by
#' copying the probe ID
#' @param replicate if TRUE probe sets assigned to multiple features
#' (separated by string provided in rep.sep) will be duplicated to
#' and each feature will get its own row
#' @param repsep separator for multiple-feature probes
#' @param verb print messages
#' @param ... arguments to the function specified in \code{avg} or for
#' function \code{\link[farms:generateExprVal.method.farms]{generateExprVal.method.farms}} if \code{farms==TRUE}
#' @export
summarizeGEOSoft <- function(data, id="ORF", avg="mean", farms=FALSE,
                             keep.empty=FALSE, replicate=FALSE, repsep=";",
                             verb=TRUE, ...) {
    ## get summarization function
    if ( mode(avg)!="function" )
        avgf <- get(avg, mode="function")
    else avgf <- avg
    
    ## get data
    dat <- data$data

    ## summarize multiple probes per gene
    ids <- data$ids

    ## 1) handle empty strings
    ## fill empty names
    empty <- ids[,id]==""
    if ( keep.empty ) {
        if ( verb )
            cat(paste("adding", sum(empty),
                      "features with empty field in column", id, "\n"))
        ids[empty,id] <- rownames(ids)[empty]
    } else if ( verb ) {
        cat(paste("discarding", sum(empty),
                  "features with empty field in column", id, "\n"))
        ids <- ids[!empty,,drop=FALSE]
        dat <- dat[!empty,,drop=FALSE]
    }

    ## 2) handle duplicates
    ## find duplicates
    uids <- ids[,id]
    dups <- duplicated(uids)
    dids <- unique(uids[dups])
    if ( verb )
        cat(paste("mapping", nrow(ids), "probes to",
                  sum(!dups), "features in column", id, "\n"))
    
    ## summarize duplicates (into first row where it occurs)
    ## TODO: for farms; should data be mean/median-centered first?
    idx <- sapply(dids, function(x) which(ids[,id]%in%x))
    for ( i in idx ) {
        if ( farms ) ## USING FARMS
            dat[i[1],] <- 2^farms::generateExprVal.method.farms(dat[i,],...)$exprs
        else ## or a simple averaging function
            dat[i[1],] <- apply(dat[i,],2,avgf,...)
    }
    ## and rm duplicate rows
    dat <- dat[!dups,]
    rownames(dat) <- uids[!dups]

    ## 3) handle duplicate features
    if ( replicate ) {
        uids <- rownames(dat)
        idx <- grep(repsep, uids)
        fts <- strsplit(rownames(dat), repsep)
        cat(paste("THIS FEATURE IS NOT YET AVAILABE;",
                  "set `replicate to FALSE\n"))
    }

    
    ## replace data in original list
    data$data <- dat # replace
    ## add info
    data$mappedto <- id
    if ( farms ) avg <- "farms::generateExprVal.method.farms"
    data$avg_function <- list(name=avg,
                             argument=paste(names(list(...)),"=",
                                            deparse(substitute(...))))
    data$empty <- sum(empty)
    data
}

## parse GFF file 
## 
## parse a GFF genome annotation file, credit: from davidTiling package via
## \url{https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html}
## @param gffFile gff file name
## @param nrows number of rows to parse
readGFF <- function(gffFile, nrows = -1) {
    cat("Reading ", gffFile, ": ", sep="")
    gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
      header=FALSE, comment.char="#", nrows = nrows,
      colClasses=c("character", "character", "character",
        "integer", "integer",
        "character", "character", "character", "character"),
      stringsAsFactors=TRUE)
    colnames(gff) = c("seqname", "source", "feature", "start", "end",
              "score", "strand", "frame", "attributes")
    cat("found", nrow(gff), "rows with classes:",
        paste(sapply(gff, class), collapse=", "), "\n")
    stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
    return(gff)
 }

#' creates a gff3-like table
#'
#' Creates a table with columns as present in gff3 genome
#' annotation files, including ;-separated attribute lists.
#' Columns with RGB colors can be used to create \code{snapgene}-specific
#' note in attributes.
#' @param tab input table of genomic features with chromosome coordinate
#' information
#' @param columns a named string vector mapping from columns in
#' \code{tab} (values) to the first 7 columns required for
#' the gff3 file (names)
#' @param attributes a named string vector mapping  from columns in
#' \code{tab} (values) to values in the attribute list; the names of
#' the vector will become the key in the attribute list
#' ("<key> = <value>") and attributes will be separated by
#' arugment \code{sep}
#' @param sep separator for attribute list
## TODO: goal is to create a gff file that can be converted
## to snapgene-genbank with gff_to_genbank.py
## note: for gff_to_genbank.py strand must be in +/-
#'@export
tab2gff <- function(tab,
                    columns=c(seqid="chr", "source"="source", type="type",
                              start="start", end="end", score="score",
                              strand="strand",phase="phase"),
                    attributes=c(ID="ID",Name="name",Alias="alias",
                                 Parent="parent",color="color"),
                    sep=";") {
    miscol <- columns[!columns%in%colnames(tab)]
    cols <- columns[columns%in%colnames(tab)]
    out <- tab[,cols]
    colnames(out) <- names(cols)
    ## TODO: add frame/phase column
    ## TODO: add score column
    ## TODO: add translation?
    ## TODO: parse attributes into
    ## 1) parse known attributes
    ## 2) collect unknown attributes
    ## 3) parse color and convert to snapgene "note"
    out
}

#' parse a GFF3 file into a table
#' 
#' parse gff3 files into tables, including conversion of all attributes
#' @param file gff file name
#' @param attrsep main separator of attribute string fields
#' @param fieldsep separator for attribute <id>=<value> couples
#' @export
gff2tab <- function(file, attrsep=";", fieldsep="=") {

    gff <- readGFF(file)

    ## parse attributes into table
    att <- strsplit(gff[,"attributes"], split = attrsep, fixed = TRUE)
 
    ## get all fields (attributes with = : strsplit returns 2 fields)
    fields <- unique(c(unlist(sapply(att, function(atts) {
        a <- strsplit(atts, split = fieldsep, fixed = TRUE)
        unlist(lapply(a, function(x) if (length(x)==2) x[1]))
    }))))
    attm <- matrix(NA, nrow=nrow(gff), ncol=length(fields))
    colnames(attm) <- fields
        
    
    ## fuse back fields without fieldsep
    att2 <- lapply(att, function(x) {
        a <- strsplit(x, split = fieldsep, fixed = TRUE)
        ## fuse back fields without =
        newa <- list()
        if ( length(a)>0 )
          for ( j in 1:length(a) ) {
              len <- length(newa)
              ## append if length is ok (= present)
              if ( length(a[[j]])==2 ) newa <- append(newa, list(a[[j]]))
              if ( length(a[[j]])==1 ) # or append to previous content
                newa[[len]][2] <- paste(newa[[len]][2],a[[j]],sep=";")
              if ( !length(a[[j]])%in%1:2) cat(paste("error\n"))
          }
        newa})
    for ( i in 1:nrow(attm) ) {
        key <- unlist(lapply(att2[[i]], function(y) y[1]))
        val <- unlist(lapply(att2[[i]], function(y) y[2]))
        attm[i,key] <- val
    }

    cbind.data.frame(gff[,colnames(gff)!="attributes"],attm,
                     stringsAsFactors=FALSE)
}

