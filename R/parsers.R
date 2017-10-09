
### A VARIETY OF BIO FILE PARSERS

## GEO MICROARRY

#' parse GO Soft archives
#'
#' parses GEO Soft microarray experiment archives into a
#' data matrix and an ID-probe mapping table
#' @param file a GEO Soft archive file (eg. GSE18902_family.soft.gz)
#' @param idcol ID columns from the platform table; NOTE: that the first
#' column in this vector 
#' @param valcol value column (currently only 1 is allowed) in the sample tables
#' @param title if TRUE a sample title will be retrieved
#' @param desc if TRUE, value descriptions in field "#<valcol>" will be retrieved
#' @export
parseGEOSoft <- function(file, idcol="ORF",valcol="VALUE", title=TRUE, desc=TRUE) {

    ## working with IDs, set this to false!
    op <- options(stringsAsFactors=F)
    
    ## parse all lines
    cnx <- gzfile(file)
    lines <- readLines(cnx)
    ## find probe-ID map
    start <- grep("^!platform_table_begin",lines)
    end <- grep("^!platform_table_end",lines)
    ids <- read.delim(cnx,skip=start, nrow=end-start-2,row.names=1)
    ids <- ids[,idcol,drop=FALSE]
    #ids[ids=="",idcol[1]] <- rownames(ids)[ids==""] 
    ## find samples
    idx <- grep("^\\^SAMPLE", lines)
    sampleids <- sub(".*= ","",lines[idx])
    sidx <- grep("^!sample_table_begin", lines)
    eidx <- grep("^!sample_table_end", lines)-1
    ## sample titles
    tit <- NULL
    if ( title ) {
        tcol <- "^!Sample_title = "
        tidx <- grep(tcol, lines)
        tit <- sub(tcol,"",lines[tidx])
        names(tit) <- sampleids
    }
    ## value description
    description <- NULL
    if ( desc ) {
        descol <- paste("^#",valcol," = ",sep="")
        vidx <- grep(descol,lines)
        description <- sub(descol,"",lines[vidx])
        names(description) <- sampleids
    }
    ## get data
    dat <- matrix(NA, nrow=nrow(ids), ncol=length(idx))
    rownames(dat) <- rownames(ids)
    colnames(dat) <- sampleids
    for ( i in 1:length(idx) ) {
        tmp <- read.delim(gzfile(file),
                          skip=sidx[i], nrow=eidx[i]-sidx[i]-1,row.names=1)
        dat[rownames(tmp),i] <- tmp[,valcol]
    }
    list(data=dat, ids=ids, title=tit, description=description)
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
