
### A VARIETY OF BIO FILE PARSERS

## GEO MICROARRY

#' parse GO Soft archives
#'
#' Parses GEO Soft microarray experiment archives into a
#' data matrix and an ID-probe mapping table. Probes can be further
#' summarized for features with \code{\link{summarizeGEOSoft}}.
#' @param file a GEO Soft archive file (eg. GSE18902_family.soft.gz)
#' @param idcol columns to retrieve from the platform table
#' @param exp number of the platform to retrieve from a GEO file
#' with multiple platforms. TODO: allow to parse
#' multiple platforms.
#' @param only.info only parse information on probes and experiments, but
#' skip the actual data
#' @param only.data skip all samples without data
#' @param valcol value column (currently only 1 is allowed) in the sample tables
#' @param title if TRUE a sample title will be retrieved
#' @param desc if TRUE, value descriptions in field "#<valcol>" will be retrieved
#' @seealso \code{\link{summarizeGEOSoft}}
#' @export
parseGEOSoft <- function(file, idcol, exp=1, only.info=FALSE,
                         only.data=TRUE,
                         valcol="VALUE", title=TRUE, desc=TRUE) {


    ## working with IDs, set this to false!
    old <- unlist(options("stringsAsFactors"))
    op <- options(stringsAsFactors=FALSE)
    
    ## parse all lines
    ## TODO: grep a  max line in file to avoid reading all lines, or
    ## avoid re-reading the file multiple times below
    cnx <- gzfile(file)
    lines <- readLines(cnx)
    ## find probe-ID map
    start <- grep("^!platform_table_begin",lines)
    end <- grep("^!platform_table_end",lines)
    if ( length(start)!=length(end) )
        stop("different number of platform_table_begin and _end tags")

    ids <- NULL
    skip.multi <- FALSE
    if ( length(start)>0 & length(end)>0 ) {

        if ( length(start)>1 ) {
            ## TODO: allow loop
            cat(paste("taking experiment", exp, "of", length(start), "\n"))
            start <- start[exp]
            end <- end[exp]

            ### TODO: skipe these from parsing below!!
            ## skip multiple tables also below, since we are not sure
            ## to which IDs they refer
            skip.multi <- TRUE
        }

        ## parse data table
        ids <- read.delim(cnx, skip=start, nrow=end-start-2, row.names=1)

        if ( missing(idcol) )
          idcol <- colnames(ids)
        ids <- ids[,idcol,drop=FALSE]
        ##ids[ids=="",idcol[1]] <- rownames(ids)[ids==""]
    } else cat(paste("no data!\n"))

    ## find samples
    idx <- grep("^\\^SAMPLE", lines)

    ## TODO: reduce to a subset, e.g. only the platform
    ## which IDs were parsed above!
    
    ## reduce by those with actual data
    cidx <- grep("^!Sample_data_row_count", lines) ## data rows present?
    samplecnt <- as.numeric(sub(".*= ","",lines[cidx]))
    if ( only.data ) 
      idx <- idx[samplecnt>0]

    ## sample IDs
    sampleids <- sub(".*= ","",lines[idx])
    
    ## sample titles
    tit <- NULL
    if ( title ) {
        tcol <- "^!Sample_title = "
        tidx <- grep(tcol, lines)
        if ( only.data ) 
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
    if ( !is.null(ids) & !only.info ) {
        dat <- matrix(NA, nrow=nrow(ids), ncol=length(idx))
        rownames(dat) <- rownames(ids)
        colnames(dat) <- sampleids
        cat(paste("parsing",length(idx), "sample series: "))
        for ( i in 1:length(idx) ) {

            cat(paste(i,","))
            tmp <- read.delim(gzfile(file),
                              skip=sidx[i], nrow=eidx[i]-sidx[i]-1,row.names=1)

            ## filter rows which are not in ids!
            miss <- which(!rownames(tmp)%in%rownames(dat))
            if ( length(miss)>0 ) {
                warning("data set ", i, ": ",
                        length(miss), " data rows not in IDs")
                tmp <- tmp[-miss,]
            }
            dat[rownames(tmp),i] <- tmp[,valcol]
        }
        cat(paste("done, parsed", nrow(dat), "data rows\n"))
      
        res <- list(data=dat, ids=ids, title=tit, description=description)
        class(res) <- "geosoft"
    } else {
        if ( !is.null(ids) )
            res <- list(ids=ids, samples=sampleids,
                        title=tit, description=description)
        else 
            res <- list(samples=sampleids, title=tit, description=description)
        class(res) <- "geosoft_info"
    }
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
#' @seealso \code{\link{parseGEOSoft}}
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
    idx <- lapply(dids, function(x) which(ids[,id]%in%x))
    names(idx) <- dids
    for ( i in idx ) {
        if ( farms ) ## USING FARMS
            dat[i[1],] <- 2^farms::generateExprVal.method.farms(dat[i,],
                                                                ...)$exprs
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
## Parses a GFF genome annotation file, credit: from davidTiling package via
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
#' @param notes as \code{arguments} but output will be
#' "note=<key>=<value>
#' @param sep separator for attribute list
#' @param source identifier of main source feature (chromosome)
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
                    notes,
                    sep="; ",
                    source=c("chromosome","plasmid")) {
    miscol <- columns[!columns%in%colnames(tab)]
    cols <- columns[columns%in%colnames(tab)]
    out <- as.data.frame(tab[,cols],stringsAsFactors = FALSE)
    colnames(out) <- names(cols)
    ## TODO: add frame/phase column
    ## TODO: add score column
    if ( !"score"%in%colnames(out) )
        out <- cbind.data.frame(out, score=rep(".", nrow(out)),
                                stringsAsFactors = FALSE)
    if ( !"phase"%in%colnames(out) ) {
        out <- cbind.data.frame(out, phase=rep(".", nrow(out)),
                                stringsAsFactors = FALSE)
        out[out[,"type"]=="CDS","phase"] <- "0"
    }
    ## TODO: parse attributes into
    ## 1) parse known attributes
    ## 2) collect unknown attributes
    ## 3) parse color and convert to snapgene "note"

    ## tag chromosomes as source
    idx <- out[,"type"]%in%source
    out[idx,"type"] <- "source"

    attributes <- attributes[attributes%in%colnames(tab)]
    atts <- matrix(NA,nrow(out),ncol=length(attributes))
    for ( i in 1:length(attributes) ) {
        colid <- attributes[i]
        colnm <- names(attributes)[i]
        isna <- is.na(tab[,colid])
        atts[!isna,i] <- paste(colnm,tab[!isna,colid],sep="=")
    }
    attl <- apply(atts, 1, function(x) paste(x[!is.na(x)],collapse=sep))
    nots <- NULL
    if ( !missing(notes) ) {
        nots <- matrix(NA,nrow(out),ncol=length(notes))
        for ( i in 1:length(notes) ) {
            colid <- notes[i]
            colnm <- names(notes)[i]
            isna <- is.na(tab[,colid])
            nots[!isna,i] <- paste0("note=\"",colnm,":",tab[!isna,colid], "\"")
        }
        notl <- apply(nots, 1, function(x) paste(x[!is.na(x)],collapse=sep))
        attl[notl!=""] <- paste(attl[notl!=""],notl[notl!=""],sep=sep)
    }
     
    out <- cbind.data.frame(out, attributes=attl,
                            stringsAsFactors = FALSE)
    ## final sort
    out <- out[,c(names(columns),"attributes")]
    out
}

#' parse a .bed format file
#'
#' UNTESTED: Parses a bed file with 5 (or more) columns into the
#' segmenTools genomic interval format.
#' @param file a file in bed format
#' @param header header names for the parsed bed file
#' @export
bed2coor <- function(file, header=c("chr","start","end","name","score")) {

    dat <- data.table::fread(file,  header=FALSE)
    ## todo: add column names if </> 5 are present
    colnames(dat) <- header

    ## 1-based starts
    dat$start <- dat$start+1

    ## numeric chromosomes
    dat$chr <- sub("chr","",dat$chr)
    if ( all(!is.na(suppressWarnings(as.numeric(dat$chr)))) )
        dat$chr <- as.numeric(dat$chr)
    
    ## numeric strand
    if ( "strand" %in% colnames(dat) ){
        str <- c("+"=1, "-"=-1)
        dat$strand <- str[as.character(dat$strand)]
    }
    
    dat
}

#' parse a GFF3 file into a table
#' 
#' Parses gff3 files into tables, including conversion of all attributes
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


## from library("Biostrings")
#' fasta sequence file parser adapted from (an older version of)
#' \code{Biostrings}
#' @param file file with sequence(s) in fasta format
#' @param checkComments remove information after a semi-colon ; in fasta
#' header
#' @param strip.descs split header at a pattern in argument \code{split.desc}
#' and store in returned structure as \code{details}
#' @param split.desc pattern to split header into description and details
#' with argument \code{strip.descs}
#' @param grepID use ID in description field as names of sequence list.
#' @export
readFASTA <- function (file, checkComments=TRUE, strip.descs=TRUE,
                       split.desc="", grepID=FALSE) 
{    
    if (is.character(file)) {
        if ( length(grep("\\.gz$",file)) )
          file <- gzfile(file)
        else
          file <- file(file, "r")
        on.exit(close(file))
    }
    else {
        if (!inherits(file, "connection")) 
            stop("'file' must be a character string or connection")
        if (!isOpen(file)) {
            open(file, "r")
            on.exit(close(file))
        }
    }
    s1 <- scan(file = file, what = "", sep = "\n", quote = "", 
        allowEscapes = FALSE, quiet = TRUE)
    if (checkComments) {
        comments <- grep("^;", s1)
        if (length(comments) > 0) 
            s1 <- s1[-comments]
    }
    descriptions <- which(substr(s1, 1L, 1L) == ">")
    numF <- length(descriptions)
    if (numF == 0) 
        stop("no FASTA sequences found")
    dp <- descriptions + 1L
    dm <- descriptions - 1L
    end <- c(dm[-1], length(s1))
    fas <- lapply(seq_len(numF), function(i) {
        desc <- s1[descriptions[i]]
        if (strip.descs) 
            desc <- substr(desc, 2L, nchar(desc))
        details <- NULL
        if ( split.desc != "" ) {
          desc <- unlist(strsplit(desc,split.desc))
          if ( length(desc)>2 ) stop("seq id not compatible with split.desc")
          details <- desc[2]
          desc <- desc[1]
        }
        if (end[i] >= dp[i]) {
            seq <- paste(s1[dp[i]:end[i]], collapse = "")
        }
        else {
            warning("record \"", desc, "\" contains no sequence")
            seq <- ""
        }
        ## replace trailing // in sequence
        seq <- gsub("/","",seq)
        list(desc = desc, details=details, seq = seq)
    })
    if ( grepID ) {
        
        desc <- lapply(fas, function(x) {unlist(strsplit(trimws(x$desc), " "))})
        ids <- unlist(lapply(desc, function(x) x[1]))
        ids <- sub(",$","", ids)
        names(fas) <- ids
    }
    fas
}
