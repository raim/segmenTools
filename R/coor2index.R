## GENOME UTILS


#' Tag duplicate names by increasing numbers.
#' @param names a vector of characters
#' @export
tagDuplicates <- function(names) {
    sel <- paste(names,".1",sep="")
    cnt <- 2
    while( sum(duplicated(sel)) ) {
        sel[duplicated(sel)] <- sub("\\..*",paste(".",cnt,sep=""),
                                    sel[duplicated(sel)])
        cnt <- cnt+1
    }
    sub("\\.1$","",sel)
}

#' Generate chromosome index \code{chrS} from lengths
#' @param chrL an ordered vector of chromosome lengths; where the
#' order must correspond to chromosome numbering in feature tables
#' for which chrS is used
#' @export
getChrSum <- function(chrL) c(0,cumsum(chrL))

## util to insert rows, after suggestion by user Ari B. Friedman at
## \url{https://stackoverflow.com/a/11562428}
insertRow <- function(existingDF, newrow, r) {
    if ( r==nrow(existingDF)+1 ) # insert as last row?
        existingDF <- rbind(existingDF, newrow)
    else if ( r<=nrow(existingDF) ) { # insert in between
        existingDF <- as.data.frame(existingDF,stringsAsFactors=FALSE)
        idx <- seq(r, nrow(existingDF)) # shift all by one below r
        existingDF[idx+1,] <- existingDF[idx,] 
        existingDF[r,] <- newrow
    } else
        stop("wrong index, can't be >nrow(<existing data.frame>)")
    existingDF
}
#' insert rows as specified positions
#'
#' Util to insert multiple rows at specified positions
#' in a \code{data.frame}, expanding single-row code by user
#' Ari B. Friedman at
#' \url{https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appendedlooping through new rows}
#' @param existingDF existing \code{data.frame}
#' @param newrows rows to add to \code{existingDF}
#' @param r positions in the existing data.frame at which rows are to
#' be inserted; \code{length(r)} must equal \code{nrow(newrows)}, and
#' all indices \code{r<=nrow(existingDF)+1}.
#' @export
insertRows <- function(existingDF, newrows, r ) {
    ## check that r is sorted and all <= nrow(existingDF)
    r <- sort(r) # SORT!
    if ( any(r>nrow(existingDF)+1) )
        stop("row indices must refer to existing data.frame and",
             " can not be >nrow(<existing data.frame>)+1")
    new <- existingDF
    for ( i in 1:nrow(newrows) )
        new <- insertRow(new, newrows[i,], r[i]+i-1)
    new
}

#' Splits genome features spanning annotated ends of circular chromosomes.
#' 
#' Splits genome features that span start/end coordinates of circular
#' chromosomes, and adds the downstream half with optional modification
#' of ID, and type values. Circular features are recognized here by
#' start > end, in left->right direction of genome annotation.
#' Strand information MUST NOT BE ENCODED IN start/end coordinate direction,
#' but explicitly provided via a strand column!
#' Note that only the upstream half retains
#' all column information (exceptions: see argument \code{copyCols}),
#' the downstream half will only carry information on coordinates, and
#' optionally updated feature type and ID.
#' The update will only happen if the passed table contains type and ID
#' information (see argument \code{idCols}. The split can be reversed
#' by function \code{removeCircularFeatures}.
#' @param features a list of genomic features with coordinates
#' @param chrL obligatory list of chromosome lengths, in order used
#' in chromosome column in \code{features} (see argument \code{coorCols}
#' @param coorCols ordered string vector providing the column names
#' of coordinate columns to be used; must be of length 4 and provide in
#' order: chromosome number (refering to argument \code{chrL}), start, end,
#' and strand (see argument \code{reverse})
#' @param reverse allowed indicators of reverse strand features
#' in strand column (see argument \code{coorCols})
#' @param idTag tag to add to downstream ID and type
#' @param idCols named vector of column names for feature ID, type,
#' and feature parent; note that a "parent" column will be added if not present
#' to refer the downstream half to its upstream feature, which retains
#' all other information
#' @param copyCols copy values to circular feature copy; either logical
#' \code{TRUE} to copy all columns, or a vector of column indices or names
#' @param insertRows insert the circular features below their parent,
#' if set to \code{FALSE} circular features will just be appended; this
#' saves a lot of time for large datasets
#' @seealso \code{removeCircularFeatures}
#' @export
expandCircularFeatures <- function(features, chrL, 
                                   coorCols=c("chr","start","end","strand"),
                                   reverse=c("-",-1),
                                   idTag="-circ2", idCols=c(ID="ID",type="type",
                                                            parent="parent"),
                                   copyCols=FALSE,
                                   insertRows=TRUE) {

    ## chromosome index - revert from chrL
    
    ## add parent column if not present
    if ( idCols["ID"]%in%colnames(features) &
        !idCols["parent"] %in% colnames(features) ) {
        features <- cbind(features,parent=rep(NA,nrow(features)))
    }

    ## filter
    if ( typeof(copyCols)=="logical" ) {
        if ( copyCols )
          copyCols <- colnames(features)
    } else if ( typeof(copyCols)=="integer" )
      copyCols <- colnames(features)[copyCols]
    
    ## get all coordinates
    start <- features[,coorCols[2]] # "start"
    end <- features[,coorCols[3]] # "end"
    strand <- features[,coorCols[4]] # "strand"
    rev <- strand%in%reverse
    
    ## get circular
    circ <- start>end # ASSUMES ORDERED START/END
    cidx <- which(circ) # index in original

    if ( sum(circ)==0 )
      return(features)
    
    ## copy and rename (ID_circ#, type: type_circular, parent: ID)
    cfeat <- as.data.frame(matrix(NA,ncol=ncol(features),nrow=length(cidx)))
    colnames(cfeat) <- colnames(features)
    ## copy requested columns
    cfeat[,copyCols] <- features[cidx,copyCols]
    ## set up type
    if ( idCols["ID"]%in%colnames(features) ) {
        cfeat[,idCols["parent"]] <- features[cidx,idCols["ID"]]
        cfeat[,idCols["ID"]] <- paste(features[cidx,idCols["ID"]],idTag,sep="")
    }
    if ( idCols["type"]%in%colnames(features) )
      cfeat[,idCols["type"]] <- paste(features[cidx,idCols["type"]],
                                      idTag,sep="")
    
      
    crev <- rev[cidx]

    ## set up coordinates
    ## c("chr","start","end","strand")
    cfeat[,coorCols] <- features[cidx,coorCols]

    ## reverse coordinates
    ## copy: end becomes chromosome length
    cfeat[crev,coorCols[3]] <- chrL[cfeat[crev,coorCols[1]]]
    ## original: start becomes 1
    features[circ&rev,coorCols[2]] <- 1

    ## forward coordinates
    ## copy: start becomes 1
    cfeat[!crev,coorCols[2]] <- 1
    ## original: end becomes chromosome length
    features[circ&!rev,coorCols[3]] <- chrL[features[circ&!rev,coorCols[1]]]

    ## insert below original & return
    ## TODO: find faster version via ID mapping!
    ## TODO: without insertRows table seems to have an empty first line!?
    if ( insertRows ) 
        res <- insertRows(features,cfeat,cidx+1)
    else
        res <- rbind(features,cfeat)
    res
}

#' NOT WORKING - Undo \code{expandCircularFeatures}
#' searches for circular features by the \code{idTag} added
#' to ID and type columns in \code{expandCircularFeatures},
#' and maps downstream coordinates back to original features.
#' @param features  list of genomic features with coordinates
#' @param coorCols ordered string vector providing the column names
#' of coordinate columns to be used; must be of length 4 and provide in
#' order: chromosome number (refering to argument \code{chrL}), start, end,
#' and strand 
#' @param idTag tag used for tagging downstream halves
#' @param idCols named vector of column names for feature ID, type,
#' and feature parent; note that a "parent" column will be removed if
#' it is (a) empty and (b) argument \code{rmParent==TRUE}
#' @param rmParent rm the parent column
#' @seealso \code{expandCircularFeatures}
#' @export
removeCircularFeatures <- function(features,
                                   coorCols=c("chr","start","end","strand"),
                                   idTag="-circ2",
                                   idCols=c(ID="ID",type="type",
                                     parent="parent"),
                                   rmParent=TRUE) {
    idCols <- idCols[idCols%in%colnames(features)]
    if ( length(idCols)==0 )
        stop("no columns present to scan for idTag, use argument idCols")
    cidx <- grep(idTag, features[,idCols[1]])
    
}

#' convert chromosome coordinates to continuous index
#' @param features a table of chromosome features that must contain
#' the chromosome number (option \code{chrCol}), one or more chromosome
#' positions (option \code{cols}) and strand information (column
#'  \code{strandCol}).
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information; simply the cumulative lengths of ordered chrosomes,
#' see function \code{\link{getChrSum}}
#' @param chrMap a vector of chromosome names using \code{features}' chromosome
#' column, in the same order as \code{chrS}
#' @param cols name of the columns giving coordinates that will be mapped
#' to continuous index
#' @param chrCol name of the column that gives the chromosome number
#' @param strandCol name of the column that gives forward/reverse strand
#' information
#' @param reverse a vector of possible reverse strand indicators
#' @param circular suppresses re-sorting to start < end for circular chromosomes
#' @export
coor2index <- function(features, chrS, chrMap,
                       cols=c("start","end","coor"),
                       chrCol="chr", strandCol="strand",
                       reverse=c("-",-1), circular=FALSE) {

    ## coordinate columns
    cols <- cols[cols%in%colnames(features)]
    ## strand column - if not present, infer from start>end
    if ( strandCol%in%colnames(features) ) {
        strand <- as.character(features[,strandCol])
    } else {
        strand <- rep("+", nrow(features))
        ## if start/end are available, infer from from start>end
        if ( sum(c("start","end")%in%colnames(features))==2 )
          strand[features[,"start"]>features[,"end"]] <- "-"
    }
    ## re-order start>end; only for non-circular chromosomes
    ## TODO: add circular info to chrS
    if ( sum(c("start","end")%in%colnames(features))==2 & !circular ) {
        rev <- features[,"start"]>features[,"end"]
        ends <- features[rev,"start"]
        features[rev,"start"] <- features[rev,"end"]
        features[rev,"end"] <- ends
    }
    ## chromosome of each feature
    chr <- features[,chrCol]
    ## map chromosomes to index
    ## TODO: automate, if chromosomes are not numeric!?
    if ( !missing(chrMap) ) {
        chrIdx <- 1:length(chrMap)
        names(chrIdx) <- chrMap
        chr <- chrIdx[as.character(chr)]
    }

    if ( any(!is.numeric(chr)) )
        stop("chromosomes must be a numeric index; use chromosome name map with argument `chrMap'!")

    ## check for missing chromosome info and issue warning
    ## remember and also set chr to NA below
    ## TODO: check for missing info in other functions as well
    nachr <- numeric()
    if ( any(is.na(chr)) ) {
        nachr <- which(is.na(chr) )
        warning("some chromosomes are not available (NA)")
    }
    
    ## convert to index
    for ( col in cols ) {
        features[,col] <- features[,col]+chrS[chr]
        minus <- strand%in%reverse
        features[minus,col] <- features[minus,col]+max(chrS)
    }
    features[,chrCol] <- 1
    if ( length(nachr)>0 )
        features[nachr,chrCol] <- NA
    ## TODO: map so that start < end??
    
    features 
}

#' Simple version of \code{\link{index2coor}} for single values
#' @param pos the continuous index position that will be mapped to
#' chromosome coordinates
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @param strands forward/reverse strand indicators
#' @export
idx2coor <- function(pos, chrS, strands=c(1,-1)) {
  coor <- cbind(chr=rep(1,length(pos)),coor=pos,strand=rep(NA,length(pos)))
  for ( i in 1:(length(chrS)-1) ) {
    ## frw strand
    current <- pos>chrS[i] & pos<=chrS[i+1]
    coor[current,"coor"] <- pos[current] - chrS[i]
    coor[current,"chr"] <- i
    coor[current,"strand"] <- strands[1]
    ## rev strand
    current <- pos>(chrS[i]+max(chrS)) & pos<=(chrS[i+1]+max(chrS))
    coor[current] <- pos[current] - chrS[i] - max(chrS)
    coor[current,"chr"] <- i
    coor[current,"strand"] <- strands[2]
  }
  coor  
}
#' get the chromosome from continuous index
#' @param idx index position for which chromosome information is reported
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @return returns the chromosome number
#' @export
idx2chr <- function(idx,chrS) {
    chr <- sapply(idx,function(x) which(chrS>=x)[1]-1)
    if ( any(is.na(chr)) )
        chr[is.na(chr)] <- sapply(idx[is.na(chr)],function(x) # reverse strand
            which((chrS+max(chrS))>=x)[1]-1)
    chr
}
#' get the strand from continuous index
#' @param idx index position for which strand information is reported
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @return returns the strand
#' @export
idx2str <- function(idx,chrS)
    ifelse(idx > max(chrS),-1,1)



#' convert continuous index to chromosome coordinates (reverse of
#' \code{\link{coor2index}})
#' @param features a table of chromosome features that must contain
#' the chromosome number (option \code{chrCol}), one or more chromosome
#' positions (option \code{cols}) and strand information (column
#'  \code{strandCol}).
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @param chrMap a vector of chromosome names, in the same order as
#' \code{chrS}; if provided chromosome index will be mapped back to
#' chromosome name
#' @param cols names of the columns giving coordinates that will be mapped
#' to continuous index
#' @param chrCol name of the column that gives the chromosome number
#' @param relCol relative position mapping left/right -> upstream/downstream,
#' depending on strand
#' @param strandCol name of the column that gives forward/reverse strand
#' information
#' @param strands forward/reverse strand indicators
#' @export
index2coor <- function(features, chrS, chrMap,
                       cols=c("start","end","coor"),
                       chrCol="chr", strandCol="strand", relCol,
                       strands=c(1,-1)) {

    cols <- cols[cols%in%colnames(features)]

    ## add relative position column:
    ## left -> upstream/downstream, right -> downstream/upstream
    cpcols <- cols
    rel2factor <- FALSE # stores wether a relative position column was factor
    if ( !missing(relCol) ) {
        if ( relCol%in%colnames(features) )  {
            cpcols <- c(cpcols, relCol)
            ## CONVERT TO CHARACTER
            if ( class(features[,relCol])=="factor" ) {
                features[,relCol] <- as.character(features[,relCol])
                rel2factor <- TRUE                
            }
        } else
            warning("relative position column 'relCol' passed as, ",relCol,
                    "but not present in columns.")
    }
    orig <- features[,cpcols,drop=FALSE]
    
    ## add chromosome and strand columns, if not present
    if ( !chrCol%in%colnames(features) )
        features <- cbind(chr=rep(NA,nrow(features)),features)
    if ( !strandCol%in%colnames(features) )
        features <- cbind(features,strand=rep(NA,nrow(features)))
    
    ## remap values back to original coordinates
    for ( i in 1:(length(chrS)-1) ) {
        ## forward strand
        current <- orig[,cols[1]]>chrS[i] & orig[,cols[1]]<=chrS[i+1]
        for ( col in cols )
            features[current,col] <- orig[current,col] - chrS[i]
        features[current,chrCol] <- i
        features[current,strandCol] <- strands[1]
        ## relative position mapping left/right -> upstream/downstream
        if ( !missing(relCol) ) {
            tmpcol <- orig[current,relCol]
            tmpcol <- gsub("left","upstream", tmpcol)
            tmpcol <- gsub("right","downstream", tmpcol)
            features[current,relCol] <- tmpcol
        }

        ## reverse strand
        current <- orig[,cols[1]]>(chrS[i]+max(chrS)) &
            orig[,cols[1]]<=(chrS[i+1]+max(chrS))
        for ( col in cols )
            features[current,col] <- orig[current,col] - chrS[i] - max(chrS)
        features[current,chrCol] <- i
        features[current,strandCol] <- strands[2]
        ## relative position mapping left/right -> downstream/upstream
        if ( !missing(relCol) ) {
            tmpcol <- orig[current,relCol]
            tmpcol <- gsub("left","downstream", tmpcol)
            tmpcol <- gsub("right","upstream", tmpcol)
            features[current,relCol] <- tmpcol
        }
    }
    ## positions as factor
    if ( rel2factor)
        features[,relCol] <- factor(features[,relCol])

    ## return to chromosome names
    if ( !missing(chrMap) ) {
        chrIdx <- 1:length(chrMap)
        names(chrIdx) <- chrMap
        features[,chrCol] <- chrMap[features[,chrCol]]
    }
    features
}


#' switches the strand information (reverse<->forward) of genomic
#' features with continuously indexed chromosome coordinates
#' @param features genomic features with continuously indexed
#' chromosome coordinates
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @param cols names of the columns holding the continuous index
#' @export
switchStrand <- function(features,chrS, cols=c("start","end","coor")) {
  cols <- cols[cols%in%colnames(features)]
  orig <- features[,cols]
  for ( col in cols ) {
    ## forward -> reverse
    current <- orig[,col] <= max(chrS)
    features[current,col] <- orig[current,col] + max(chrS)
    ## reverse -> forward
    current <- orig[,col] > max(chrS)
    features[current,col] <- orig[current,col] - max(chrS)
  }
  features
}


#' align genome data at specified coordinates (e.g. TSS)
#' @param coors genome positions (chromosome, coordinate, strand)
#' at which data will be aligned
#' (TODO: allow start/end coors and set NA if beyond)
#' @param data genome data to be aligned; NOTE, that currently this
#' is required to be fully expanded matrix covering each chromosome position,
#' i.e. \code{nrow(data)==max(chrS)}
#' (TODO: allow non-expanded data)
#' @param dst upstream/downstream length to be aligned
## (TODO: allow different upstream and downstream ranges)
## (TODO: allow individual ranges)
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in \code{data}, derived from chromosome length
#' information, see function \code{\link{getChrSum}}
#' @param coorCols ordered string vector providing the column names
#' of coordinate columns to be used; must be of length 3 and provide in
#' order: chromosome number (refering to argument \code{chrS}), position 
#' and strand (see argument \code{reverse})
#' @param reverse a vector of possible reverse strand indicators, all other
#' values in the strand column will be taken as forward strand!
## TODO: generalize for not fully expanded data w/o chrS
## TODO: allow different downstream and upstream ranges
#' @export
alignData <- function(coors, data, dst=500, chrS,
                      coorCols=c(chr="chr", position="coor", strand="strand"),
                      reverse=c("-",-1)) {

    ## get coordinates
    starts <- as.numeric(coors[,coorCols["position"]])
    chrs <- as.numeric(coors[,coorCols["chr"]])
    strands <- as.character(coors[,coorCols["strand"]])

    ## catch wrong coordinates, eg. due to use of coor2index(coors)
    if ( any(starts>max(chrS) ) ) 
        stop(sum(starts>max(chrS)),
             " start coordinates are beyond chromosome length in index `chrS`")
    
    ## catch wrong data dimension
    if ( nrow(data)!=max(chrS) )
        stop("`data` rows (", nrow(data),
             ") do not cover the full chromosome length in index `chrS`")
    
    ## add chromosome lengths to get the index (row) in data
    starts <- chrS[chrs] + starts
    
    ## ranges in full data
    rng <- t(apply(t(starts), 2, function(x) (x-dst):(x+dst)))
    
    ##  reverse for reverse strand
    rng <- cbind(!strands%in%reverse,rng)
    rng <- t(apply(rng, 1, function(x) {

        if (x[1]==1) return(x[2:length(x)])
        else return(rev(x[2:length(x)]))

    }))

    ## cut chromosome ends
    ## TODO: implement circular chromsomes!
    ## TODO: add warning, eg. if icoors were passed negative
    ## strand values are beyond chromosome ends
    rng <- cbind(min=chrS[chrs],max=chrS[chrs+1],rng)
    rng <- t(apply(rng, 1, function(x) {

        rm <- x[3:length(x)] <= x[1] | x[3:length(x)] >= x[2];
        x <- x[3:length(x)]; x[rm] <- NA;return(x)

    }))
  
    ## split off coordinate columns, if not separately supplied
    firstcol <- 1
    if ( sum(c("chr","coor")%in%colnames(data))==2 )
        firstcol <- 3
    
    ## get aligned data for each data column
    geneData <- list()
    for ( i in firstcol:ncol(data) ) 
        geneData <- append(geneData,
                           list(t(apply(rng, 1, function(x) data[x,i]))))
    names(geneData) <- colnames(data)[firstcol:ncol(data)]
    
    ## relative coordinates as colnames
    xax <- -dst:dst
    geneData <- lapply(geneData, function(x) {colnames(x) <- xax; x})
    
    ## copy rownames
    geneData <- lapply(geneData, function(x) {rownames(x) <- rownames(coors);x})
    
    return(geneData)
}


## TODO: alignment on relative segment length
alignData_relative <- function(coors, data, dst=500, chrS,
                               coorCols=c(chr="chr", start="start", end="end",
                                          strand="strand"),
                               reverse=c("-",-1)) {

    ## TODO: sort start<ends?
    
    starts <- as.numeric(coors[,coorCols["start"]])
    ends <- as.numeric(coors[,coorCols["end"]])
    chrs <- as.numeric(coors[,coorCols["chr"]])
    strands <- as.character(coors[,coorCols["strand"]])

    ## get length and relative coordinates for each segment -1 to 2 rel.length
    ## then summarize in bins 
    
    ## add chromosome lengths to get direct index 
    starts <- chrS[chrs] + starts
    ## TODO: should rev.strand be shifted by one?
    rng <- t(apply(t(starts), 2, function(x) (x-dst):(x+dst)))
    rng <- cbind(!strands%in%reverse,rng)
    ##  reverse for reverse strand
    rng <- t(apply(rng,1,function(x)
        if (x[1]==1) return(x[2:length(x)])
        else return(rev(x[2:length(x)]))))
    
    ## cut chromosome ends
    ## TODO: implement circular chromsomes!
    ##chr <- feats[,"chr"] ## THIS NOT PASSED!?
    rng <- cbind(min=chrS[chrs],max=chrS[chrs+1],rng)
    rng <- t(apply(rng,1,function(x) {
        rm <- x[3:length(x)] <= x[1] | x[3:length(x)] >= x[2];
        x <- x[3:length(x)]; x[rm] <- NA;return(x)}))
    
    ## get data!
    ## split off coordinate columns, if not separately supplied
    firstcol <- 1
    if ( sum(c("chr","coor")%in%colnames(data))==2 )
        firstcol <- 3
    
    geneData <- list()
    for ( i in firstcol:ncol(data) ) 
        geneData <- append(geneData,
                           list(t(apply(rng, 1, function(x) data[x,i]))))
    
    names(geneData) <- colnames(data)[firstcol:ncol(data)]
    
    ## relative coordinates as colnames
    xax <- -dst:dst
    geneData <- lapply(geneData, function(x) {colnames(x) <- xax; x})
    
    ## copy rownames
    geneData <- lapply(geneData, function(x) {rownames(x) <- rownames(coors);x})
    
    return(geneData)
}

#' export internal coordinate format to bed file format
#' @param coor genome coordinates
#' @param file optional name for bed output file
#' @param coors column names of coordinate columns
#' @param reverse reverse strand characters
#' @param name column name for bed file name column (column 4)
#' @param score column name for bed file score column (column 5)
#' @param prefix prefix to be added to the name and score columns (column 4,5),
#' used by segmenTools interface to bedtools for unique names.
#' @param verb verbosity level, 0: silent
#' @seealso bed2coor
#' @export
coor2bed <- function(coor, file,
                     coors=c(chr="chr", start="start", end="end",
                             strand="strand"),
                     reverse=c("-",-1), name, score, prefix, verb=1) {

    ## add missing columns
    if ( missing(name) )  {

        rm <- which(colnames(coor)=="name") # rm existing name column
        if ( length(rm)>0 ) coor <- coor[,-rm]
        
        coor <- cbind(coor, name=paste0("id",1:nrow(coor)))
        name <- "name"
    }
    if ( missing(score) ) {
        rm <- which(colnames(coor)=="score") # rm existing score column
        if ( length(rm)>0 ) coor <- coor[,-rm]
        
        coor <- cbind(coor, score=rep(0, nrow(coor)))
        score <- "score"
    }


    ## sort start/end
    starts <- apply(coor[,coors[c("start","end")]], 1, min)
    ends <- apply(coor[,coors[c("start","end")]], 1, max)

    coor[,coors["start"]] <- starts
    coor[,coors["end"]] <- ends

    ## convert strand column
    reverse <- coor[,"strand"] %in% reverse
    coor[ reverse,coors["strand"]] <- "-"
    coor[!reverse,coors["strand"]] <- "+"

    ## order by starts
    coor <- coor[order(coor[,coors["start"]]), ]
    coor <- coor[order(coor[,coors["chr"]]), ]

    ## convert starts to 0-based
    ## NOTE: end is non-inclusive in bed-format, and
    ## thus not corrected!
    coor[,"start"] <- coor[,"start"]-1

    ## add prefix
    if ( !missing(prefix) ) {
        coor[,name] <- paste0(prefix, coor[,name])
        coor[,score] <- paste0(prefix, coor[,score])
    }

    ## chromosomes must begin with chr
    coor[,coors["chr"]] <- paste0("chr", sprintf("%02d",coor[,coors["chr"]]))
    
    bed <- coor[,c(coors[c("chr","start","end")], name, score,coors["strand"])]

    if( !missing(file) ) {
        if ( verb>0 )
            cat(paste("writing bed file:", file,"\n"))
        write.table(x=bed, file=file, sep="\t",
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    invisible(bed)
}

