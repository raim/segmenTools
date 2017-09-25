## GENOME UTILS

#' Generate chromosome index \code{chrS} from lengths
#' @param chrL an ordered vector of chromosome lengths; where the
#' order must correspond to chromosome numbering in feature tables
#' for which chrS is used
#' @export
getChrS <- function(chrL) c(0,cumsum(chrL))

## util to insert rows, by user Ari B. Friedman at
## https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended
insertRow <- function(existingDF, newrow, r) {
    existingDF <- as.data.frame(existingDF,stringsAsFactors=FALSE)
    existingDF[seq(r+1,nrow(existingDF)+1),] <-
        existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
}
#' Util to insert multiple rows at specified positions
#' in a \code{data.frame}, expanding single-row code by user
#' Ari B. Friedman at
#' \url{https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appendedlooping through new rows}
#' @param existingDF existing \code{data.frame}
#' @param newrows rows to add to \code{existingDF}
#' @param r positions at which rows are to be inserted, \code{length(r)}
#' must equal \code{nrow(newrows)}, and \code{r<=nrow(existingDF)}
#' @export
insertRows <- function(existingDF, newrows, r ) {
    new <- existingDF
    for ( i in 1:nrow(newrows) )
        new <- insertRow(new, newrows[i,], r[i]+i-1)
    new
}

#' Splits genome features crossing ends of circular chromosomes.
#' Splits genome features that cover start/end coordinates of circular
#' chromosomes, and adds the downstream half with optional modification
#' of ID, and type values. Note that only the upstream half retains
#' all column information (exceptions: see argument \code{copyCols}),
#' the downstream half will only carry information on coordinates, and
#' optionally updated feature type and ID.
#' The update will only happen if the passed table contains type and ID
#' information (see argument \code{idCools}. The split can be reversed
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
#' @seealso \code{removeCircularFeatures}
#' @export
expandCircularFeatures <- function(features, chrL, 
                                   coorCols=c("chr","start","end","strand"),
                                   reverse=c("-",-1),
                                   idTag="-circ2", idCols=c(ID="ID",type="type",
                                          parent="parent"),copyCols=FALSE) {

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
    ## set up type
    if ( idCols["ID"]%in%colnames(features) ) {
        cfeat[,idCols["parent"]] <- features[cidx,idCols["ID"]]
        cfeat[,idCols["ID"]] <- paste(features[cidx,idCols["ID"]],idTag,sep="")
    }
    if ( idCols["type"]%in%colnames(features) )
      cfeat[,idCols["type"]] <- paste(features[cidx,idCols["type"]],
                                      idTag,sep="")
    ## copy requested columns
    cfeat[,copyCols] <- features[cidx,copyCols]
    
      
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
    insertRows(features,cfeat,cidx+1)
}

#' NOT WORKING - Undo \code{expandCircularFeatures}
#' searches for circular features by the \code{idTag} added
#' to ID and type columns in \code{expandCircularFeatures},
#' and maps downstream coordinates back to original features.
#' @param features  list of genomic features with coordinates
#' @param idTag tag used for tagging downstream halves
#' @param idCols named vector of column names for feature ID, type,
#' and feature parent; note that a "parent" column will be removed if
#' it is (a) empty and (b) argument \code{rmParent==TRUE}
#' @param rmParent the column parent
#' @seealso \code{expandCircularFeatures}
#' @export
removeCircularFeatures <- function(features, idTag="-circ2", idCols=c(ID="ID",type="type",parent="parent"), rmParent=TRUE) {
}

#' convert chromosome coordinates to continuous index
#' @param features a table of chromosome features that must contain
#' the chromosome number (option \code{chrCol}), one or more chromosome
#' positions (option \code{cols}) and strand information (column
#'  \code{strandCol}).
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information, see function \code{\link{getChrS}}
#' @param cols name of the columns giving coordinates that will be mapped
#' to continuous index
#' @param chrCol name of the column that gives the chromosome number
#' @param strandCol name of the column that gives forward/reverse strand
#' information
#' @param reverse a vector of possible reverse strand indicators
#' @param circular suppresses re-sorting to start < end for circular chromosomes
#' @export
coor2index <- function(features, chrS,
                       cols=c("start","end","coor"),
                       chrCol="chr", strandCol="strand",
                       reverse=c("-",-1), circular=FALSE) {

    ## coordinate columns
    cols <- cols[cols%in%colnames(features)]
    ## strand column - if not present, infer from start>end
    if ( strandCol%in%colnames(features) ) {
        strand <- features[,strandCol]
    } else {
        strand <- rep("+", nrow(features))
        strand[features[,"start"]>features[,"end"]] <- "-"
    }
    ## re-order start>end; only for non-circular chromosomes
    ## TODO: add circular info to chrS
    if ( sum(c("start","end")%in%colnames(features))==2 &
         !circular ) {
        rev <- features[,"start"]>features[,"end"]
        ends <- features[rev,"start"]
        features[rev,"start"] <- features[rev,"end"]
        features[rev,"end"] <- ends
    }
    ## convert to index
    for ( col in cols ) {
        features[,col] <- features[,col]+chrS[features[,chrCol]]
        minus <- strand%in%reverse
        features[minus,col] <- features[minus,col]+max(chrS)
    }
    features[,chrCol] <- 1
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
#' positions (option \code{coorCols}) and strand information (column
#'  \code{strandCol}).
#' @param chrS the chromosome index, indicating the start position
#' of each chromosome in the continuous index, derived from chromosome length
#' information
#' @param cols names of the columns giving coordinates that will be mapped
#' to continuous index
#' @param chrCol name of the column that gives the chromosome number
#' @param relCol relative position mapping left/right -> upstream/downstream,
#' depending on strand
#' @param strandCol name of the column that gives forward/reverse strand
#' information
#' @param strands forward/reverse strand indicators
#' @export
index2coor <- function(features, chrS,
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
    if ( rel2factor)
         features[,relCol] <- factor(features[,relCol])
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

