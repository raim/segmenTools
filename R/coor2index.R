## GENOME UTILS

#' generate chromosome index \code{chrS} from and ordered (!)
#' chromosome length table
#' @param chrIdx an ordered chromosome length table; feature
#' files for these chromosomes must have a column "chr" which
#' used the row number in chrIdx as its chromosome identifier
#' @param lcol numeric or character indicating the  number
#' or name of the column with chromosome lengths
#' @export
getChrS <- function(chrIdx, lcol="length")
    c(0,cumsum(chrIdx[,lcol]))

## TODO: for overlap analysis and browser, detect
## overlaps in circular chromosomes; generate two genes for each
## end, re-tag original as gene_parent_circular
expandCircularFeatures <- function(features, chrS, chrC) {
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

