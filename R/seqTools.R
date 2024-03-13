
### DNA/RNA/AA SEQUENCE UTILS
### 20180423: copied from genomeBrowser_utils.R

#' calculate reverse complement of RNA or DNA
#'
#' Returns the reverse complement of RNA/DNA bases in an input string.
#' If \code{reverse==FALSE} it just reports complementary bases in the
#' same order.
#' @param sq a, RNA or DNA string
#' @param type "RNA" or "DNA"
#' @param na string to use for non-recognized nucleotides
#' @param reverse revert to sequence
#'@export
revcomp <- function(sq, type="DNA", na="N", reverse=TRUE) {
    ## TODO: instead use base R chartr function
  if ( type=="DNA" ) bp <- c(A="T",T="A",C="G",G="C")
  else if ( type=="RNA" ) bp <- c(A="U",U="A",C="G",G="C")
  sq <- unlist(strsplit(sq,""))
  rv <- bp[sq]
  rv[is.na(rv)] <- na
  if ( reverse ) rv <- rev(rv)
  paste(rv,collapse="")
}

#' calculates letter frequences in ranges of strings
#'
#' calculates letter frequences in ranges of strings, e.g.,
#' nucleotide sequences. The input is a sequence object
#' as returned by \code{readFASTA}.
#' @param seq a sequence object as returned by \code{readFASTA}
#' @param ranges a table provided ranges in each string 
#' @param coors coordinate names in ranges table
#' @export
range2nt <- function(seq, ranges, coors=c(chr="chr",start="start",end="end") ) {

    ranges <- ranges[,coors]
    nt <- list(rep(NA), nrow(ranges))
    sqs <- lapply(seq, function(x) unlist(strsplit(x$seq,"")))

    ## count all letters in each 
    t(apply(ranges, 1, function(x) {
        x <- unlist(x)
        sq <- sqs[[x[coors["chr"]]]][x[coors["start"]]:x[coors["end"]]]
        table(sq)/length(sq)
    }))
 }

#' calculate local nucleotide content
#'
#' calculate local nucleotide content in moving windows
#' and returns genomeData matrix
#' TODO: move readFasta here and document!
#' @param seq a sequence object as returned by \code{readFasta}
#' @param window width of the moving window
#' @param step step size
#' @param abc sequence letters for which averages are calculated
#' @param skew calculate skew for letters 1 and 2 in argument \code{skew}
#' (e.g. \code{skew=c("G","C"), skew=(G-C)/(G+C)})
#' @param cum set TRUE to calculate cumulative skew after Grigoriev 1998 NAR
#' @param circular treat sequence as circular
#' @param verb print progress messages
#'@export
seq2nt <- function(seq, window=100, step=1, abc, skew="",
                   cum=FALSE, circular=FALSE, verb=FALSE) {

  ntprofiles <- NULL
  for ( i in 1:length(seq) ) {

    if ( verb ) cat(paste("chromosome", i,", "))

    sq <- unlist(strsplit(seq[[i]]$seq,""))
    if ( missing(abc) ) 
      abc <- sort(unique(sq))

    if ( skew!="" )
      abc <- unlist(strsplit(skew,""))

    if ( verb ) cat(paste("scanning", paste(abc,collapse=","),":"))
    
    ## convert sequence to 0/1 table
    ntl <- matrix(NA, ncol=length(abc), nrow=length(sq))
    colnames(ntl) <- abc
    for ( nuc in abc )
      ntl[,nuc] <- as.numeric(sq == nuc)

    len <- ifelse(circular,length(sq),length(sq)-window)
    steps <- seq(1,len,step)
    pos <- rep(NA, len)
    if ( skew=="" ) {
      nt <- matrix(NA, ncol=length(abc), nrow=length(steps))
      colnames(nt) <- abc
    }else{
      nt <- matrix(NA, ncol=1, nrow=length(steps))
      colnames(nt) <- skew
    }
    for ( j in 1:length(steps) ) {
      if ( verb & j%%round(length(steps)/5,0) == 0)
        cat(paste(" ", round(j/length(steps),2)))
      start <- steps[j]
      bin <- start:(start+window)
      ## wind edges
      if ( circular ) {
        bin[bin>length(sq)] <- bin[bin>length(sq)] - length(sq)
      }
      ## cut all remaining - NOTE: there shouldn't be any!
      if ( max(bin)>length(sq) )
        cat("DEVEL WARNING: bin shouldn't contain positions beyond seqlen\n")
      bin <- bin[bin<=length(sq)]

      if ( length(bin)==1 ) {
        pos[j] <- bin
      }else {
        pos[j] <- bin[floor(length(bin)/2)]
      }
      if ( skew != "" ) {
        ## calculate skew (.e.g. (G-C)/(G+C)
        nt1 <- sum(ntl[bin,abc[1]])
        nt2 <- sum(ntl[bin,abc[2]])
        nt[j,] <- (nt1 - nt2)/(nt1 + nt2)
      }else{
        ## calculate frequency of each nucleotide in window
        for ( nuc in abc )
          nt[j,nuc] <- mean(ntl[bin,nuc])
      }
    }
    if ( verb ) cat("\n")
    ## cumulative GC-skew, norm. for window and seq. lengths
    ## see Grigoriev 1998 NAR
    if ( cum ) nt[,1] <- cumsum(nt[,1])*window/length(sq)
    ntprofiles <- rbind(ntprofiles,
                        cbind(chr=rep(i, length(pos)),coor=pos,nt)[order(pos),])
  }
  ntprofiles
}

#' find genomic coordinates of amino acid position in a protein on genome
#' TODO: this does not account for introns
#'
#' TODO: to solve this generally, the amino acid sequence of the protein
#' must be blasted against the genome
#' @param start start position of the ATG start codon
#' @param strand coding strand 
#' @param aa position of the amino acid in the protein
#'@export
findAACodon <- function(start, aa, strand) {
    offset <- (aa-1)*3 ## 3 positions per amino acide; start at one lower
    if ( as.character(strand)%in%c("1","+","+1") )
      pos <- (start-1) + offset
    else if ( as.character(strand)%in%c("-1","-") )
      pos <- (start+1) - offset
    pos
}

#' mutate positions in a string
#' @param seq a string to be mutated
#' @param mutations a list of strings of the format 'A3T'
#' where the A at position 3 of the string should be changed
#' to a T
#' @export
mutatePositions <- function(seq, mutations) {

    str <- strsplit(seq,"")[[1]]
    fromto <- do.call(cbind,strsplit(mutations, "[0-9]+"))
    at <- as.numeric(gsub("\\*","",gsub("[A-Z]","", mutations)))
    if ( any(duplicated(at)) )
        stop("duplicated mutation loci are not supported")
    if ( any(is.na(at)) )
        stop("location not found")
    for ( j in 1:ncol(fromto) ) {
        if ( str[at[j]]!=fromto[1,j] )
            stop("mutated position ", j, " is not as indicated.",
                 "searched: ",fromto[1,j], ", but found: ", str[at[j]])
        else str[at[j]] <- fromto[2,j]
    }
    paste(str, collapse="")
}

## TODO:
## di/tri-nucleotide parameter profile
## position weight matrix 
