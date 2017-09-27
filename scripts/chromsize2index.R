#!/usr/bin/Rscript
# chromsize2index.R

## converts a simple chromosome size file to the sequence index file
## used in tataProject; the chromosome size file must have the
## form id\tlength, and should give chromosomes in the correct order
## (except for yeast where the chrI convention with roman literals
## can be translated to the correct order)

## use, e.g., in tataProject/yeast:
# $TATADIR/r/chromsize2index.R --inf $TATADIR/yeast/originalData/Yeast_S288C.chromsizes --out $TATADIR/yeast/chromosomes/sequenceIndex_R64-1-1_20110208.csv --yeast

yeast <- FALSE

args  <- R.utils::commandArgs(excludeReserved=TRUE, asValues=TRUE)
for ( i in 2:length(args) ) {
  arg <- names(args)[i]
  ## fix(?) since R 3.0.0 ?
  if ( is.null(args[[arg]]) ) args[arg] <- TRUE
  assign(arg, as.character(args[arg]))
}

if ( !exists("inf",mode="character") ) {
  cat("no size chromsome size file provided, use --inf=<file>\n",file=stderr())
  if ( !interactive() ) quit(save="no")
}
if ( !exists("out",mode="character") ) out <- ""

yeast <- as.logical(yeast)


chr <- read.table(inf,sep="\t")
num <- 1:nrow(chr)
ids <- as.character(chr[,1])

## catch yeast genomes
if ( yeast ) {
  
  ## get roman literals
  num <- as.numeric(as.roman(sub("chr","",chr[,1])))

  ## order
  chr <- chr[order(num),]
  num <- sort(num)
  ## replace mitochondrial genome index
  num[17] <- 17
  ids[17] <- "chrMito"
} 

idx <- cbind(num, ids, chr[,2])
colnames(idx) <- c("#ID","name","length")

if ( out=="" ) {
  write.table(idx, file=stdout(), quote=FALSE, row.names=FALSE, sep="\t") 
} else {
  write.table(idx, file=out, quote=FALSE, row.names=FALSE, sep="\t")
}
  
                                        
