#!/usr/bin/Rscript

## reads *.bed genome files from an input directory and
## a chromosome index file and expands the bed files to
## full coordinate files (chr\tcoor\tvalue)

verb <- FALSE

args  <- R.utils::commandArgs(excludeReserved=TRUE, asValues=TRUE)
for ( i in 2:length(args) ) {
  arg <- names(args)[i]
  ## fix(?) since R 3.0.0 ?
  if ( is.null(args[[arg]]) ) args[arg] <- TRUE
  assign(arg, as.character(args[arg]))
}
if ( !exists("idx",mode="character") ) {
  cat(paste("no chromosome index file provided, use --idx=<file>\n"))
  if ( !interactive() ) quit(save="no")
}
## input (should contain bed files) and output directories
if ( !exists("dir",mode="character") ) dir <- "."
if ( !exists("out",mode="character") ) out <- dir
if ( !exists("pat",mode="character") ) pat <- ""
verb <- as.logical(verb)

## create output directory
dir.create(out)


## read chromosome index file
## should contain for each chromosome ID\tlength
chrIdx <- read.table(idx, sep="\t",header=FALSE)

## bed files
files <- Sys.glob(file.path(dir,paste(pat,"*.bed",sep="")))


for ( i in 1:length(files) ) {
  
  if ( verb ) cat(paste("reading", files[i], "\n"))

  ## bed file data
  bed <- read.table(files[i],sep="\t")

  ## vector to hold data from bed file
  chrDat <- matrix(NA, ncol=3, nrow=sum(chrIdx[,2]))

  ## expand bed to full genome length
  oldchr <- "none"
  for ( j in 1:nrow(bed) ) {
    chr <- as.character(bed[j,1])
    if ( chr!=oldchr & verb ) {
      cat(paste("chromosome", chr, "\n"))
      oldchr <- chr
    }
    chrNum <- which(as.character(chrIdx[,1])==chr)
    chrPos <- ifelse(chrNum==1,0,cumsum(chrIdx[,2])[chrNum -1])
    start <- bed[j,2]+1 + chrPos
    end <- bed[j,3]+ chrPos
    chrDat[start:end, 3] <- bed[j,4]
   }

  ## coordinates
  for ( j in 1:nrow(chrIdx) ) {
    start <- ifelse(j==1,1,cumsum(chrIdx[,2])[j-1]+1)
    end <- cumsum(chrIdx[,2])[j]
    chrDat[start:end,1:2] <-  cbind(rep(j,chrIdx[j,2]),1:chrIdx[j,2])
  }
  file <- sub(".*/","",files[i])
  colnames(chrDat) <- c("chr","coor",file)

  ## write file
  file <- file.path(out, sub(".bed", ".csv", file))
  if ( verb ) cat(paste("writing", file, "\n"))
  write.table(chrDat,file, sep="\t", quote=FALSE, na="", row.names=FALSE)
  
  ## rm to save memory (required?)
  rm(chrDat)
  rm(bed)
}
  
# End of file
