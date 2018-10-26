#!/usr/bin/env Rscript

## TODO - bugs
## NULL at begnning when just redirecting to stdout

## This script
## processes circular features that have been calculated for 
## artificially extended circular chromosomes;
## 1) deletes features within the extension (without checking for existence!)
## 2) detects features that span start/end coordinates and corrects their end

## NOTE: the result file will indicate circular features by start > end!
## strand information must be present in strand column!

## TODO: input refers either to extended chromosome, which will
## lead to deletion and coordinate fixing, OR input has already
## fixed coordinates (ext=0), with start > end indicating start/end
## spanning features

## requires input:
## the chromosome index of the circularized chromosomes,
## and the extension (bp) that was used for circularization,


library(segmenTools)
## nicer timestamp
time <- function() format(Sys.time(), "%Y%m%d %H:%M:%S")
## messages
msg <- function(x) cat(x, file=msgfile)

debug <- FALSE

## PARSE OPTIONS
suppressPackageStartupMessages(library(optparse))
option_list <- list(
    ## GENOME OPTIONS
    make_option(c("--chrfile"), type="character", default="",
                help="chromosome index file, providing a sorted list of chromosomes and their lengths in column 3 [default %default]"),
    make_option(c("-e", "--ext"), type="integer", default=5000, 
                help="size of circularization extension"),
    ## QUERY OPTIONS
    make_option(c("-i", "--infile"), type="character", default="", 
                help="infile: set of chromosomal segments"),
    make_option(c("--id"), type="character", default="ID", 
                help="ID column name; ignored if missing!"),
    make_option(c("--only.split"), action="store_true", default=FALSE,
                help="the input has already fixed coordinates, with start > end
indicating not the strand, but features spanning the end; only option --split will be used to split features in upstream/downstream halves at the start/end"),
    make_option(c("--fix.neg"), action="store_false", default=TRUE,
                help="subtract chromosome length from features with negative start coordinates"),
    ## OUTPUT OPTIONS
    make_option(c("--split"), action="store_true", default=FALSE,
                help="splits features spanning start/end in two and adds a postfix `-circ2' and a parent column entry to the downstream end "),
    make_option(c("--sort.rows"), action="store_true", default=FALSE,
                help="3' half of split circular children will be inserted below their 5' parent; can be slow, use --sort.rows to activate"),
    make_option(c("-o", "--outfile"), type="character", default="", 
                help="file name to write circularized version of infile"),
    make_option(c("-v", "--verb"), type="integer", default=1, 
                help="0: silent, 1: main messages, 2: warning messages"))

## get command line options
opt <- parse_args(OptionParser(option_list=option_list))

## TODO: process comma-separated list arguments

## promote options to main environment 
for ( i in 1:length(opt) ) {
    arg <- names(opt)[i]
    assign(arg, opt[[arg]])
}
## select file pointer for output and messages
if ( outfile=="" ) {
    outfile <- stdout()
    msgfile <- stderr()
} else {
    msgfile <- stdout()
}

if ( interactive() & debug ) {
    chrfile <- "~/data/introns/prokaryotes_2018/allcontigs/allcontigs.idx"
    infile <- "~/data/introns/prokaryotes_2018/allcontigs/gI_hits.tab"
    ext <- 5000
    chrmap <- TRUE
}

## get chromosome index and extension
if ( verb>0 )
    msg(paste("Loading chromosome index file:", chrfile, "\t\n"))
cf <- read.table(chrfile,sep="\t",header=FALSE, stringsAsFactors=FALSE, quote="")
chrL <- as.numeric(cf[,ncol(cf)]) ## length
chrS <- c(0,cumsum(chrL)) ## index of chr/pos = chrS[chr] +
## name information in first column
chrMap <- cf[,1]
names(chrS) <- c("",chrMap)
names(chrL) <- chrMap # chromosome length hash

## real length: remove extension
chrRL <- chrL - ext

## load segments
if ( verb>0 ) msg(paste("Loading infile:", infile, "\t\n"))
input <- read.table(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

## ONLY HANDLE SEGMENTS WITH GENOME INFO!!
missing <- which(!input$chr%in%names(chrL))
if ( length(missing)>0 ) {
    ## TODO: allow storing them, but operations not requiring chrL
    msg(paste("WARNING:", length(missing),
              "features with chromosomes not listed in chromosome index file",
              "these will be lost!",
              paste(unique(input$chr[missing]),collapse=";"), "\n"))
    #missdat <- input[missing,]
    input <- input[-missing,]
}

## test numeric character of start/end
if ( !is.integer(input$start) ) {
    idx <- which(is.na(as.integer(input$start)))
    if ( length(idx) )
        stop("non-integer start columns in lines ", paste(idx,collapse=";"))
    else input$start <- as.numeric(input$start)
}
if ( !is.integer(input$end) ) {
    idx <- which(is.na(as.numeric(input$end)))
    if ( length(idx) )
        stop("non-integer end columns in lines ", paste(idx,collapse=";"))
    else input$end <- as.numeric(input$end)
}


if ( !only.split ) {
    
    ## order start and end
    rev <- input[,"start"]>input[,"end"]
    ends <- input[rev,"start"]
    input[rev,"start"] <- input[rev,"end"]
    input[rev,"end"] <- ends
    
    
    ## scan segments for coordinates longer then chromosome
    ## 1) rm segments within the extension (optional check whether
    ##    they do exist in beginning of chromosome)
    rm <- ( input[,"end"]>chrRL[input[,"chr"]] &
            input[,"start"]>chrRL[input[,"chr"]] )
    
    if ( verb>0 )
        msg(paste("removing", sum(rm), "hit(s) within extension",ext,"\n"))
    input <- input[-which(rm),]
    
    ## 2) detect hits that span start/end,
    span <- ( input[,"end"]>chrRL[input[,"chr"]] &
              input[,"start"]<=chrRL[input[,"chr"]] )

    ## correct end coordinate, start will be > end!
    input[span,"end"] <- input[span,"end"] - chrRL[input[span,"chr"]]

    ## 2a) detect partial hits that are covered by the span hits
    ## simple heuristic: same end?
    ## TODO: make this saver
    dupls <- c()
    for ( sp in which(span) )
        dupls <- c(dupls, which( input[,"chr"]==input[sp,"chr"] &
                                 input[,"end"]==input[sp,"end"] ))
    covered <- rep(FALSE, length(span))
    covered[dupls] <- TRUE
    

    if ( verb>0 )
        msg(paste("fixing", sum(span), "hit(s) spanning start/end\n"))
    if ( verb>0 )
        msg(paste("removing", sum(covered &!span),
                  "hit(s) covered by end hit\n"))
    
    result <- input[-which(covered & !span),]
} else {
    ## assumes that features are already ordered with start<end
    ## and start>end implying features spanning circular ends
    result <- input
    span <- input[,"start"]>input[,"end"]
}

## 3) detect and correct negative values
if ( fix.neg ) {
    negs <- result[,"start"] < 0
    if ( verb>0 )
        msg(paste("fixing", sum(negs), "hit(s) with negative start coordinate\n"))
    if ( sum(negs)>0 )
        result[negs,"start"] <- chrRL[result[negs,"chr"]] +  result[negs,"start"]
}

## 4) split hits
if ( split | only.split ) {
    if ( verb>0 )
        msg(paste("splitting", sum(span), "hit(s) spanning start/end\n"))
    result <- expandCircularFeatures(result, chrL=chrRL, insertRows=sort.rows, copyCols=TRUE, idCols=c(ID=id,type="type",parent="parent"))
}

## TODO: allow adding unknown chromosomes back?
#if ( length(missing)>0 ) {
#}

if ( verb>0 )
    msg(paste("DONE. WRITING RESULTS\t",time(),"\n",sep=""))
if ( verb>0 )
    msg(paste("Writing result:", outfile, "\t\n"))

write.table(result, file=outfile, sep="\t", quote=FALSE,row.names=FALSE, na="")

if ( verb>0 )
    msg(paste("writing\t", nrow(result), " segments\n",sep=""))
