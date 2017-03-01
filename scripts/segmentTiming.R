
library(segmenTools) # for segment classes

## computation timing file was constructed by grepping entries in the
## runSegmentier <log files>, as recorded by stdout redirection:
## grep "SEGMENT TYPE\|elapsed\|Timeseries" <logfiles> | sed 's/.*SEGM/SEGM/g;s/.*elapsed/elapsed/g;s/.*Timeseries/Timeseries/' > timing.csv
timing.file <- "timing.csv"
dat <- read.table(timing.file, sep="\t", header=FALSE,fill=TRUE)


## get pre-segments
tsidx <- which(dat[,1]=="Timeseries N")
nn <- length(tsidx)
tsidx <- c(tsidx,nrow(dat)+1)

## get only values where the range between two Timeseries entries
## indicates proper run: each "Timeseries" must be followed by
## pairs of "SEGMENT TYPE" and "elapsed, sec" entries
lst <- lapply(1:nn, function(x) {
    from <- tsidx[x]+1
    to <- tsidx[x+1]-1
    sq <- NULL
    if ( from<=to & (to-from+1)%%2==0)
        sq <- seq(from,to,by=1)
    sq})

wrong <- which(unlist(lapply(lst, is.null)))
if ( length(wrong)>0 )
    lst <- lst[-wrong]

## check if the entries in column 1 are as expected
## takes LONG
check <- FALSE
if ( check ) {
    chck <- lapply(lst, function(x) {
        data.frame(length=rep(dat[x[1]-1,1],length(x)/2),
                   type=dat[x[seq(1,length(x)-1,2)],1],
                   time=dat[x[seq(2,length(x),2)],1])
    })
    chck <- do.call(rbind, chck)
    apply(chck,2,unique)
}

sglens <- lapply(lst, function(x) {
    data.frame(length=rep(as.numeric(as.character(dat[x[1]-1,2])),length(x)/2),
               type=           dat[x[seq(1,length(x)-1,2)],2],
               time=as.numeric(as.character(dat[x[seq(2,length(x),2)],2])))
})

## mean time per segment length
mn.time <- unlist(lapply(sglens, function(x) mean(x[,3])))
mn.len <- unlist(lapply(sglens, function(x) mean(x[,1])))
## all segments table!
sglens <- do.call(rbind, sglens)

## generate class table and analyze by class
types <- as.character(unique(sglens[,2]))
color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

typ.col <- color_hue(length(types))
names(typ.col) <- types


xlim <- c(900,max(sglens[,1]))
ylim <- c(min(sglens[sglens[,3]>0,3]),max(sglens[,3]))
ylim <- c(1/60,2e2)

## segmentation classes
cltab <- getSegmentClassTable(types)

for ( cl in colnames(cltab) ) {
    cat(paste("plotting by", cl, "\n"))
    clk <- unique(cltab[,cl])
    clk.cols <- color_hue(length(clk))
    names(clk.cols) <- clk
    png(paste("calculation_timing_by",cl,".png",sep=""),
        res=200,width=5,height=3,units="in")
    par(mai=c(.6,.6,.1,.1),mgp=c(1.5,.5,0))
    plot(mn.len,mn.time,col=NA,xlab="segment length, bp",ylab="calculation time, min",log="xy",xlim=xlim,ylim=ylim)
    for ( k in clk ) {
        t <- rownames(cltab)[cltab[,cl]==k]
        ft <- sglens[,2]%in%t
        points(sglens[ft,1],sglens[ft,3]/60,col=paste(clk.cols[as.character(k)],"11",sep=""),pch=4,cex=.3)
    }
    legend("topleft",names(clk.cols),col=clk.cols,pch=4)
    dev.off()
}

## segment classes
png("calculation_timing.png",res=200,width=5,height=3,units="in")
par(mai=c(.6,.6,.1,.1),mgp=c(1.5,.5,0))
plot(mn.len,mn.time,col=NA,xlab="segment length, bp",ylab="calculation time, min",log="xy",xlim=xlim,ylim=ylim)
for ( t in types ) {
    ft <- sglens[,2]==t
    points(sglens[ft,1],sglens[ft,3]/60,col=paste(typ.col[t],"11",sep=""),pch=4,cex=.3)
}
dev.off()


## cumulative timing
cumtime <- cumsum(sort(sglens[,3]))
png("calculation_timing_cumulative.png",res=200,width=5,height=2.5,units="in")
par(mai=c(.6,.6,.1,.1),mgp=c(1.5,.5,0))
plot(cumtime/60/60,type="l",ylab="cumul. calc. time, hours",xlab="sorted jobs",log="y")
legend("right", paste("total:",round(max(cumtime)/60/60/24), "days"))
dev.off()
