
library(segmenTools) # for segment classes

timing.file <- "timing.csv"
data <- read.table(timing.file, sep="\t", header=FALSE,fill=TRUE)

dat <- data
tsidx <- which(dat[,1]=="Timeseries N")

nn <- length(tsidx)
tsidx <- c(tsidx,nrow(dat)+1)
## TODO: whats the error, why is "try" required?
lst <- lapply(1:nn, function(x) {
    from <- tsidx[x]+1
    to <- tsidx[x+1]-1
    sq <- NULL
    if ( from<=to & (to-from+1)%%2==0)
        sq <- seq(from,to,by=1)
    sq})
wrong <- which(unlist(lapply(lst, is.null)))

tsidx <- tsidx[-wrong]
lst <- lst[-wrong]

chck <- lapply(lst, function(x) {
    data.frame(length=rep(dat[x[1]-1,1],length(x)/2),
               type=dat[x[seq(1,length(x)-1,2)],1],
               time=dat[x[seq(2,length(x),2)],1])
})
chck <- do.call(rbind, chck)
apply(chck,2,unique)

sglens <- lapply(lst, function(x) {
    data.frame(length=rep(as.numeric(dat[x[1]-1,2]),length(x)/2),
               type=dat[x[seq(1,length(x)-1,2)],2],
               time=as.numeric(dat[x[seq(2,length(x),2)],2]))
})
sglens <- do.call(rbind, sglens)

types <- as.character(unique(sglens[,2]))
color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

typ.col <- color_hue(length(types))
names(typ.col) <- types

png("calculation_timing.png",res=200,width=5,height=3,units="in")
par(mai=c(.6,.6,.1,.1),mgp=c(1.5,.5,0))
plot(1,xlab="segment length, bp",ylab="calculation time, min",log="",col=sglens[,2],xlim=range(sglens[,1]),ylim=range(sglens[,3]/60))
for ( t in types ) {
    ft <- sglens[,2]==t
    points(sglens[ft,1],sglens[ft,3]/60,col=typ.col[t])
}
dev.off()
