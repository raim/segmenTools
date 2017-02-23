
timing.file <- "timing.csv"
data <- read.csv(timing.file, sep="\t", header=FALSE)

## by segment length
dat <- data[data[,1]!="SEGMENT TYPE",]
tsidx <- which(dat[,1]=="Timeseries N")
sglen <- dat[tsidx,2]
sglen <- as.numeric(levels(sglen))[sglen]

times <- as.numeric(levels(dat[,2]))[dat[,2]]

nn <- length(tsidx)
tsidx <- c(tsidx,nrow(dat))
## TODO: whats the error, why is "try" required?
lst <- lapply(1:nn, function(x) {
    from <- tsidx[x]+1
    to <- tsidx[x+1]-1
    sq <- NA
    if ( from<=to )
        sq <- seq(from,to,by=1)
    sq})

times <- lapply(lst, function(x) times[x])
mn.times <- unlist(lapply(times, mean))
plot(sglen,mn.times/60,xlab="segment length, bp",ylab="calculation time, min")

## TODO: by class!
