library(segmenTools)
library(viridis)

## TODO: add to function as an example

x <- rnorm(1000)
png("example_selectColors.png", width=400, height=300)
colors <- segmenTools::selectColors(x, colf=viridis::viridis, mx=2, mn=-1)
dev.off()

## same colors
cols <- segmenTools::num2col(x, limits=c(-1,2), colf=viridis::viridis)
