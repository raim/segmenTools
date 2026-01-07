library(segmenTools)

## DOCUMENT LOG/ARCSINH FUNCTIONALITY

x <- -1000:1000

plotdev('arcsinh_log_comparison', width = 3.5, height = 3.5)
par(mai=rep(.5, 4), mgp = c(1.3, .3, 0), tcl = -.25)
dense2d(x, log10(x), axes=FALSE, ylim=c(-10,10), ylab=NA)
axis(1)
axis(2)
abline(h=-10:10, col=8, lwd=.1)
mtext(expression(log[10](x)), 2, 1.3)

par(new=TRUE)
dense2d(x, ash(x), axes=FALSE, ylim=c(-10,10), xlab=NA, ylab=NA,
        colf=grey.colors)
axis(4)
mtext(expression(arcsinh(x)), 4, 1.3)
legend('topleft', expression(arcsinh(x), log[10](x)),
       col=c(viridis::viridis(5)[4], 'gray'), pch=19, inset = c(.05, -.15),
       xpd = TRUE)
dev.off()


plotdev('arcsinh_log_logaxis', width = 3.5, height = 3.5)
par(mai=rep(.5, 4), mgp = c(1.3, .3, 0), tcl = -.25)
dense2d(x, log10(x), axes=FALSE,  ylab=NA)
axis(1)
logaxis(2)
#abline(h=-10:10, col=8, lwd=.1)
mtext(expression(log[10](x)), 2, 1.3)
dev.off()

plotdev('arcsinh_log_ashaxis', width = 3.5, height = 3.5)
par(mai=rep(.5, 4), mgp = c(1.3, .3, 0), tcl = -.25)
dense2d(x, ash(x), axes=FALSE, ylab=NA,
        colf=grey.colors)
axis(1)
ashaxis(2)
mtext(expression(arcsinh(x)), 2, 1.3)
dev.off()

plotdev('arcsinh_log_lplot', width = 3.5, height = 3.5)
par(mai=rep(.5, 4), mgp = c(1.3, .3, 0), tcl = -.25)
lplot(x, x, log = 'x', ash = 'y', xlab = 'log10-scaled x',
      ylab = 'arscinh-scaled x', ylim=c(0, ash(max(x))))
logaxis(1)
ashaxis(2)
dev.off()

plotdev('arcsinh_log_plot', width = 3.5, height = 3.5)
par(mai=rep(.5, 4), mgp = c(1.3, .3, 0), tcl = -.25)
lplot(log10(x), ash(x), xlab = 'log10-scaled x',
      ylab = 'arscinh-scaled x', ylim=c(0, ash(max(x))))
axis(1)
axis(1, at=seq(0,10,.5), labels=FALSE)
axis(2)
ashaxis(4)
logaxis(3)
dev.off()

