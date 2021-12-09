### PLOT UTILS
#' Switch between plot devices
#' @param file.name file name without suffix (.png, etc)
#' @param type plot type: png, jpeg, eps, pdf, tiff or svg
#' @param width figure width in inches
#' @param height figure height in inches
#' @param res resolution in ppi (pixels per inch), only for 'png' and 'tiff'
#' @param bg background color
#' @export
plotdev <- function(file.name="test", type="png", width=5, height=5, res=100,
                    bg="white") {
  file.name <- paste(file.name, type, sep=".")
  if ( type == "png" )
    grDevices::png(file.name, width=width, height=height, units="in", res=res)
  if ( type == "eps" )
    grDevices::postscript(file.name, width=height, height=width,paper="special",
                          horizontal = FALSE, onefile = FALSE)
  if ( type == "pdf" )
    grDevices::pdf(file.name, width=width, height=height, bg=bg)
  if ( type == "tiff" )
    grDevices::tiff(file.name, width=width, height=height, units="in", res=res)
  if ( type == "svg" )
    grDevices::svg(file.name, width=width, height=height)
  if ( type == "jpeg" )
    grDevices::jpeg(file.name, width=width, height=height, units="in", res=res)
}


#' Add sub-figure label to corners of plots
#'
#' copied from code by January Weiner 3rd at
#' \url{https://www.r-bloggers.com/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/},
#' originally published at
#' \url{https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/}
#' @param text figure label
#' @param region add to a figure, plot or device
#' @param pos position key, as in \code{legend}
#' @param cex label size
#' @param ... additional parameters for call to \code{text}
#' @export
figlabel <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)

  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}




