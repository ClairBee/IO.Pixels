
#' Plot given coordinates
#'
#' Plot pixel map
#' @param px Coordinates of pixels to be plotted
#' @export
#' @examples
#' pixel.plot(which(pw.m[,,"white"] < 15000, arr.ind = T))
#' 
pixel.plot <- function(px, xlim = c(0,2048), ylim = c(0,2048), pch = 15, panels = F, cex = 0.4,
                       main = "", xlab = "", ylab = "", ...) {
    
    plot(px[,1:2], asp = T, 
         xlim = xlim, ylim = ylim, pch = pch, main = main, xlab = xlab, ylab = ylab, cex = cex, ...)
}


#' Colours for plotting bad pixels
#' 
#' Support function: colour scheme for given pixel types
#' @return List of colours
#' @export
#' 
px.cols <- function() {
    c("magenta3", "black", "purple", "red", "orange", "gold", 
       "blue", "skyblue", "green1", "green3", "cyan1", "cyan3")
}