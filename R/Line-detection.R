
#' Identify bad rows or columns
#' 
#' Perform a 2d convolution over an image to enhance lines of bright or dim pixels.
#' @param im 2d image array to be convolved.
#' @param k.size Size of (square) kernel to use. Default is 5.
#' @param horizontal Boolean: look for horizontal lines (T) instead of vertical (F)? Default is F, look for vertical lines.
#' @param dim Boolean: look for lines of dim pixels (T) instead of bright (F)? Default is F, enhance bright lines.
#' @return Convolved image.
#' @export
#' @examples
#' zz <- convolve.lines(pw.m[,,"black", "160430"], k.size = 5)
#' 
convolve.lines <- function(im, k.size = 5, horizontal = F, dim = F) {
    
    # define kernel to use in convolution
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    
    # if looking for horizontal lines, transpose kernel
    if (horizontal) {k <- t(k)}
    
    # if looking for dim lines, invert
    if (dim) {k <- -k}
    
    # perform convolution
    r2m(focal(m2r(im), k))
}



#' Smooth lines
#' 
#' Use a focal window to smooth line segments and remove short fragments.
#' @param im 2d image array to be convolved.
#' @param sm.size Length of (linear) kernel to use. Default is 11.
#' @param horizontal Boolean: look for horizontal lines (T) instead of vertical (F)? Default is F, look for vertical lines.
#' @param min.length Integer: minimum number of pixels to accept as a line after smoothing.
#' @return Smoothed and thresholded image
#' @export
#' 
smooth.lines <- function(im, sm.size = 11, horizontal = F, min.length = 6) {
    
    # define kernel to be used in filtering
    k <- matrix(rep(1, sm.size), ncol = 1)
    
    # if looking for horizontal lines, transpose kernel
    if (horizontal) {k <- t(k)}
    
    # perform convolution
    sm <- r2m(focal(m2r(im), k))
    
    # threshold to remove short line fragments
    threshold(sm, level = min.length - 0.5)
}