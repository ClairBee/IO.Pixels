
#' Apply shading correction
#' 
#' Calculate a shading-corrected image from a set of blanks
#' @param im Set of pixelwise means for black, grey and white acquisitions
#' 
#' @export
#' 
shading.corrected <- function(im, fix.inf = T) {
    sc <- 60000 * (im[,,"grey"] - im[,,"black"]) / (im[,,"white"] - im[,,"black"])
    
    if (fix.inf) {
        sc[is.infinite(sc)] <- NA
    }
    sc
}