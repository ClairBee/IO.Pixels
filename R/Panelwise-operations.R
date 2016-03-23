
#' Split image array into subpanels
#'
#' Take a single image layer and return 32-layer labelled subpanel array to allow easier computation over panels.
#' @param data Single-layer 1996x1996 image array
#' @param left.pad Integer numer of pixels cropped from left-hand edge of image. Default is 2
#' @param upper.pad Integer numer of pixels cropped from upper edge of image. Default is 20
#' @export
#' @examples
#' load.images(150828, "black")
#' pw.m <- pixelwise.mean(b.150828)
#' 
#' zz <- subpanels(pw.m)
#' pixel.image(zz[,,"U8"], break.levels = sd.levels(b.150828))
#' 
subpanels <- function(data, left.pad = 2, upper.pad = 20) {
    
    # only works on single layer - unwieldy when larger
    
    panel.names <- apply(cbind(c(rep("U", 16), rep("L", 16)),
                               rep(c(1:16), 2)), 
                         1, paste, collapse = "")
    
    # pad original image with NA to bring to 2048 x 2048 square
    n <- array(NA, dim = c(2048, 2048))
    n[left.pad + c(1:1996), 2048-1996-upper.pad + c(1:1996)] <- data
    
    # create array to hold split data
    m <- array(dim = c(128, 1024, 32), 
               dimnames = list(dimnames(data)[[1]], dimnames(data)[[2]], panel.names))

    # extract panels
    for (i in 1:16) {
        m[ , , i] <- n[(128 * (i-1)) + c(1: 128), 1025:2048]        # upper row
        m[ , , i + 16] <- n[(128 * (i-1)) + c(1: 128), 1:1024]      # lower row
    }
    return(m)
}


#' Merge subpanels into image array
#' 
#' Take a subpanel image array and recombine into single panel image
#' @param data 128x1024x32 array of subpanel images
#' @param left.pad Integer numer of pixels cropped from left-hand edge of image. Default is 2
#' @param upper.pad Integer numer of pixels cropped from upper edge of image. Default is 20
#' @export
#' @examples
#' 
#' zz <- split.panels(pw.m)
#' qq <- join.panels(zz)
#' all(pw.m == qq)               # TRUE
join.panels <- function(data, left.pad = 2, upper.pad = 20) {
    m <- array(dim = c(2048, 2048), dimnames = list(dimnames(data)[[1]], dimnames(data)[[2]]))
    
    # extract panels
    for (i in 1:16) {
        m[(128 * (i-1)) + c(1: 128), 1025:2048] <- data[ , , i]         # upper row
        m[(128 * (i-1)) + c(1: 128), 1:1024] <- data [ , , i + 16]      # lower row
    }
    
    # return valid 1996x1996 region
    # (can't simply remove NA - after smoothing, NA may exist within panel)
    return(m[left.pad + c(1: 1996), 2048-1996-upper.pad + c(1:1996)])
}