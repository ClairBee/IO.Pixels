
#' Identify pixels covered by spots on screen
#' 
#' Find dim spots of a given size or larger in a single bright image, and return their x and y coordinates
#' @param im Image array containing an image with x-ray exposure, in which non-responsive pixels are to be found.
#' @param smooth.span Smoothing span to be used for Lowess. Default is 1/15 (this is a lower proportion than the standard default, since the aim is to 'skim over' any deviations from a smooth curve, not to fit them)
#' @param min.diam Integer: minimum diameter of a dim spot, in pixels (used to set size of circular structuring element for morphological closing). Default is 5.
#' @param midline Numeric: y-coordinate of horizontal panel midline, dividing electrically separates panels. Default is 1024.5.
#' @return Matrix of XY coordinates of pixels identified as probable screen spots
#' @export
#' @examples
#' bp <- data.frame(screen.spots(pw.m[,,"white", "160430"]), type = "screen.spot")
#'
screen.spots <- function(im, min.diam = 5, smooth.span = 1/15, midline = 1024.5, edge.crop = 10, coords = T) {
    
    # offset correction for upper vs lower panels (if midline exists)
    if (!is.na (midline)) {
        up <- apply(im[, floor(midline) + c(1:100)], 1, median, na.rm = T)
        lp <- apply(im[, floor(midline) + c(0:-100)], 1, median, na.rm = T)
        
        im[, ceiling(midline):dim(im)[[2]]] <- im[, ceiling(midline):dim(im)[[2]]] - (up - lp)
    }
    
    # Lowess smoothing over all columns
    # (faster than Loess & easier to apply in this format)
    res <- im - do.call("rbind", lapply(apply(im, 1, lowess, f = smooth.span), "[[", 2))
    med.res <- median(res, na.rm = T)
    
    # truncate residual values at median
    tr <- res
    tr[res > med.res] <- med.res
    
    # morphological closing
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    cl <- closing(tr, sk)
    
    # pad with median residual value
    cl[is.infinite(cl) | is.na(cl)] <- med.res
    
    # threshold at 1 SD below median
    th <- threshold(cl, method = "literal", level = - sd(res, na.rm = T))
    
    # remove any pixels that lie at edge of active area
    th[, apply(which(!is.na(res), arr.ind = T), 2, range)[,2]] <- 1 # residual artefact at edge of smoothed area
    th[c(1 + c(0:edge.crop), dim(th)[1] - c(0:edge.crop)), ] <- 1    # panel edges
    th[, c(1 + c(0:edge.crop), dim(th)[1] - c(0:edge.crop))] <- 1    # panel edges
    
    # enlarge by eroding with same structuring element
    exp <- erode(th, sk)
    
    if (coords) {return(which(exp == 0, arr.ind = T))} else {return(exp)}
}

