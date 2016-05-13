
#' Get cutoffs for bright/dim pixels
#' 
#' For a single image array, return a list of thresholds to categorise very bright, bright, slightly bright, slightly dim, dim, and very dim pixels. Edge pixels are ignored for the purpose of this calculation.
#' @param im Two-dimensional image array containing pixel values to be thresholded
#' @return List of named numbers to use as cutpoints
#' @export
#' @examples
#' img <- pw.m[,,"white", "160314"]
#' table(findInterval(img, unlist(bad.px.limits(img))))
#' 
bad.px.limits <- function(im, inner = "Johnson") {
    
    inner <- tolower(inner)
    
    # remove edge pixels from threshold calculation
    im <- im[11:1985, 11:1985]
    
    JF <- JohnsonFit(im, moment = "quant")
    
    bright.s <- switch(inner,
                       "johnson" = qJohnson(0.9999, JF),
                       "sd2" = median(im) + 2 * sd(im))
    
    bright.v <- median(im) + abs((max(im) - median(im)) / 2)
    bright <- median(im) + abs((bright.v - median(im)) / 2)

    dim.s <-     switch(inner,
                        "johnson" = qJohnson(0.0001, JF),
                        "sd2" = median(im) - 2 * sd(im))
    
    dim.v <- median(im) - abs((median(im) - min(im)) / 2)
    dim <- median(im) - abs((median(im) - dim.v) / 2)

    list(dv = dim.v, dm = dim, ds = dim.s, bs = bright.s, bm = bright, bv = bright.v)
}


#' Identify bright and dim pixels
#' 
#' Use the thresholds returned by \link{\code{bad.px.limits}} to classify bright and dim pixels in an image.
#' @param im Two-dimensional image array containing pixel values to be classified
#' @param cropped Boolean: crop edge pixels before thresholding? Default is T.
#' @return Data frame containing coordinates of bright and dim pixels, with a factor indicating the classification
#' @export
#' @examples
#' bp <- get.dim.bright.px(res[,,"white", "160314"])
#' table(bp$type)
#' 
get.dim.bright.px <- function(im) {
    lim <- bad.px.limits(im)
    
    # check each category in turn
    bp <- rbind(data.frame(which(im > lim$bv, arr.ind = T), type = "v.bright"),
                data.frame(which(im > lim$bm, arr.ind = T), type = "bright"),
                data.frame(which(im > lim$bs, arr.ind = T), type = "s.bright"),
                data.frame(which(im < lim$dv, arr.ind = T), type = "v.dim"),
                data.frame(which(im < lim$dm, arr.ind = T), type = "dim"),
                data.frame(which(im < lim$ds, arr.ind = T), type = "s.dim"))
    
    # remove duplicates (most severe category will be retained)
    bp <- bp[!duplicated(bp[,1:2]),]
    return(bp)
}


#' Get bad pixel map for a given date
#' 
#' Use thresholding and absolute values to classify bad pixels in a single acquisition, using the thresholding limits set in \link{\code{bad.bx.limits}}.
#' @details Requires that the pixelwise mean array \code{pw.m} and a residual array \code{res} be available in the global environment.
#' @param dt Integer or string giving acquisition date (format yymmdd)
#' @export
#' @examples
#' bp <- bad.pixels("160314")
#' 
bad.px <- function(dt) {
    dt <- toString(dt)
    
    qq <- screen.spots.xy(dt)
    
    bp <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                             rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                    ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright", "s.bright", "screen spot", "edge", "v.dim", "dim", "s.dim"))),
                         c("row", "col", "type")),
                data.frame(no.response(dt), type = "no response"),
                data.frame(which(pw.m[, , "black", dt] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(pw.m[, , "white", dt] == 0, arr.ind = T), type = "dead"),
                screen.spots.xy(dt),
                get.dim.bright.px(res[, , "grey", dt]),
                get.dim.bright.px(res[, , "white", dt]))
    
    bp <- bp[order(bp$type),]
    bp[!duplicated(bp[,1:2]),]
}