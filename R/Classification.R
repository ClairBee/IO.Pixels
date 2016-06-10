
#' Return coordinates of edge pixels
#' 
#' For a given image array, return the coordinates of pixels lying within a given distance of the edge
#' @details Edge pixels are applied to the outer edge of the image array, regardless of padding.
#' @param im Single-layer image array (needed to obtain dimensions)
#' @param edge.width Integer: number of pixels to treat as panel edges. Default is 10.
#' @return Data frame containing coordinates of edge pixels
#' @export
#'  
edge.px <- function(im, edge.width = 10) {
    
    im.width <- dim(im)[1]; im.height <- dim(im)[2]
    
    x <- c(sort(rep(1:edge.width, im.height)),                               # left edge
           sort(rep((im.width - edge.width + 1):im.width, im.height)),       # right edge
           rep(1:im.width, edge.width * 2))                                  # top & bottom edges
    
    y <- c(rep(1:im.height, edge.width * 2),
           sort(rep(1:edge.width, im.width)),
           sort(rep((im.height - edge.width + 1):im.height, im.width)))
    
    data.frame(row = x, col = y)
}


#' Identify pixels with no x-ray response
#' 
#' Identify pixels whose behaviour in an illuminated image is the same as in a dark image (ie. pixels that exhibit no response to source)
#' @param bright.image Image array containing an image with x-ray exposure, in which non-responsive pixels are to be found.
#' @param dark.image Image array containing a dark image (no x-ray exposure), used to calibrate 'normal' behaviour of pixels when not exposed to an x-ray source.
#' @param limit Quantile of black pixel population to be used as cutoff for 'normal' behaviour. Default is 0.0005, capturing central 99.9% of distribution.
#' @return Matrix of coordinates of non-responsive pixels
#' @export
#' @examples
#' qq <- no.response(160314, limit = 0.001)
no.response <- function(bright.image, dark.image, limit = 0.0005, exclude.edge = 10) {

    nc <- ncol(bright.image)
    
    bn <- qJohnson(c(limit, 1-limit), JohnsonFit(dark.image[!is.na(dark.image)]))
    
    un <- rbind(which(matrix(findInterval(bright.image, bn), ncol = nc) == 1, arr.ind = T))
    
    # filter out any pixels that fall in excluded edge region
    un[un[,1] > exclude.edge & 
           un[,2] > exclude.edge & 
           un[,1] < ncol(bright.image) - exclude.edge &
           un[,2] < nrow(bright.image) - exclude.edge, ]
    
    return(un)
}


# Get bright & dim pixels
#' For a single image array, return a list of thresholds to categorise very bright, bright, slightly bright, slightly dim, dim, and very dim pixels. Edge pixels are ignored for the purpose of this calculation.



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
bad.px.limits <- function(im, inner = "Johnson", edge.width = 10) {
    
    edge.width <- edge.width + 1
    x <- ncol(im); y <- nrow(im)
    inner <- tolower(inner)
    
    # remove edge pixels from threshold calculation
    im <- im[edge.width:(x - edge.width), edge.width:(y - edge.width)]
    
    # remove NA pixels
    im <- im[!is.na(im)]
    
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



#' Identify locally bright pixels
#' 
#' Identify locally bright pixels in an image
#' @details Requires existence of image array \code{pw.m} and median-differenced array \code{md}
#' @param dt Integer or string (format yymmdd) giving date to identify bright pixels in.
#' @return Data frame containing coordinates of dead pixels, with a factor indicating the classification
#' @export
#' 
locally.bright.px <- function(dt, setting) {
    
    dt <- toString(dt)
    setting <- tolower(setting)
    th <- mad(pw.m[,,setting, dt], na.rm = T) * 2
    
    if (length(which(md[,, setting, dt] > th)) == 0) {
        return(NULL)
    } else {
        return(data.frame(which(md[,, setting, dt] > th, arr.ind = T), type = "l.bright"))
    }
}


#' Identify locally dim pixels
#' 
#' Identify locally dim pixels in an image
#' @details Requires existence of image array \code{pw.m} and median-differenced array \code{md}
#' @param dt Integer or string (format yymmdd) giving date to identify dim pixels in.
#' @return Data frame containing coordinates of dead pixels, with a factor indicating the classification
#' @export
#' 
locally.dim.px <- function(dt, setting) {
    
    dt <- toString(dt)
    setting <- tolower(setting)
    th <- mad(pw.m[,,setting, dt], na.rm = T) * - 2
    
    if (length(which(md[,, setting, dt] < th)) == 0) {
        return(NULL)
    } else {
        return(data.frame(which(md[,, setting, dt] < th, arr.ind = T), type = "l.dim"))
    }
}