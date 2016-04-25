
#' Identify dim spots on scree
#'
#' Find dim spots of a given size or larger, in a single white image
#' @details Requires that \code{pw.m} be available in the global environment, as by running \link{\code{load.pixel.means()}}
#' @param dt Integer/string specifying image date to be imported
#' @param smooth.span Smoothing span to be used for Lowess. Default is 1/5
#' @param min.diam Integer: minimum diameter of a dim spot, in pixels (used to set size of circular structuring element for morphological closing). Default is 5.
#' @param edge.width Integer: dim spots whose midpoint falls within this distance of the detector edge will be removed from the list of screen defects.
#' @export
#' @examples
#' load.pixel.maps()
#' ds <- screen.spots(160314)
#' 
screen.spots <- function(dt, smooth.span = 1/5, min.diam = 5, edge.width = 10) {
    dt <- toString(dt)
    im <- pw.m[,,"white", dt]
    
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")

    # apply lowess smoothing
    smoo <- lowess.per.column(im, span = smooth.span)
    res <- im - smoo
    
    # flatten further by setting brighter pixels to mean value 
    res[res > mean(res)] <- mean(res)
    
    # dilate resulting image
    dilated <- dilate(res, sk)
    
    # erode resulting image (complete morphological closing)
    eroded <- erode(dilated, sk)
    
    # use k-means thresholding to identify spots
    # use 1-thresholded value to assign 1 to spots, 0 to background
    dim <- array(1, dim = c(1996, 1996)) - threshold(eroded, method = "kmeans")
    
    #------------------------------------------------------------------------------
    # convert to raster & clump values to identify individual spots
    # (values need to be reordered to maintain image orientation)
    blobs <- clump(raster(t(dim[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5), dir = 4)
    
    # summarise & remove any clusters that are too close to an edge
    # check location of cluster (looking for edge clusters)
    sc <- ddply(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))),
                           id = getValues(blobs)[!is.na(getValues(blobs))]),
                .(id), summarise, 
                xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y),
                xm = round(mean(x),0), ym = round(mean(y),0))
    
    sc$n.id <- sc$id
    sc$n.id[(sc$ym >= 1996 - edge.width) | (sc$ym <= edge.width) | 
                (sc$xm >= 1996 - edge.width) | (sc$xm <= edge.width)] <- NA

    blobs <- subs(blobs, sc[,c("id", "n.id")])
    
    # return array with numbered blobs
    return(blobs)
}