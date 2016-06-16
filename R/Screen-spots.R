
#' Identify pixels covered by spots on screen
#' 
#' Find dim spots of a given size or larger in a single bright image, and return their x and y coordinates
#' @param bright.image Image array containing an image with x-ray exposure, in which non-responsive pixels are to be found.
#' @param smooth.span Smoothing span to be used for Lowess. Default is 1/5 (this is a lower proportion than the standard default, since the aim is to 'skim over' any deviations from a smooth curve, not to fit them)
#' @param min.diam Integer: minimum diameter of a dim spot, in pixels (used to set size of circular structuring element for morphological closing). Default is 5.
#' @param midline Integer: y-coordinate of horizontal panel midline, dividing electrically separates panels. Per-column Lowess smoothing stops at the midline.
#' @param enlarge Boolean: should identified spots be dilated (enlarged) by the structuring element to smooth the edges and cover the largest area possible? default is T.
#' @param auto.threshold Boolean: should screen spots be thresholded using k-means clustering (T) or by thresholding at the mean value - 3sd (F)? Default is T (auto-thresholding).
#' @return Matrix of XY coordinates of pixels identified as probable screen spots
#' @export
#' @examples
#' bp <- data.frame(screen.spots(pw.m[,,"white", "160430"]), type = "screen.spot")
#'
screen.spots <- function(bright.image, smooth.span = 1/5, min.diam = 5, midline = 992.5, enlarge = F, auto.threshold = T, ignore.edges = 10) {
    
    im.dims <- dim(bright.image)
    
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    # apply lowess smoothing
    # don't use default smoothing in this case - don't want smoother to be drawn into dips, so use lower proportion
    smoo <- lowess.per.column(bright.image, midline = midline, span = smooth.span)
    
    res <- bright.image - smoo
    
    # flatten further by setting brighter pixels to mean value 
    res[res > mean(res)] <- mean(res)
    
    # dilate resulting image
    dilated <- dilate(res, sk)
    
    # erode resulting image (complete morphological closing)
    eroded <- erode(dilated, sk)
    
    # use k-means thresholding to identify spots
    # use 1-thresholded value to assign 1 to spots, 0 to background
    if (auto.threshold) {
        dim.sp <- array(1, dim = im.dims) - threshold(eroded, method = "kmeans")
    } else {
        dim.sp <- array(1, dim = im.dims) - threshold(eroded, mean(eroded) - 3*sd(eroded))
    }
    
    # if spot identification should be as large as possible, 
    if (enlarge) {dim.sp <- dilate(dim.sp, sk)}
    
    # remove any spots whose midpoint lies within the edge region
    if (ignore.edges > 0) {
        blobs <- clump(m2r(dim.sp), dir = 4)
        
        sc <- ddply(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))),
                               id = getValues(blobs)[!is.na(getValues(blobs))]),
                    .(id), summarise, 
                    x.midpoint = round(mean(x),0), y.midpoint = round(mean(y),0),
                    size = length(x))
        
        sc$n.id <- sc$id
        sc$n.id[(sc$y.midpoint >= im.dims[2] - ignore.edges) | (sc$y.midpoint <= ignore.edges) | 
                    (sc$x.midpoint >= im.dims[1] - ignore.edges) | (sc$x.midpoint <= ignore.edges)] <- NA
        blobs <- subs(blobs, sc[,c("id", "n.id")])
        
        dim.sp <- bpx2im(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))), type = 1))
    }
    
    if (sum(dim.sp == 1) == 0) {
        return(NULL)
    } else {
        return(which(dim.sp == 1, arr.ind = T))
    }
}



# identify individual screen spots from list of coordinates
#' 
#' Convert a list of screen spot coordinates into an image array and summary table
#' @param ss.coords Matrix of XY coordinates of identified screen spots, as returned by \code{\link{screen.spots}}
#' @param im.dims Dimensions of image on which to display coordinates
#' @param ignore.edges Integer: width of border around panel within which pixels are assumed to be dim because of their proximity to the edge, rather than because of a spot on the screen. Spots whose midpoint falls within this border will be removed. Default is 10.
#' @return List containing an image array of the specified dimensions, and a data frame summarising each spot's midpoint and size.
#' @export
#' @examples
#' image(screen.spot.clusters(screen.spots(pw.m[,,"white", "160430"]))$spots)
#' 
screen.spot.clusters <- function(ss.coords, im.dims = c(1996, 1996), ignore.edges = 10) {
    
    ss.im <- bpx2im(data.frame(ss.coords, type = 1), im.dim = im.dims)
    
    #------------------------------------------------------------------------------
    # convert to raster & clump values to identify individual spots
    # (values need to be reordered to maintain image orientation)
    blobs <- clump(m2r(ss.im), dir = 4)
    
    # summarise & remove any clusters that are too close to an edge
    # check location of cluster (looking for edge clusters)
    sc <- ddply(data.frame(xyFromCell(blobs, which(!is.na(getValues(blobs)))),
                           id = getValues(blobs)[!is.na(getValues(blobs))]),
                .(id), summarise, 
                x.midpoint = round(mean(x),0), y.midpoint = round(mean(y),0),
                size = length(x))
    
    sc$n.id <- sc$id
    sc$n.id[(sc$y.midpoint >= im.dims[2] - ignore.edges) | (sc$y.midpoint <= ignore.edges) | 
                (sc$x.midpoint >= im.dims[1] - ignore.edges) | (sc$x.midpoint <= ignore.edges)] <- NA
    
    blobs <- subs(blobs, sc[,c("id", "n.id")])
    
    # return array with numbered blobs
    return(list(spots = blobs, summary = sc[!is.na(sc$n.id),1:4]))
}