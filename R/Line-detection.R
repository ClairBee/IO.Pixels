
#' Identify bad rows or columns
#' 
#' Perform a 2d convolution over an image to enhance lines of bright or dim pixels.
#' @param im 2d image array to be convolved.
#' @param k.size Size of (square) kernel to use. Default is 5.
#' @param horizontal Boolean: look for horizontal lines (T) instead of vertical (F)? Default is F, look for vertical lines.
#' @param dim.lines Boolean: look for lines of dim pixels (T) instead of bright (F)? Default is F, enhance bright lines.
#' @return Convolved image.
#' @export
#' @examples
#' zz <- convolve.lines(pw.m[,,"black", "160430"], k.size = 5)
#' 
convolve.lines <- function(im, k.size = 5, horizontal = F, dim.lines = F) {
    
    # define kernel to use in convolution
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)), rep(k.size - 1, k.size), rep(-1,k.size * floor(k.size / 2))), ncol = k.size)
    
    # if looking for horizontal lines, transpose kernel
    if (horizontal) {k <- t(k)}
    
    # if looking for dim lines, invert
    if (dim.lines) {k <- -k}
    
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



#' Summarise dimensions of line segments identified
#' 
#' After identifying potential bright/dim lines/rows using \link{\code{find.lines}}, produce a table of summary information about each column found.
#' @param im 2d image array after convolution, thresholding and smoothing
#' @param horizontal Boolean: summarise over rows (T) or columns (F)? Default is columns.
#' @return Data frame containing maximum and minimum row identified as part of a line segment, and the range and proportion of each column covered by line segments.
#' @export
#' 
summarise.lines <- function(im, horizontal = F, midline = 992.5) {
    
    xy <- data.frame(xyFromCell(m2r(im), which(getValues(m2r(im)) > 0)))
    
    if (horizontal) {
        df <- ddply(xy, .(row = y), summarise, xmin = min(x), xmax = max(x),
                    filled = length(x), range = max(x) - min(x) + 1,
                    cover = filled / range)
    } else {
        df <- rbind(ddply(xy[xy$y > midline,], .(col =x), summarise, panel = "U", ymin = min(y), ymax = max(y),
                          filled = length(x), range = max(y) - min(y) + 1,
                          cover = filled / range),
                    ddply(xy[xy$y < midline,], .(col = x), summarise, panel = "L", ymin = min(y), ymax = max(y),
                          filled = length(x), range = max(y) - min(y) + 1,
                          cover = filled / range))
    }

    return(df)
}


#' Filter lines to retain only long ones
#' 
#' @param im 2d image array after convolution, thresholding and smoothing
#' @param horizontal Boolean: summarise over rows (T) or columns (F)? Default is columns.
#' @param filter.at List of numeric filtering parameters: 'cover' is minimum proportion of line that must be filled to qualify as a long line (default: 0.5), 'filled' is minimum length (default: 20). 
#' @return image array with lines marked with separate indices for identification
#' @export
#' 
filter.lines <- function(im, filter.at = list(cover = 0.5, filled = 20), edges = 20, horizontal = F, midline = 992.5) {
    
    df <- summarise.lines(im, horizontal = horizontal)
    
    df <- df[df$cover > filter.at$cover & df$filled > filter.at$filled,]
    
    # if line ends are within 10px of panel centre or panel edge, reset to the centre/edge
    # (panel edges are cropped by convolution)
    df$ymin[df$ymin %in% (ceiling(midline) + c(0:9))] <- ceiling(midline)
    df$ymin[df$ymin %in% (1 + c(0:9))] <- 1
    
    df$ymax[df$ymax %in% (floor(midline) - c(0:9))] <- floor(midline)
    df$ymax[df$ymax %in% (nrow(im) - c(0:9))] <- nrow(im)
    
    new.im <- array(0, dim = dim(im))
    
    if (nrow(df) > 0) {
        for (i in 1:nrow(df)) {
            new.im[df$col[i], df$ymin[i]:df$ymax[i]] <- i
        }
    }

    # ignore any lines that are within 20px of a panel edge
    
    
    return(new.im)
}



#' Find lines in image
#' 
#' Use a combination of convolution & thresholding to identify consistently bright or dim columns or rows in an image
#' @param im 2d image array to be convolved.
#' @param k.size Size of (square) kernel to use. Default is 5.
#' @param threshold.at Numeric value at which to threshold convolved data. Default is 5500, which will highlight rows that are approximately 300 grey values higher than both of their neighbours in the direction of interest.
#' @param horizontal Boolean: look for horizontal lines (T) instead of vertical (F)? Default is F, look for vertical lines.
#' @param dim.lines Boolean: look for lines of dim pixels (T) instead of bright (F)? Default is F, enhance bright lines.
#' @param sm.size Length of (linear) kernel to use. Default is 11.
#' @param min.length Integer: minimum number of pixels to accept as a line after smoothing.
#' @return Image array with lines marked
#' @export
#' 
find.lines <- function(im, k.size = 5, threshold.at = 5500, sm.size = 11, min.length = 6, midline = 992.5,
                       filter.at = list(cover = 0.5, filled = 20),
                       horizontal = F, dim.lines = F) {
    
    # perform specified convolution
    conv <- convolve.lines(im, k.size = k.size, horizontal = horizontal, dim.lines = dim.lines)
    
    # threshold resulting image
    th <- threshold(conv, level = threshold.at)
    
    # smooth & repair lines
    sm <- smooth.lines(th, sm.size = sm.size, horizontal = horizontal, min.length = min.length)
    
    # filter out short line segments and return image of long segments identified
    lines <- filter.lines(sm, filter.at = filter.at, midline = midline)
    
    return(lines)
}