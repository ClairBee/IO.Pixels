
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
#' @param min.segment Integer: minimum number of pixels to accept as a line after smoothing.
#' @return Smoothed and thresholded image
#' @export
#' 
smooth.lines <- function(im, sm.size = 11, horizontal = F, min.segment = 6) {
    
    # define kernel to be used in filtering
    k <- matrix(rep(1, sm.size), ncol = 1)
    
    # if looking for horizontal lines, transpose kernel
    if (horizontal) {k <- t(k)}
    
    # perform convolution
    sm <- r2m(focal(m2r(im), k))
    
    # threshold to remove short line fragments
    threshold(sm, level = min.segment - 0.5)
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
    
    # convert to raster & clump to get line IDs
    lines <- clump(m2r(im), dir = 4)
    
    xy <- data.frame(xyFromCell(lines, which(!is.na(getValues(lines)))),
                     id = getValues(lines)[!is.na(getValues(lines))])
    
    if (horizontal) {
        # labels are incorrect, but much easier to work with same column names
        df <- ddply(xy, .(id, col = y), summarise, 
                    ymin = min(x), ymax = max(x), length = length(x))
    } else {
        df <- rbind(ddply(xy[xy$y > midline,], .(id, col = x), summarise, panel = "U",
                          ymin = min(y), ymax = max(y), length = length(x)),
                    ddply(xy[xy$y < midline,], .(id, col = x), summarise, panel = "L",
                          ymin = min(y), ymax = max(y), length = length(x)))
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
filter.lines <- function(im, edges = 10, min.line = 10, horizontal = F, midline = 992.5, trim.ends.by = 2) {
    
    df <- summarise.lines(im, horizontal = horizontal, midline = midline)
    
    err <- "Still to be recoded... no horizontal lines when checking manually, so haven't prioritised automation"
    # remove line segments within 10 columns of panel edge
    if (horizontal) {
        df <- df[df$col > edges & df$col < dim(im)[2] - edges,]
    } else {
        df <- df[df$col > edges & df$col < dim(im)[1] - edges,]
    }
    
    # remove line segments less than 10 in length
    # could relate to kernel parameters - what is likely result of convolution & smoothing?
    df <- df[df$length > min.line,]
    
    # if line ends are within 10px of panel centre or panel edge, reset to the centre/edge
    # (panel edges are cropped by convolution)
    # could relate this to smoothing parameters - what is distance likely to be?

        df$ymin[df$ymin %in% (ceiling(midline) + c(0:9))] <- ceiling(midline)
        df$ymin[df$ymin %in% (1 + c(0:9))] <- 1
        
        df$ymax[df$ymax %in% (floor(midline) - c(0:9))] <- floor(midline)
        df$ymax[df$ymax %in% (nrow(im) - c(0:9))] <- nrow(im)


    # recalculate segment lengths based on new endpoints
        df$length <- df$ymax - df$ymin + 1
    
    # remove any segments less than 20px in length - UNLESS associated with longer segment
    long <- df[df$length > 20,]
    df <- df[df$col %in% long$col,]
    
    # trim line ends to remove 'smoothed' area left over from convolution
    df$ymin[!(df$ymin %in% c(ceiling(midline),1))] <- df$ymin[!(df$ymin %in% c(ceiling(midline),1))] + trim.ends.by

    df$ymax[!(df$ymax %in% c(floor(midline), dim(im)[2]))] <- df$ymax[!df$ymax %in% c(floor(midline), dim(im)[2])] - trim.ends.by

    # create new image array containing lines identified
    new.im <- array(0, dim = dim(im))
    
    if (nrow(df) > 0) {
        for (i in 1:nrow(df)) {
            new.im[df$col[i], df$ymin[i]:df$ymax[i]] <- i
        }
    }
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
#' @param min.segment Integer: minimum number of pixels to accept as a line segment after smoothing.
#' @param min.line Integer: minimum number of pixels to accept as a single line during filtering.
#' @return Image array with lines marked
#' @export
#' 
find.lines <- function(im, k.size = 5, threshold.at = 5500, sm.size = 11, min.segment = 6, midline = 992.5,
                       edges = 10, min.line = 10, trim.ends = T,
                       horizontal = F, dim.lines = F) {

    # perform specified convolution
    conv <- convolve.lines(im, k.size = k.size, horizontal = horizontal, dim.lines = dim.lines)
    
    # threshold resulting image
    th <- threshold(conv, level = threshold.at)
    
    # smooth & repair lines
    sm <- smooth.lines(th, sm.size = sm.size, horizontal = horizontal, min.segment = min.segment)
    
    # filter out short line segments and return image of long segments identified
    lines <- filter.lines(sm, edges = 10, min.line = 10, horizontal = horizontal, 
                          midline = midline, trim.ends.by = floor(k.size / 2))
    
    return(lines)
}