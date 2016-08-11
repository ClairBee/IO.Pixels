
#' Find dark lines in an image
#' 
#' Identify pixels below a certain threshold and clump with adjacent dark pixels to find long columns of dark pixels. (To find rows, the bright image should be transposed)
#' @param bright.image Bright image array in which dark pixels are to be identified
#' @param dark.limit Integer grey-value below which a pixel will be classified as dark. Default is 25000
#' @param min.length Integer minimum length of a line segment to be considered a line (rather than a cluster). Default is 5
#' @param midline Y-coordinate of horizontal midline across panel, where upper and lower subpanels are electrically separated. Default is 1024.5; if no midline exists, use midline = NA
#' @return Data frame of x and y coordinates of identified dark lines, with line type and feature type
#' @export
#' @examples
#' dc <- dark.lines(pw.m[,,"white", "130613"])         # look for dark columns
#' dr <- dark.lines(t(pw.m[,,"white", "130613"]))      # look for dark rows
#' 
dark.lines <- function(bright.image, dark.limit = 25000, min.length = 5, midline = 1024.5) {
    
    # threshold image to identify dark pixels
    dpx <- bright.image < dark.limit
    
    # clump dark pixels, get coordinates
    lines <- clump(m2r(dpx), dir = 4)
    xy <- data.frame(xyFromCell(lines, which(!is.na(getValues(lines)))), 
                     id = getValues(lines)[!is.na(getValues(lines))])
    
    # summarise by column & clump ID. Filter out short segments
    dl <- ddply(xy, .(id, col = x), summarise, ymin = min(y), ymax = max(y), length = ymax - ymin + 1)
    dl <- dl[dl$length > min.length,]
    
    if (nrow(dl) == 0) return (NULL)
    
    # return coordinates of selected lines
    x <- c(unlist(unname(apply(dl, 1, function(ll) rep(unname(ll["col"]), ll["length"])))))
    y <- c(unlist(unname(apply(dl, 1, function(ll) seq(ll["ymin"] * 1, ll["ymax"] * 1)))))
    
    df <- data.frame(x, y, type = "dark", f.type = "line.body", stringsAsFactors = F)
    
    # if no midline, use top of panel as line end
    if (is.na(midline)) midline <- nrow(bright.image)
    
    # set line ends as line.root feature type (if line ends are visible)
    for (r in 1:nrow(dl)) {
        l.end <- as.numeric(dl[r, c("ymin", "ymax")][which.max(abs(dl[r, c("ymin", "ymax")] - 1024.5))])
        if (!l.end %in% range(which(!is.na(bright.image[as.integer(dl[r, "col"]),])))) {
            df[df$x == as.numeric(dl[r, "col"]) & df$y == l.end, "f.type"] <- "line.root"
        }
    }
    
    return(df)
}



#' Find dim or bright lines in an image
#' 
#' Use convolution and filtering to identify vertical segments of bright or dim pixels.
#' @param k.size Size of (square) kernel to use. Default is 5.
#' @param threshold.at Numeric value at which to threshold convolved data. Default is 5500, which will highlight rows that are approximately 300 grey values higher than both of their neighbours in the direction of interest.
#' @param min.length Integer: minimum number of pixels to accept as a line segment after smoothing. Default is 10.
#' @param midline Y-coordinate of horizontal midline across panel, where upper and lower subpanels are electrically separated. Default is 1024.5; if no midline exists, use midline = NA
#' @return Data frame containing x and y coordinates of identified lines, and classification of each segment as either a line root, line body, or channel offset.
#' @export
#' @examples
#' 
offset.lines <- function(md, k.size = 5, threshold.at = 5500, min.length = 10, midline = 1024.5) {
    
    # convolution kernel to search for vertical lines
    k <- matrix(c(rep(-1, k.size * floor(k.size / 2)),
                  rep(k.size - 1, k.size), 
                  rep(-1,k.size * floor(k.size / 2))),
                ncol = k.size)
    
    # perform convolution
    conv <- r2m(focal(m2r(md), k))
    
    # threshold (retaining distinction between bright & dim lines)
    th <- conv > threshold.at
    th[which(conv < - threshold.at, arr.ind = T)] <- 2
    
    # find line segments by clumping adjacent pixels
    lines <- clump(m2r(th), dir = 4)
    xy <- data.frame(xyFromCell(lines, which(!is.na(getValues(lines)))), 
                     id = getValues(lines)[!is.na(getValues(lines))])
    xy$type <- th[as.matrix(xy[,1:2])]
    
    # summarise line segments per clump, column and type
    dl <- ddply(xy, .(id, x, type), summarise, ymin = min(y), 
                ymax = max(y), length = length(x), density = length / (ymax - ymin + 1))
    
    # filter out short line segments, or those which are less than half-covered
    dl <- dl[dl$length > min.length & dl$density > 0.5, ]

    # return NULL if no segments found
    if (nrow(dl) == 0) return(NULL)
    
    #-------------------------------------------------
    
    # otherwise merge any line segments that fall in same channel
    if (is.na(midline)) {
        dl <- ddply(dl, .(x), summarise, y.min = min(ymin), y.max = max(ymax))
    } else {
        dl <- rbind(ddply(dl[dl$ymin > 1024.5,], .(x), summarise, ymin = min(ymin), ymax = max(ymax)),
                    ddply(dl[dl$ymin < 1024.5,], .(x), summarise, ymin = min(ymin), ymax = max(ymax)))
    }
    
    #identify most likely breakpoint & get mean value for each segment
    
    cp <- rbind.fill(apply(dl, 1, 
                           function(ll) {
                               # extract relevant line segment
                               if (!is.na(midline)) {
                                   p.adj <- as.integer(floor(midline) * (ll["ymax"] > midline))
                                   vv <- md[as.integer(ll["x"]), c(1:floor(midline)) + p.adj]
                               } else {
                                   p.adj <- 0
                                   vv <- md[ll["x"],]
                               }
                               
                               # store offset for NA values padding edge of array
                               vv.adj <- (min(which(is.na(vv))) == 1) * max(which(is.na(vv)))
                               l <- length(vv[!is.na(vv)])
                               
                               # identify changepoints (at most one changepoint)
                               cm <- cpt.mean(vv[!is.na(vv)], class = F)[1]
                               
                               if (cm == 1) {
                                   # no changepoints found - whole-column defect
                                   data.frame(x = unlist(rep(ll["x"], l)),
                                              y = c(1:l) + p.adj + vv.adj,
                                              type = c("dim", "bright")[(mean(vv, na.rm = T) > 0) + 1],
                                              f.type = "channel",
                                              stringsAsFactors = F)
                               } else {
                                   # changepoint found - include both segments
                                   px <- rbind(data.frame(x = unlist(rep(ll["x"], cm)),
                                                          y = c(1:cm) + p.adj + vv.adj,
                                                          type = c("dim", "bright")[(mean(vv[1:(cm + vv.adj)], na.rm = T) > 0) + 1],
                                                          f.type = "line.body",
                                                          stringsAsFactors = F),
                                               data.frame(x = unlist(rep(ll["x"], l-cm)),
                                                          y = c((cm + 1):l) + p.adj + vv.adj,
                                                          type = c("dim", "bright")[(mean(vv[(cm + vv.adj):length(vv)], na.rm = T) > 0) + 1],
                                                          f.type = "line.body",
                                                          stringsAsFactors = F))
                                   px[px$y == cm + p.adj + vv.adj, "f.type"] <- "line.root"
                                   px
                               }
                           }))
    
    return(cp)
}





####################################################################################################





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
summarise.lines <- function(im, horizontal = F, midline = 1024.5) {
    
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
filter.lines <- function(im, edges = 10, min.line = 10, horizontal = F, midline = 1024.5, trim.ends.by = 2) {
    
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

    if (horizontal) {
        # snap to panel RHS
        df$ymax[df$ymax %in% (ncol(im) - c(0:9))] <- ncol(im)
    } else {
        # snap to midline
        df$ymin[df$ymin %in% (ceiling(midline) + c(0:9))] <- ceiling(midline)
        df$ymax[df$ymax %in% (floor(midline) - c(0:9))] <- floor(midline)
        
        # snap to panel top
        df$ymax[df$ymax %in% (nrow(im) - c(0:9))] <- nrow(im)
    }
    # snap to panel LHS/bottom
    df$ymin[df$ymin %in% (1 + c(0:9))] <- 1
    
    # trim line ends to remove 'smoothed' area left over from convolution
    if (horizontal) {
        df$ymin <- df$ymin + trim.ends.by
        df$ymax <- df$ymax-trim.ends.by
    } else {
        df$ymin[!(df$ymin %in% c(ceiling(midline),1))] <- df$ymin[!(df$ymin %in% c(ceiling(midline),1))] + trim.ends.by
        df$ymax[!(df$ymax %in% c(floor(midline), dim(im)[2]))] <- df$ymax[!df$ymax %in% c(floor(midline), dim(im)[2])] - trim.ends.by
    }
    
    # recalculate segment lengths based on new endpoints
    df$length <- df$ymax - df$ymin + 1
    
    # remove any segments less than 20px in length - UNLESS associated with longer segment
    long <- df[df$length > 20,]
    df <- df[df$col %in% long$col,]
    
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
find.lines <- function(im, k.size = 5, threshold.at = 5500, sm.size = 11, min.segment = 6, midline = 1024.5,
                       edges = 10, min.line = 10, trim.ends = T,
                       horizontal = F, dim.lines = F) {

    # remove padding
    ar <- apply(which(!is.na(im), arr.ind = T), 2, range)
    org.dims <- list(im = dim(im), midline = midline)
    im <- im[ar[1,"row"]:ar[2,"row"], ar[1,"col"]:ar[2,"col"]]
    midline <- midline - (ar[1,"col"] - 1)
    
    # perform specified convolution
    conv <- convolve.lines(im, k.size = k.size, horizontal = horizontal, dim.lines = dim.lines)
    
    # threshold resulting image
    th <- threshold(conv, level = threshold.at)
    
    # smooth & repair lines
    sm <- smooth.lines(th, sm.size = sm.size, horizontal = horizontal, min.segment = min.segment)
    
    # filter out short line segments and return image of long segments identified
    lines <- filter.lines(sm, edges = edges, min.line = min.line, horizontal = horizontal, 
                          midline = midline, trim.ends.by = floor(k.size / 2))
    
    # re-pad image
    l <- array(dim = org.dims$im)
    l[ar[1,"row"]:ar[2,"row"], ar[1,"col"]:ar[2,"col"]] <- lines
    
    return(l)
}