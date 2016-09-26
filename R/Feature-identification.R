
#' Get screen spots in linear/white images
#' 
#' Screen spot identification in linear/white images
#' @export
#' 
screen.spots <- function(im, smooth.span = 1/15, min.diam = 5, midline = 1024.5, match.width = 100) {
    
    # support function to adjust offset across midline
    adjust.midline.offset <- function(cc, midline = 1024.5, match.width = 100) {
        os.u <- median(cc[ceiling(midline) + c(0:match.width)], na.rm = T)
        os.l <- median(cc[floor(midline) - c(0:match.width)], na.rm = T)
        cc[ceiling(midline):length(cc)] <- cc[ceiling(midline):length(cc)] - (os.u - os.l)
        return(cc)
    }
    
    # get white image & residuals from linear regression
    w <- im[,,"white"]
    linear <- fit.w.lm(im)
    
    # define kernel for morphological operations
    sk <- shapeKernel(c(min.diam, min.diam), type = "disc")
    
    #----------------------------------------------------------
    
    if(!is.na(midline)) {
        l.adj <- t(apply(linear, 1, adjust.midline.offset, midline = midline, match.width = match.width))
        w.adj <- t(apply(w, 1, adjust.midline.offset, midline = midline, match.width = match.width))
    } else {
        l.adj <- linear
        w.adj <- w
    }
    
    # find spots in linear residuals
    {
        # lowess residuals along columns
        l.res <- l.adj - t(apply(l.adj, 1, function(cc) lowess(cc, f = smooth.span)$y))
        l.res.high <- l.res.low <- l.res
        
        # dim spots in linear residuals
        {
            l.res.low[l.res.low > 0 | is.na(l.res.low)] <- 0        # filter low residuals
            l.dim <- closing(l.res.low, sk)                         # morphological closing
            l.dim.th <- (l.dim < -mad(l.res, na.rm = T)) * 1        # threshold
            l.dim.th <- dilate(l.dim.th, sk)                        # slightly enlarge all spots
            l.dim.px <- which(l.dim.th > 0, arr.ind = T)            # get pixel coordinates
        }
        
        # bright regions in linear residuals
        {
            l.res.high[l.res.high < 0 | is.na(l.res.high)] <- 0
            l.bright <- opening(l.res.high, sk)
            l.bright.th <- (l.bright > mad(l.res, na.rm = T)) * 1
            l.bright.th <- dilate(l.bright.th, sk)
            l.bright.px <- which(l.bright.th > 0, arr.ind = T)
        }
    }
    
    # find spots in white image
    {
        # offset correction & lowess residuals along columns
        w.res <- w.adj - t(apply(w.adj, 1, function(cc) lowess(cc, f = smooth.span)$y))
        
        # filter residuals - retain only dim regions
        w.res.low <- w.res
        w.res.low[w.res.low > 0 | is.na(w.res.low)] <- 0
        
        w.dim <- closing(w.res.low, sk)                         # morphological closing
        w.dim.th <- (w.dim < -mad(w.res, na.rm = T)) * 1        # threshold
        w.dim.th <- dilate(w.dim.th, sk)                        # slightly enlarge all spots
        w.dim.px <- which(w.dim.th > 0, arr.ind = T)            # get pixel coordinates
    }
    
    list(nl.dim = l.dim.px, nl.bright = l.bright.px, w.dim = w.dim.px)
}


#' Label identified screen spots in pixel list
#' 
#' Given a pixel list and list of identified screen spots, return a pixel list with defects belonging to screen spots identified.
#' @param px Data frame containing pixel coordinates to be grouped
#' @param spot.list List of matrices of coordinates of identified screen spots, as returned by \link{\code{screen.spots}}
#' @return Original pixel list, with column 'f.type' now showing identified screen spots
#' @export
#' 
label.screen.spots <- function(px, spot.list) {
    
    spots <- unique(data.frame(do.call("rbind", sp[sapply(sp, nrow) > 0]), s.spot = T))
    
    # identify individual spots
    cc <- clump(px2r(spots), dir = 4)
    
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    xy <- merge(xy, px, by = c(1:2), all.x = T)
    
    # filter out any spots where no component pixels were identified as defects
    xy.s <- ddply(xy, .(id), summarise, defects = sum(!is.na(type)), chk = sum(is.na(type)), size = length(x))
    xy <- xy[xy$id %in% xy.s$id[xy.s$defects > 0],]
    
    npx <- merge(px, data.frame(xy[,1:2], s.spot = T), by = c(1:2), all.x = T)
    
    npx$f.type[npx$s.spot & npx$type %in% c("dim", "l.dim", "nl.dim", "spot.dim", NA)] <- "s.spot"
    
    return(npx[colnames(npx) %in% c(colnames(px), "f.type")])
}


#' Find column features in pixel map
#' 
#' Applies a linear convolution filter to identify vertical linear features in the pixel map.
#' @param px Data frame containing pixel coordinates to be grouped
#' @param kernel.size Edge length of square kernel to be used in convolution. Default is 5 (aims to pick up line segments of length 5 or greater) 
#' @param midline Location of panel midline, which electrically separates the upper and lower panels, and so cannot be crossed by column defects. Default is 1024.5; if the detector has no midline separation, enter NA.
#' @param min.length Minimum length of segment to retain after convolution. Default is 3 * kernel size.
#' @return Data frame px, with column 'f.type' labelling columns found as feature type 'line.c'
#' @export
#' 
find.columns <- function(px, kernel.size = 5, midline = 1024.5, min.length = 3 * kernel.size) {
    
    # if kernel size is even sides, add 1
    if(kernel.size %% 2 == 0) {kernel.size <- kernel.size + 1}
    
    # define kernel to use in convolution
    k <- matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
                  rep(kernel.size - 1, kernel.size), 
                  rep(-1, kernel.size * floor(kernel.size / 2))),
                ncol = kernel.size)
    
    # if midline given, split data & apply twice (row/col labels are inverted!)
    if(!is.na(midline)) {
        px.list <- list(px[px$col > midline,],
                        px[px$col < midline,])
    } else {
        px.list <- list(px)
    }
    
    lines <- lapply(px.list, 
                    function(c.px) {
                        
                        # convert px to raster, convolve with kernel, convert back to matrix
                        rr <- r2m(focal(px2r(c.px, im.dim = c(2048, 2048)), k))
                        
                        # find candidate lines
                        cand <- data.frame(which(rr >= 2 * kernel.size, arr.ind = T))
                        cc <- ddply(cand, .(x = row), summarise,
                                    size = length(row), len = max(col) - min(col) + 1, cover = size / len)
                        cc <- cc[cc$size > min.length, ]
                        
                        # for each candidate line, set pixels on that columns to feature type "line"
                        if (is.na(midline)) {
                            r.rng <- 1:max(cand$col)
                        } else {
                            if(max(cand$col) > midline) {
                                r.rng <- ceiling(midline):max(cand$col)     # rows in upper panel
                            } else {
                                r.rng <- 1:floor(midline)                   # rows in lower panel
                            }
                        }
                        data.frame(x = rep(cc$x, each = length(r.rng)), 
                                   r.rng = rep(r.rng, length(cc$x)),
                                   col.type = "line.c")
                    })
    lines <- rbind.fill(lines[sapply(lines, nrow) > 0])
    px.new <- merge(px, lines, by = c(1:2), all.x = T)
    px.new$f.type[px.new$col.type == "line.c"] <- "line.c"
    
    # return only original columns, with new "f.type" column if not already present
    return(px.new[, colnames(px.new) %in% c(colnames(px), "f.type")])
}


#' Find row features in pixel map
#' 
#' Applies a linear convolution filter to identify horizontal linear features in the pixel map.
#' @param px Data frame containing pixel coordinates to be grouped
#' @param kernel.size Edge length of square kernel to be used in convolution. Default is 5 (aims to pick up line segments of length 5 or greater) 
#' @param min.length Minimum length of segment to retain after convolution. Default is 3 * kernel size.
#' @return Data frame px, with column 'f.type' labelling columns found as feature type 'line.c'
#' @export
#' 
find.rows <- function(px, kernel.size = 5, min.length = 3 * kernel.size) {
    
    # if kernel size is even sides, add 1
    if(kernel.size %% 2 == 0) {kernel.size <- kernel.size + 1}
    
    # define kernel to use in convolution
    k <- t(matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
                    rep(kernel.size - 1, kernel.size), 
                    rep(-1, kernel.size * floor(kernel.size / 2))),
                  ncol = kernel.size))
    
    # convert px to raster, convolve with kernel, convert back to matrix
    rr <- r2m(focal(px2r(px, im.dim = apply(px[,1:2], 2, max) + kernel.size * c(2,2)), k))
    
    # find candidate lines
    cand <- data.frame(which(rr >= 2 * kernel.size, arr.ind = T))
    cc <- ddply(cand, .(x = col), summarise,
                size = length(row), len = max(row) - min(row) + 1, cover = size / len)
    cc <- cc[cc$size > min.length, ]
    
    c.rng <- 1:max(px$row)
    
    rows.found <- data.frame(c.rng = rep(c.rng, length(cc$x)),
                             x = rep(cc$x, each = length(c.rng)), 
                             row.type = "line.r")
    
    px.new <- merge(px, rows.found, by = c(1:2), all.x = T)
    px.new$f.type[px.new$row.type == "line.r"] <- "line.r"
    
    # return only original columns, with new "f.type" column if not already present
    return(px.new[,colnames(px.new) %in% c(colnames(px), "f.type")])
}


#' Find dense regions in pixel map
#' 
#' Apply a density filter to the pixel map and identify the densest regions.
#' @param px Data frame containing pixel coordinates to be examined
#' @param th.u Upper density threshold: any region with density above this threshold will be marked. Default is 0.5
#' @param th.l Lower density threshold: any region with density above this threshold, if contiguous with density above the upper threshold, will be marked. Ensures that the edges of dense regions are also identified. Default is 0.05
#' @param area Size of area around each pixel to consider when calculating the density at that pixel. Default is 11 (square extending 1mm out from centre pixel in each direction)
#' @param dilate.by Size of structuring element to use in smoothing and dilating the edges of the identified area, to ensure that all edge pixels are included. Default is 21.
#' @return Original pixel list, with column 'f.type' populated with factor 'dense.region' where appropriate.
#' @export
#' 
dense.regions <- function(px, th.u = 0.5, th.l = 0.05, area = 11, dilate.by = 21) {
    
    # filter out any features already identified
    if("f.type" %in% colnames(px)) {
        fpx <- px[is.na(px$f.type),]
    } else {
        fpx <- px
    }
    
    # define kernel
    k <-  matrix(rep(1 / area^2, area^2), ncol = area)
    
    # convolve image with kernel
    px.density <- r2m(focal(px2r(fpx, im.dim = (apply(px[,1:2], 2, max) + c(2, 2) * area)), k))
    th.density <- (px.density > th.l) * 1
    
    # get clumps of higher-density pixels
    cc <- clump(m2r(th.density))
    
    cand <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                       id = getValues(cc)[!is.na(getValues(cc))])
    cand$density <- px.density[as.matrix(cand[,1:2])]
    
    blocks <- ddply(cand, .(id), max.d = max(density), summarise)
    
    if ((dilate.by == 0) | is.na(dilate.by)) {
        qq <- as.matrix(cand[cand$id %in% blocks$id[blocks$max.d >= th.u], c("x", "y")])
    } else {
        th.density[as.matrix(cand[cand$id %in% blocks$id[blocks$max.d < th.u], c("x", "y")])] <- 0
        th.density[is.na(th.density)] <- 0
        
        sk <- shapeKernel(c(dilate.by, dilate.by), type = "disc")
        zz <- dilate(th.density, sk)
        
        qq <- which(zz == 1, arr.ind = T)
    }
    
    # label all free pixels within affected regions
    npx <- merge(px, data.frame(qq, dense = T), by = c(1:2), all.x = T)
    npx$f.type[npx$dense & is.na(npx$f.type)] <- "dense.region"
    
    return(npx[,colnames(npx) %in% c(colnames(px), "f.type")])
}


#' Find and label clusters of pixels
#' 
#' Cluster pixels and identify them, labelling outer pixels as cluster bodies and the brightest central pixel as the root.
#' @param px Data frame containing pixel coordinates to be examined
#' @param midline Location of panel midline, which electrically separates the upper and lower panels, and so cannot be crossed by column defects. Default is 1024.5; if the detector has no midline separation, enter NA.
#' @return Original pixel list, with column 'f.type' populated with factor 'cl.body' or 'cl.root' where appropriate.
#' @export
#'
find.clusters <- function(px, midline = 1024.5) {
    
    # filter out any features already identified
    if("f.type" %in% colnames(px)) {
        fpx <- px[is.na(px$f.type),]
    } else {
        fpx <- px
    }
    
    # convert free px to raster & get clumps of pixels
    cc <- clump(px2r(fpx, im.dim = apply(px[,1:2], 2, max)), dir = 4)
    
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    # filter out single pixels, label larger groups as clusters
    xy <- merge(xy, count(xy, "id"), all = T)[,c("x", "y", "id", "freq")]
    xy <- xy[xy$freq > 1,]
    xy$cl.type <- "cl.body"
    
    # ---------------------------------------------------------------------------
    # FIND CLUSTER 'ROOTS'
    
    xy <- merge(xy, px, by = c(1:2), all.x = T)
    
    # iterate over clusters to identify most severe defect
    roots <- lapply(unique(xy$id),
                    function(i) {
                        cand <- xy[xy$id == i,]
                        
                        # find most severe defect type in cluster
                        cand <- cand[cand$type == min(cand$type),]
                        
                        # if multiple candidates, get closest to horizontal midpoint
                        if (nrow(cand) > 1) {
                            
                            mp <- mean(cand$x)
                            
                            # if midpoint is between two pixels, shift left
                            mp <- mp - ((mp - as.integer(mp) != 0) * 0.5)
                            
                            cand$xd <- abs(cand$x - mp)
                            cand <- cand[cand$xd == min(cand$xd),]
                            
                            # still multiple candidates? Find pixel closest to readout sensor
                            if (nrow(cand) > 1) {
                                
                                # identify edge closest to readout sensor for this cluster
                                if (is.na(midline)) {
                                    edge <- max(px$col)
                                } else {
                                    if (max(cand$y) < midline) {
                                        edge <- 0
                                    } else {
                                        edge <- max(px$col)
                                    }
                                }
                                cand$yd <- abs(cand$y - edge)
                                cand <- cand[cand$yd == min(cand$yd),]
                            }
                        }
                        data.frame(x = cand$x, y = cand$y, cl.type = "cl.root")
                    })
    
    npx <- merge(px, xy[,c("x", "y", "cl.type")], by = c(1:2), all.x = T)
    npx <- merge(npx, rbind.fill(roots), by = c(1:2), all.x = T)
    npx$f.type[npx$cl.type.x == "cl.body"] <- "cl.body"
    npx$f.type[npx$cl.type.y == "cl.root"] <- "cl.root"
    
    return(npx[colnames(npx) %in% c(colnames(px), "f.type")])
}

