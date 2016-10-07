
# SCREEN SPOTS                                                                                  ####

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
    linear <- fit.gv.lm(im)
    
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
    
    # filter to check if spots found in multiple lists
    sp <- spot.list[sapply(spot.list, nrow) > 0]
    
    # if no screen spots found, return original pixel list
    if(length(sp) == 0) {
        return(px)          # no pixels to label
    }
    
    spots <- unique(data.frame(do.call("rbind", sp), s.spot = T))
    
    # identify individual spots
    cc <- clump(px2r(spots), dir = 4)
    
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    xy <- merge(xy, px, by = c(1:2), all.x = T)
    
    # filter out any spots where no component pixels were identified as defects
    xy.s <- ddply(xy, .(id), summarise, defects = sum(!is.na(type)), chk = sum(is.na(type)), size = length(x))
    xy <- xy[xy$id %in% xy.s$id[xy.s$defects > 0],]
    
    if(nrow(xy) == 0) {
        return(px)          # no pixels to label
    }
    
    npx <- merge(px, data.frame(xy[,1:2], s.spot = T), by = c(1:2), all.x = T)
    
    npx$f.type[npx$s.spot & npx$type %in% c("dim", "l.dim", "nl.dim", "spot.dim", NA)] <- "s.spot"
    
    return(npx[colnames(npx) %in% c(colnames(px), "f.type")])
}

####################################################################################################

# MULTI-PIXEL DEFECTS                                                                           ####

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
find.columns <- function(px, kernel.size = 5, midline = 1024.5, min.length = 3 * kernel.size, max.sep = 100) {
    
    # filter out any features already identified
    if("f.type" %in% colnames(px)) {
        fpx <- px[is.na(px$f.type),]
    } else {
        fpx <- px
        px$f.type <- NA
    }
    
    # if kernel size is even sides, add 1
    if(kernel.size %% 2 == 0) {kernel.size <- kernel.size + 1}
    
    # define kernels to use in convolution
    # find single columns
    k <- matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
                  rep(kernel.size - 1, kernel.size), 
                  rep(-1, kernel.size * floor(kernel.size / 2))),
                nrow = kernel.size)
    
    # find double columns
    k3 <- matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
                   rep((kernel.size - 1) / 3, kernel.size * 3), 
                   rep(-1, kernel.size * floor(kernel.size / 2))),
                 nrow = kernel.size)
    
    # find triple columns
    k5 <- matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
                   rep((kernel.size - 1) / 5, kernel.size * 5), 
                   rep(-1, kernel.size * floor(kernel.size / 2))),
                 nrow = kernel.size)
    
    # if midline given, split data & apply twice (row/col labels are inverted!)
    if(!is.na(midline)) {
        px.list <- list(fpx[fpx$col > midline,],
                        fpx[fpx$col < midline,])
    } else {
        px.list <- list(fpx)
    }
    
    lines <- lapply(px.list, 
                    function(c.px) {
                        
                        # convert px to raster, convolve with kernel, convert back to matrix
                        rr <- r2m(focal(px2r(c.px, im.dim = c(2048, 2048)), k, pad = T, padValue = 0))
                        rr3 <- r2m(focal(px2r(c.px, im.dim = c(2048, 2048)), k3, pad = T, padValue = 0))
                        rr5 <- r2m(focal(px2r(c.px, im.dim = c(2048, 2048)), k5, pad = T, padValue = 0))
                        
                        # find candidate lines
                        cand <- data.frame(which(rr >= 2 * kernel.size | rr3 >= 2 * kernel.size |
                                                     rr5 >= 2 * kernel.size, arr.ind = T))
                        cc <- ddply(cand, .(x = row), summarise, lower = min(col), upper = max(col),
                                    size = length(row), len = max(col) - min(col) + 1, cover = size / len)
                        cc <- cc[cc$size > min.length, ]
                        
                        if(nrow(cc) > 0) {
                            
                            # set extension to each column
                            if (is.na(midline)) {
                                cc$end <- cc$upper + max.sep
                                cc$start <- cc$lower - max.sep
                            } else {
                                if(max(cand$col) > midline) {
                                    cc$end <- sapply(cc$lower, function(lim) max(lim - 100, ceiling(midline)))
                                    cc$start <- cc$upper + max.sep
                                } else {
                                    cc$end <- min(cc$lower - max.sep, 0)
                                    cc$start <- sapply(cc$upper, function(lim) min(lim + 100, floor(midline)))
                                }
                            }
                            
                            zz <- apply(cc, 1, 
                                        function(ll) {
                                            ll <- unlist(ll)
                                            data.frame(x = ll["x"], y = c(ll["end"]:ll["start"]), row.names = NULL)
                                        })

                            # for each candidate line, set pixels on that column to feature type "line"
                            if(class(zz) == "list") {
                                data.frame(do.call("rbind", zz),
                                                   col.type = "line.c")
                            } else {
                                data.frame(zz, col.type = "line.c")
                            }
                        }
                    })
    
    lines <- lines[!sapply(sapply(lines, nrow), is.null)]
    
    if(length(lines) == 0) {return(px)}
    
    if(length(lines) == 2) {lines <- rbind.fill(lines[sapply(lines, nrow) > 0])}
    
    px.new <- merge(px, lines, by = c(1:2), all.x = T)
    px.new$f.type[px.new$col.type == "line.c" & is.na(px.new$f.type)] <- "line.c"
    
    # return only original columns, with new "f.type" column if not already present
    return(px.new[, colnames(px.new) %in% c(colnames(px), "f.type")])
}


#' Summarise columns by width
#' 
#' Group all adjacent column defects in pixel map and summarise each column defect by width
#' @param px Data frame containing pixel coordinates to be grouped
#' @param midline Location of panel midline, which electrically separates the upper and lower panels, and so cannot be crossed by column defects. Default is 1024.5; if the detector has no midline separation, enter NA.
#' @return data.frame summarising multi-column features identified by panel, location and width
#' @export
#' 
count.columns <- function(px, midline = 1024.5) {
    
    fpx <- px[px$f.type == "line.c",]
    
    if (nrow(fpx) == 0) {return(NULL)}
    
    if (is.na(midline)) {
        px.list <- list(fpx)
    } else {
        px.list <- list(upper = fpx[fpx$col > midline,],
                        lower = fpx[fpx$col < midline,])
        px.list <- px.list[sapply(px.list, nrow) > 0]
    }
    
    lines <- lapply(px.list, 
                    function(c.px) {
                        
                        # fill in broken lines
                        c.fill <- ddply(c.px, .(row), summarise, start = min(col), end = max(col))
                        cf <- apply(c.fill, 1, 
                                    function(fcol) {
                                        fcol <- unlist(fcol)
                                        data.frame(x = fcol["row"], 
                                                   y = c(fcol["start"] : fcol["end"]), row.names = NULL)
                                    })
                        cf <- do.call("rbind", cf)
                        
                        # clump adjacent pixels together
                        cc <- clump(px2r(cf))
                        xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                                         id = getValues(cc)[!is.na(getValues(cc))])
                        
                        ll <- ddply(xy, .(id), summarise,
                                    width = max(x)-min(x)+1,
                                    start = min(x), range = max(x))
                        ll$range[ll$width > 1] <- apply(ll[ll$width > 1,3:4], 1, paste, collapse = ":")
                        ll[order(ll$start),c("range", "width")]
                    })
    
    # label upper/lower panel, remove any tables with no columns found 
    if (length(lines) > 1) {
        lines <- lapply(names(lines[sapply(lines, nrow) > 0]), 
                        function(panel) data.frame(panel, lines[[panel]], stringsAsFactors = F))
    }
    
    lines <- rbind.fill(lines)
    return(lines)
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
    
    # filter out any features already identified
    if("f.type" %in% colnames(px)) {
        fpx <- px[is.na(px$f.type),]
    } else {
        fpx <- px
        px$f.type <- NA
    }
    
    # if kernel size is even sides, add 1
    if(kernel.size %% 2 == 0) {kernel.size <- kernel.size + 1}
    
    # define kernel to use in convolution
    k <- t(matrix(c(rep(-1, kernel.size * floor(kernel.size / 2)), 
                    rep(kernel.size - 1, kernel.size), 
                    rep(-1, kernel.size * floor(kernel.size / 2))),
                  ncol = kernel.size))
    
    # convert px to raster, convolve with kernel, convert back to matrix
    rr <- r2m(focal(px2r(fpx, im.dim = apply(px[,1:2], 2, max) + kernel.size * c(2,2)), k))
    
    # find candidate lines
    cand <- data.frame(which(rr >= 2 * kernel.size, arr.ind = T))
    cc <- ddply(cand, .(x = col), summarise,
                size = length(row), len = max(row) - min(row) + 1, cover = size / len)
    cc <- cc[cc$size > min.length, ]
    
    if (nrow(cc) == 0) {
        return(px)
    }
    
    c.rng <- 1:max(px$row)
    
    rows.found <- data.frame(c.rng = rep(c.rng, length(cc$x)),
                             x = rep(cc$x, each = length(c.rng)), 
                             row.type = "line.r")
    
    px.new <- merge(px, rows.found, by = c(1:2), all.x = T)
    px.new$f.type[px.new$row.type == "line.r" & is.na(px$f.type)] <- "line.r"
    
    # return only original columns, with new "f.type" column if not already present
    return(px.new[,colnames(px.new) %in% c(colnames(px), "f.type")])
}


#' Find dense regions in pixel map
#' 
#' Apply a density filter to the pixel map and identify the densest regions.
#' @param px Data frame containing pixel coordinates to be examined
#' @param offset Vector giving x- and y- offsets of image. To be used in padding the image to avoid edge effects
#' @param th.u Upper density threshold: any region with density above this threshold will be marked. Default is 0.5
#' @param th.l Lower density threshold: any region with density above this threshold, if contiguous with density above the upper threshold, will be marked. Ensures that the edges of dense regions are also identified. Default is 0.05
#' @param area Size of area around each pixel to consider when calculating the density at that pixel. Default is 11 (square extending 1mm out from centre pixel in each direction)
#' @param dilate.by Size of structuring element to use in smoothing and dilating the edges of the identified area, to ensure that all edge pixels are included. Default is 21.
#' @return Original pixel list, with column 'f.type' populated with factor 'dense.region' where appropriate.
#' @export
#' 
dense.regions <- function(px, th.u = 0.5, area = 11, th.l = 2/area, dilate.by = 21) {
    
    # filter out any features already identified
    if("f.type" %in% colnames(px)) {
        fpx <- px[is.na(px$f.type),]
    } else {
        fpx <- px
        px$f.type <- NA
    }
    
    # define kernel
    k <-  matrix(rep(1 / area^2, area^2), ncol = area)
    
    # convert px to matrix & trim offset regions
    padded <- array(0, dim = apply(px[,1:2], 2, max))
    padded[as.matrix(fpx[,1:2])] <- 1
    
    px.lim <- apply(which(padded > 0, arr.ind = T), 2, min)
    padded <- padded[px.lim[1]:nrow(padded), px.lim[2]:ncol(padded)]
    
    # convolve image with kernel
    px.density <- r2m(focal(m2r(padded), k, na.rm = T))
    th.density <- (px.density > th.l) * 1
    
    # get clumps of mid-density pixels
    cc <- clump(m2r(th.density))
    
    cand <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                       id = getValues(cc)[!is.na(getValues(cc))])
    cand$density <- px.density[as.matrix(cand[,1:2])]
    
    if (nrow(cand) == 0) {return(px)}
    
    blocks <- ddply(cand, .(id), summarise, max.d = max(density))
    
    if ((dilate.by == 0) | is.na(dilate.by)) {
        qq <- as.matrix(cand[cand$id %in% blocks$id[blocks$max.d >= th.u], c("x", "y")])
    } else {
        th.density[as.matrix(cand[cand$id %in% blocks$id[blocks$max.d < th.u], c("x", "y")])] <- 0
        th.density[is.na(th.density)] <- 0
        
        sk <- shapeKernel(c(dilate.by, dilate.by), type = "disc")
        zz <- dilate(th.density, sk)
        
        qq <- which(zz == 1, arr.ind = T)
    }
    
    # adjust coordinates for trimmed regions
    qq[,1] <- qq[,1] + px.lim[1] - 1
    qq[,2] <- qq[,2] + px.lim[2] - 1
    
    if(nrow(qq) == 0) {return(px)}
    
    # label all free pixels within affected regions
    npx <- merge(px, data.frame(qq, dense = T), by = c(1:2), all.x = T)
    npx$f.type[npx$dense & is.na(npx$f.type)] <- "dense.region"
    
    return(npx[, colnames(npx) %in% colnames(px)])
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


#' Count clusters by size
#' 
#' Count all clusters of given size in pixel map. Also checks main direction of spread (horizontal or vertical) 
#' @param px Data frame containing pixel coordinates to be examined
#' @param check.dir Boolean: instead of producing summary table, should paired t-test of mean dimensions of clusters be performed? Default is F.
#' @return Data.frame summarising sizes of clusters found by dominant direction of spread
#' @export
#' 
count.clusters <- function(px, check.dir = F) {
    c.px <- px[px$f.type %in% c("cl.root", "cl.body"),]
    
    cc <- clump(px2r(c.px), dir = 4)
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    cl <- ddply(xy, .(id), summarise,
                xm = mean(x), ym = mean(y), size = length(x), 
                height = max(y) - min(y) + 1, width = max(x) - min(x) + 1,
                diff = height - width,
                type = c("single", "double", "small", "medium", "large", "mega", "error")[findInterval(size, c(0, 1.5, 2.5, 3.5, 9.5, 17.5, 36.5, 99999))],
                dir = c("flat", "round", "tall")[findInterval(diff, c(-3000, -1.5, 1.5, 3000))])
    
    if (check.dir) {
        return(t.test(cl$width, cl$height, paired = T))
    }

    cl$dir[cl$size <= 4 & cl$diff >= 1] <- "tall"
    cl$dir[cl$size <= 4 & cl$diff <= -1] <- "flat"
    
    clc <- ddply(cl, .(type, dir), summarise, freq = length(type), o = min(size))
    clc[order(clc$o),1:3]
}

####################################################################################################

