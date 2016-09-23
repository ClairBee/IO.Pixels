
#' Basic bad pixel map
#' 
#' @export
#' 
basic.bpx <- function(acq, screen.spots = T) {
    
    # calculate thresholds
    th <- apply(acq[,,c("black", "grey")], 3, function(im) {
        med <- median(im, na.rm = T)
        c(v.dim = med * 0.5, dim = med * 0.75,
          bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
    })
    
    px <- list(edge = edge.px(acq, edge.width = 40),
               no.resp = no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]),
               hot = which(acq[,,"black"] == 65535, arr.ind = T),
               dead = which(acq[,,"grey"] == 0, arr.ind = T),
               line = which(find.lines(acq[, , "black"], midline = 1024.5) > 0, arr.ind = T),
               v.bright = which(acq[, , "black"] > th["v.bright", "black"] | 
                                    acq[, , "grey"] > th["v.bright", "grey"], arr.ind = T),
               bright = which(acq[, , "black"] > th["bright", "black"] |
                                  acq[, , "grey"] > th["bright", "grey"], arr.ind = T),
               v.dim = which(acq[, , "black"] < th["v.dim", "black"] | 
                                 acq[, , "grey"] < th["v.dim", "grey"], arr.ind = T), 
               dim = which(acq[, , "black"] < th["dim", "black"] | 
                               acq[, , "grey"] < th["dim", "grey"], arr.ind = T))
    
    if (screen.spots) {
        px$screen.spot <- screen.spots(acq[,,"white"], enlarge = T, ignore.edges = 40)
    }
    
    # remove empty types & combine into data frame
    px <- px[sapply(px, length) > 0]
    px <- rbind.fill(lapply(names(px), function(df) data.frame(px[[df]], type = toString(df))))
    
    Cat <- c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "edge", "line.d", "screen.spot", "v.dim", "dim")
    px$type <- ordered(px$type, levels = Cat)
    
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    return(px)
}



#' Linearity-based pixel map
#' 
#' @export
#' 
nl.bpx <- function(acq, w.th = 1000) {
    
    # calculate thresholds
    th <- apply(acq[,,c("black", "grey")], 3, function(im) {
        med <- median(im, na.rm = T)
        c(v.dim = med * 0.5, dim = med * 0.75,
          bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
    })
    
    wlm <- fit.w.lm(acq)$df
    
    px <- rbind(data.frame(which(acq[,,"black"] > th["bright", "black"], arr.ind = T), type = "warm"),
                data.frame(which(acq[,,"grey"] > th["bright", "grey"], arr.ind = T), type = "bright"),
                data.frame(setNames(wlm[abs(wlm$res) > w.th, 1:2], nm = c("row", "col")), type = "nonlinear"),
                data.frame(which(acq[,,"grey"] - acq[,,"black"] < 500 & acq[,,"grey"] > acq[,,"black"], arr.ind = T), type = "stuck"),
                data.frame(which(acq[,,"black"] < th["dim", "black"], arr.ind = T), type = "cool"),
                data.frame(which(acq[,,"grey"] < th["dim", "grey"], arr.ind = T), type = "dim"))
    
    px$type <- ordered(px$type, levels = c("stuck", "warm", "bright", "nonlinear", "cool", "dim"))
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    
    return(px)
}


####################################################################################################


#' Identify 'bad' pixels using quantiles of Lowess residuals
#'
#' Load a specified image and apply Lowess smoothing along columns, then cut residuals at specified quantiles to classify 'bad' pixels
#' @details Requires that pixelwise mean files are already loaded, as by \link{\code{load.pixel.maps}}
#' @param img.date Integer/string specifying image date to be imported
#' @param batch String specifying colour batch to be imported (black, grey or white)
#' @param span Smoothing span to be used for Lowess. Default is 1/15
#' @param lq Lower quartile at which to start looking for 'outliers'. Default is 0.001
#' @param uq Upper quartile at which to start looking for 'outliers'. Default is 0.999
#' @param dist Multiples of IQR away from lq/uq at which a point is classified as 'bad'.
#' @param details Boolean: return lists of identified 'bad pixels' (T) or only summary (F)? Default is F.
#' @export
#' @examples
#' load.pixel.maps()
#' zz <- rbind(bp.lowess.quantiles(160314, "white"),
#'             bp.lowess.quantiles(151015, "white"),
#'             bp.lowess.quantiles(150828, "white"))
#' 
bp.lowess.quantiles <- function(img.date, batch, span = 1/15, lq = 0.001, uq = 0.999, dist = 1.5, details = F) {
    img.date <- toString(img.date)
    b <- tolower(substring(batch,1,1))
    
    im <- get(paste0("pw.", b))[,,img.date]
    
    smoo <- lowess.per.column(im, span = span)
    res <- im - smoo
    
    ul <- quantile(res, uq) + (dist * IQR(res))
    ll <- quantile(res, lq) - (dist * IQR(res))
    
    # convert parameters to strings for summary export
    if (missing(lq)) {lq <- "0.001"} else {lq <- format(as.list(match.call())$lq, scientific = F)}
    if (missing(uq)) {uq <- "0.999"} else {uq <- format(as.list(match.call())$uq, scientific = F)}
    if (missing(span)) {f <- "1/15"} else {f <- format(as.list(match.call())$span, scientific = F)}
    
    summ <- data.frame(method = "lowess-quantiles", img.date, batch, 
                       params = paste(f, paste(lq, uq, dist, sep = ", "), sep = "; "),
                       lower = ll, upper = ul,
                       n.low = length(which(res < ll)),
                       n.high = length(which(res > ul)),
                       res.mean = mean(res),
                       res.sd = sd(res),
                       res.mad = mad(res),
                       res.med = median(res),
                       row.names = NULL)
    
    if (details) {
        return(list(summ,
                    low = which(res < ll, arr.ind = T),
                    high = which(res > ul, arr.ind = T)))
    } else {
        return(summ)
    }
}



#' Identify 'bad' pixels using quantiles of Lowess residuals
#'
#' Load a specified image and apply Lowess smoothing along columns, then fit a Johnson model and use its quantiles to extract residuals
#' @details Requires that pixelwise mean files are already loaded, as by \link{\code{load.pixel.maps}}
#' @param img.date Integer/string specifying image date to be imported
#' @param batch String specifying colour batch to be imported (black, grey or white)
#' @param span Smoothing span to be used for Lowess. Default is 1/15
#' @param lq Lower quartile at which to start looking for 'outliers'. Default is 0.001
#' @param uq Upper quartile at which to start looking for 'outliers'. Default is 0.999
#' @param method String: fit Johnson distribution by quantiles or by matching moments? Default is "quant"
#' @param details Boolean: return lists of identified 'bad pixels' (T) or only summary (F)? Default is F.
#' @export
#' @examples
#' load.pixel.maps()
#' zz <- rbind.fill(bp.lowess.quantiles(160314, "white"),
#'                  bp.lowess.johnson(160314, "white"))
#'                  
bp.lowess.johnson <- function(img.date, batch, span = 1/15, lq = 0.001, uq = 0.999, method = "quant", details = F) {
    img.date <- toString(img.date)
    b <- tolower(substring(batch,1,1))
    
    im <- get(paste0("pw.", b))[,,img.date]
    
    smoo <- lowess.per.column(im, span = span)
    res <- im - smoo
    
    if (method == "quant") {
        parms <- JohnsonFit(c(res), moment = "quant")
    } else {
        parms <- JohnsonFit(c(res), moment = "find")
    }
    
    ll <- qJohnson(lq, parms)
    ul <- qJohnson(uq, parms)
    
    # convert parameters to strings for summary export
    if (missing(lq)) {lq <- "0.001"} else {lq <- format(as.list(match.call())$lq, scientific = F)}
    if (missing(uq)) {uq <- "0.999"} else {uq <- format(as.list(match.call())$uq, scientific = F)}
    if (missing(span)) {f <- "1/15"} else {f <- format(as.list(match.call())$span, scientific = F)}
    
    summ <- data.frame(method = "lowess-Johnson", img.date, batch, 
                       params = paste(f, paste(lq, uq, sep = ", "), method, sep = "; "),
                       lower = ll, upper = ul,
                       n.low = length(which(res < ll)),
                       n.high = length(which(res > ul)),
                       res.mean = mean(res),
                       res.sd = sd(res),
                       res.mad = mad(res),
                       res.med = median(res),
                       fitted.params = paste(lapply(c(parms)[1:4], round, 3), c(parms)[5], collapse = ", "),
                       row.names = NULL)
    
    if (details) {
        return(list(summ,
                    low = which(res < ll, arr.ind = T),
                    high = which(res > ul, arr.ind = T), 
                    model = parms))
    } else {
        return(summ)
    }
}




#' Identify 'bad' pixels using quantiles of residuals after parametric fitting
#'
#' Load a specified image and fit parametric circular spot and panelwise gradients, then cut residuals at specified quantiles to classify 'bad' pixels
#' @details Requires that pixelwise mean files are already loaded, as by \link{\code{load.pixel.means}}
#' @param img.date Integer/string specifying image date to be imported
#' @param batch String specifying colour batch to be imported (black, grey or white)
#' @param spot.order Integer (or NA if no spot model is required) specifying order of polynomial model to fit to circular spot.
#' @param panel.terms String specifying formula to apply in linear regression over panel coordinates.
#' @param lq Lower quartile at which to start looking for 'outliers'. Default is 0.001
#' @param uq Upper quartile at which to start looking for 'outliers'. Default is 0.999
#' @param dist Multiples of IQR away from lq/uq at which a point is classified as 'bad'.
#' @param robust Boolean: use robust model fitting (\link{\code{rlm}})? Default is F (use \link{\code{lm}} for regression).
#' @param details Boolean: return lists of identified 'bad pixels' (T) or only summary (F)? Default is F.
#' @export
#' @examples
#' load.pixel.means()
#' zz <- rbind(bp.parametric.quantiles(160314, "white", ),
#'             bp.lowess.quantiles(160314, "white"))
#' 
bp.parametric.quantiles <- function(img.date, batch, 
                                    spot.order = 2, panel.terms = "x + y", 
                                    lq = 0.001, uq = 0.999, dist = 1.5, 
                                    robust = F, details = F) {
    img.date <- toString(img.date)
    im <- pw.m[,, batch, img.date]
    
    # fit circular model (if specified)
    if (!is.na(spot.order)) {
        circ.lm <- spot.lm(im, order = spot.order, robust = robust)
        c.par <- paste0(toString(spot.order), "-spot")
        circ.res <- matrix(circ.lm$residuals, ncol = 1996)
    } else {
        circ.res <- im
        c.par = "no spot"
    }
    
    # fit per-panel model (if specified)
    if (!is.na(panel.terms)) {
        panel.lm <- panel.lm(circ.res, panel.terms, robust)
        p.par <- paste0("panel ", panel.terms)
        res <- circ.res - panel.lm$fitted.values
    } else {
        res <- circ.res
        p.par <- "no panels"
    }
    
    # cutoff points for 'bad' pixels, based on residual quantiles
    ul <- quantile(res, uq) + (dist * IQR(res))
    ll <- quantile(res, lq) - (dist * IQR(res))
    
    healthy <- res[which(res >= ll & res <= ul, arr.ind = T)]
    
    summ <- data.frame(method = paste(c.par, p.par, "quantiles", sep = "; "),
                       img.date, batch,
                       params = paste(lq, uq, dist, sep = ", "), 
                       lower = ll, upper = ul, 
                       n.low = length(which(res < ll)),
                       n.high = length(which(res > ul)),
                       res.mean = mean(res), 
                       res.sd = sd(res), res.mad = mad(res), res.med = median(res), 
                       healthy.mean = mean(healthy), healthy.sd = sd(healthy), 
                       healthy.mad = mad(healthy), healthy.med = median(healthy),
                       row.names = NULL)
    
    if (details) {
        return(list(details = summ, 
                    low = which(res < ll, arr.ind = T), 
                    high = which(res > ul, arr.ind = T)))
    } else {
        return(summ)
    }
}


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
#' @param limit Quantile of black pixel population to be used as cutoff for 'normal' behaviour. Default is 0.0005, capturing central 99.9percent of distribution.
#' @return Matrix of coordinates of non-responsive pixels
#' @export
#' @examples
#' qq <- no.response(160314, limit = 0.001)
no.response <- function(bright.image, dark.image, limit = 0.0005, exclude.edge = 10) {
    
    nc <- ncol(bright.image)
    
    bn <- qJohnson(c(limit, 1-limit), JohnsonFit(dark.image[!is.na(dark.image)]))
    
    un <- rbind(which(matrix(findInterval(bright.image, bn), ncol = nc) == 1, arr.ind = T))
    
    # filter out any pixels that fall in excluded edge region
    un[un[,1] <= exclude.edge, ] <- NA
    un[un[,1] >= ncol(bright.image) - exclude.edge, ] <- NA
    
    un[un[,2] <= exclude.edge, ] <- NA
    un[un[,2] >= nrow(bright.image) - exclude.edge,] <- NA
    
    return(un[!is.na(un[,1]),])
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



#' Identify clusters of pixels
#' 
#' Identify clusters of adjacent pixels of each type in a bad pixel map.
#' @param bpm Bad pixel map: a matrix or data frame containing coordinates of bad pixels and column "type" defining bad pixel type ("hot", "bright", "dim", "dead").
#' @return Data frame containing bad pixel coordinates and type, along with ID and size of cluster to which each is assigned.
#' @export
#' @examples
#' zz <- cluster.by.type(bpm)
#' 
cluster.by.type <- function(bpm) {
    
    # ensure that coordinates are unique
    bpm <- bpm[!duplicated(bpm[,1:2]),]
    
    # remove edge pixels from bad pixel map: redundant to cluster these
    px <- bpm[bpm$type != "edge",]
    
    # create matrix to populate raster
    px.vals <- array(dim = c(1996, 1996))
    px.vals[as.matrix(px[,1:2])] <- px$type
    
    # convert to raster
    r <- raster(t(px.vals[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5)
    
    # cluster bad pixels by type
    xy <- data.frame()
    for (val in unique(getValues(r)[!is.na(getValues(r))])) {
        
        # mask all but one type of bad pixel
        tmp <- r; 
        tmp[tmp != val] <- NA
        
        # clump remaining bad pixels together
        cc <- clump(tmp, dir = 4)
        
        xy <- rbind(xy,
                    merge(data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                                     type = val,
                                     id = getValues(cc)[!is.na(getValues(cc))]),
                          setNames(data.frame(val, freq(cc)), c("type", "id", "count"))))
    }
    xy$type <- factor(levels(px$type)[xy$type], levels = levels(px$type))
    
    # rationalize cluster IDs (otherwise IDs are reused for each type)
    xy <- transform(xy, id = match(apply(xy[, c("type", "id")], 1, paste, collapse = "-"),
                                   unique(apply(xy[, c("type", "id")], 1, paste, collapse = "-"))))
    
    # recombine with excluded edge pixels to retain full bad pixel map
    if (nrow(px) != nrow(bpm)) {
        xy <- rbind(xy,
                    setNames(data.frame(bpm[bpm$type == "edge",], 
                                        "id" = max(xy$count) + 1,
                                        "count" = nrow(bpm[bpm$type != "edge",])),
                             c("x", "y", "type", "id", "count")))
    }
    
    # return data frame of pixel coords, pixel type, cluster size & id
    return(xy[,c("x", "y", "type", "id", "count")])
}


#' Enlarge screen spots to include any adjacent dim pixels
#' 
#' Extend any screen spots in an image to include any adjacent dim pixels detected.
#' @param bpm Bad pixel map with categories assigned to each pair of coordinates, as generated by \link{\code{bad.px}}
#' @return Updated ad pxiel map
#' @export
#' 
enhance.spots <- function(bpm) {
    
    # filter to retain only slightly dim pixels and identified screen spots
    px <- bpm[bpm$type %in% c("s.dim", "screen spot"),]
    
    # identify clusters of each type of bad pixel
    org.clumps <- cluster.by.type(px)
    
    # now need to identify clusters containing either/both types of bad pixel
    
    # create matrix to populate raster
    px.vals <- array(dim = c(1996, 1996))
    px.vals[as.matrix(px[,1:2])] <- 1
    
    # convert to raster
    r <- raster(t(px.vals[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5)
    
    # clump pixels of both types together, get cluster sizes
    cc <- clump(r, dir = 4)
    new.clumps <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                             id = getValues(cc)[!is.na(getValues(cc))])
    new.clumps <- merge(new.clumps, count(new.clumps, "id"))
    
    matched <- merge(org.clumps[, c("x", "y", "type", "count")],
                     new.clumps[, c("x", "y", "freq", "id")])
    
    # assign a unique category to each new cluster identified
    # (if cluster contains screen spot, whole cluster becomes a screen spot)
    cat <- setNames(matched[order(matched$type), c("id", "type")], c("id", "new.type"))
    cat <- cat[!duplicated(cat$id),]
    
    # assign new category to each pixel
    matched <- merge(matched, cat)
    
    # recombine into complete bad pixel map
    if (nrow(px) != nrow(bpm)) {
        # rename columns in matched data to allow rbind
        df <- rbind(bpm[!bpm$type %in% c("s.dim", "screen spot"),],
                    setNames(matched[, c("x", "y", "new.type")], c("row", "col", "type")))
    }
    
    return(df)
}


#' Identify superclusters of pixels
#' 
#' Identify superclusters of bad pixels: small clusters of adjacent bad pixels of different types.
#' @param cl Bad pixel map or clustered bad pixel map: a matrix or data frame containing coordinates of bad pixels and column "type" defining bad pixel type, with cluster ID and count if already clustered by \link{\code{cluster.by.type}}
#' @param Vector of bad pixel types to exclude from superclusters. Default is to exclude screen spots and edge pixels.
#' @details Only clusters of size 9 or less will be considered as part of a supercluster, with the exception of bright or dim columns or rows identified. To exclude bright or dim lines, add those types to the 'excl' parameter. 
#' @return Data frame containing bad pixel coordinates and type, along with ID and size of clusters and superclusters to which each is assigned, and a shape classification based on these.
#' @export
#' @examples
#' # should both give same result
#' zz <- superclusters(cluster.by.type(bp))
#' zz <- superclusters(bp)
#' 
superclusters <- function(cl, excl = c("")) {
    
    # if data isn't already clustered, do so now
    if (is.null(cl$id) | is.null(cl$count)) {
        cl <- cluster.by.type(cl)
    }
    
    # exclude screen spots and any other categories to be excluded
    px <- cl[!(cl$type %in% c("screen spot", "edge", excl)),]
    
    # exclude clusters larger than 9 pixels: only interested in small clusters here
    # however, retain columns or rows of bright pixels (if not already excluded by user)
    # generally, removes areas of slightly dim pixels in corners/edges of image
    px <- px[(px$count <= 9) | (px$type %in% c("line.b", "line.d")),]
    
    # create matrix to populate raster
    px.vals <- array(dim = c(1996, 1996))
    px.vals[as.matrix(px[,1:2])] <- px$type
    
    # convert to raster
    r <- raster(t(px.vals[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5)
    
    # get clumps of pixels
    sc.cc <- clump(r, dir = 4)
    sc.xy <- merge(data.frame(xyFromCell(sc.cc, which(!is.na(getValues(sc.cc)))),
                              sc.id = getValues(sc.cc)[!is.na(getValues(sc.cc))]),
                   setNames(data.frame(freq(sc.cc)), c("sc.id", "sc.count")))
    
    # combine into single table of bad pixels with cluster & supercluster information
    xy <- merge(px, sc.xy, all = T)
    
    # assign supercluster classification to each pixel
    xy$sc.type <- factor("other", levels = c("singleton", "cluster", "supercluster", "other"))
    xy$sc.type[xy$count <= 9] <- "cluster"
    xy$sc.type[xy$count == 1] <- "singleton"
    xy$sc.type[xy$sc.count > xy$count] <- "supercluster"
    
    # recombine with previously excluded pixels to give complete bad pixel map
    if (nrow(cl) != nrow(px)) {
        xy <- rbind(xy,
                    data.frame(cl[(cl$type %in% c("screen spot", "edge", excl)) |
                                      ((cl$count > 9) & (!cl$type %in% c("line.b", "line.d"))),], 
                               sc.id = NA, sc.count = NA, sc.type = "other"))
    }
    
    # rationalise supercluster IDs
    xy[xy$sc.type != "supercluster", c("sc.id", "sc.count")] <- NA
    xy[xy$sc.type == "supercluster",] <- transform(xy[xy$sc.type == "supercluster",], 
                                                   sc.id = match(sc.id, unique(sc.id)))
    return(xy)
}


#' Classify bad pixels according to the objects to which they belong.
#' 
#' Given a bad pixel map, classify singletons, clusters and superclusters according to their composition.
#' @param bpm Bad pixel map
#' @return Updated bad pixel map
#' @export
#' 
classify.states <- function(bpm) {
    
    bpm <- merge(bpm, 
                 data.frame(sc.id = unique(bpm$sc.id[!is.na(bpm$sc.id)]),
                            sc.comp = sapply(split(levels(bpm$type)[(which(apply(as.matrix(table(bpm$sc.id, bpm$type)), 1, ">", 0))-1) %% 11 + 1],
                                                   which(apply(as.matrix(table(bpm$sc.id, bpm$type)), 1, ">", 0)) %/% 11),
                                             paste, collapse = ", ")),
                 all = T, sort = F)[, c(colnames(bpm), "sc.comp")]
    
    # single pixels & small clusters are marked with type; slightly bright or dim clusters/singles are ignored.
    bpm$cat[bpm$sc.type == "singleton"] <- apply(cbind("s", levels(bpm$type)[bpm$type]), 1, paste, collapse = "-")[bpm$sc.type == "singleton"]
    bpm$cat[bpm$sc.type == "cluster"] <- apply(cbind("c", levels(bpm$type)[bpm$type]), 1, paste, collapse = "-")[bpm$sc.type == "cluster"]
    bpm$cat[bpm$sc.type %in% c("singleton", "cluster") & (bpm$type %in% c("s.bright", "s.dim"))] <- "marginal"
    
    # screen spots, edge pixels and large areas of bright or dim pixels are marked
    bpm$cat[bpm$sc.type == "other"] <- "other"
    bpm$cat[bpm$type == "screen spot"] <- "screen spot"
    bpm$cat[bpm$type == "edge"] <- "edge"
    
    # superclusters, separated by severity of damage
    bpm$cat[grep("bright", bpm$sc.comp)] <- "sc-bright"
    bpm$cat[grep("hot", bpm$sc.comp)] <- "sc-hot"
    bpm$cat[grep("no response", bpm$sc.comp)] <- "sc-no-response"
    bpm$cat[grep("dead", bpm$sc.comp)] <- "sc-dead"
    
    bpm$cat <- ordered(bpm$cat, levels = unique(c("sc-dead", "sc-no-response", "sc-hot", "sc-bright",
                                                  "c-dead", "c-no response", "c-hot", "c-v.bright", "c-bright", "c-v.dim", "c-dim",
                                                  "s-dead", "s-no response", "s-hot", "s-v.bright", "s-bright", "s-v.dim", "s-dim",
                                                  "screen spot", "edge", "other", "marginal"), levels(bpm$cat)))
    
    bpm <- bpm[,-which(colnames(bpm) == "sc.comp")]
    return(bpm)
}



#' Create array of all pixelwise means
#'
#' Load all pixelwise mean summaries (already created and stored as .rds objects) into a single labelled array.
#' @details RDS objects for import can be downloaded from https://github.com/ClairBee/Pixels/tree/master/02_Objects/images
#' @param fpath Path to top level of stored .rds objects. Default is "./02_Objects/images"
#' @param acq.list Vector containing acquisitions to import. If left blank, will import every file in the target folder named "pwm-xxxxxx.rds".
#' @export
#' 
load.pixel.means <- function(fpath = "./02_Objects/images", acq.list) {
    
    # remove trailing blackslash, if necessary
    fpath <- gsub("/$", "", fpath)
    
    if (missing (acq.list)) {
        acq.list <- gsub(".rds", "", gsub("pwm-", "", list.files(fpath, pattern = "pwm-[a-z, A-Z, 0-9]+\\.rds$", full.names = F)))
    }
    abind(sapply(acq.list, function(nm) readRDS(paste0(fpath, "/pwm-", nm, ".rds")), 
                 simplify = F),
          along = 4)
}

#' Create array of all pixelwise SDs
#'
#' Load all pixelwise standard deviation matrices into a single labelled array.
#' @param acq.list Vector containing acquisitions to import. If left blank, will import every file in the target folder named "pwm-xxxxxx.rds".
#' @param fpath Path to top level of stored .rds objects. Default is "./02_Objects/images"
#' @export
#' 
load.pixel.sds <- function(acq.list, fpath = "./02_Objects/stdevs") {
    if (missing (acq.list)) {
        acq.list <- gsub(".rds", "", gsub("pwsd-", "", list.files(fpath, pattern = "pwsd-[a-z, A-Z, 0-9]+\\.rds$", full.names = F)))
    }
    abind(sapply(acq.list, function(nm) readRDS(paste0(fpath, "/pwsd-", nm, ".rds")), 
                 simplify = F),
          along = 4)
}



####################################################################################################

#' Load one day's acquisitions from one batch
#'
#' Import one day's images from a single batch into an array. If the function is called by the user, the data is loaded into an object named automatically by the function. If called by another function, will return an array.
#' @param img.date Date of images to import, in format yymmdd.
#' @param batch Specify batch to import: black, grey or white.
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @param x Width of image. Default is 1996.
#' @param y Height of image. Default is 1996.
#' @param z Number of images to import. Default is 20.
#' @export
#' @examples
#' load.images(150828, "black")   # will create array b.150828
#' 
#' 
load.images <- function(img.date, batch, fpath = "/home/clair/Documents/Pixels/Image-data/", x = 1996, y = 1996, z = 20) {
    
    img.date <- toString(img.date)
    obj.nm <- paste0(substring(batch,1,1), ".",img.date)
    global = F
    
    # was function called directly or within another function?
    if (environmentName(parent.frame()) == "R_GlobalEnv") {
        global = T
        
        # check if object already exists: don't overwrite
        if (exists(obj.nm)) {
            err <- paste("An object named",obj.nm,"already exists. If you want to reload the data, please delete it first.")
            stop(err)
        }            
    }
    
    # convert date to filenm format, set up array
    rev.date <- paste0(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2))
    m <- array(dim = c(x, y, z))
    
    # read each file in turn, assign values into array
    # DON'T iterate over file names - images will not be in the correct order
    for (i in 1:z) {
        filenm <- paste0(fpath,img.date,"/",batch,"/",batch,"_",rev.date,"_",i,".tif")
        tmp <- readTIFF(filenm, as.is = T)
        
        # if imported array has more than 2 dimensions, keep only the first
        if (length(dim(tmp)) == 2) {
            tmp <- tmp
        } else {
            tmp <- tmp[,,rep(1, length(dim(tmp)) - 2)]
        }
        
        # flip & transpose matrix to match orientation of original image
        m[,,i] <- t(tmp[nrow(tmp):1,,drop=FALSE])
    } 
    
    # if called from global evironment, automatically create & name object
    # otherwise, return object so that it can be assigned by calling function
    if (global) {
        assign(obj.nm, m, envir = .GlobalEnv)
    } else {
        return(m)
    }
}


#' Load one day's acquisitions from all batches
#'
#' Import one day's images from into  1996x1996x20x3 array.
#' @param img.date Date of images to import, in format yymmdd.
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @export
#' @examples
#' img.150828 <- load.daily(150828)
#' 
#' 
load.daily <- function(img.date, fpath = "/home/clair/Documents/Pixels/Image-data/") {
    m <- array(dim = c(1996, 1996, 20, 3),
               dimnames = list(NULL, NULL, NULL, c("black", "grey", "white")))
    
    m[,,,1] <- load.images(img.date, "black")
    m[,,,2] <- load.images(img.date, "grey")
    m[,,,3] <- load.images(img.date, "white")
    return(m)
}



#' Load xml profiles of one day's acquisitions from one batch
#'
#' Support function: import xml profiles of one day's images into a single data frame
#' @param img.date Date of images to import, in format yymmdd.
#' @param batch Specify batch to import: black, grey or white.
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing all xml profile data for the specified image set.
#' @export
#' @examples
#' b.150828.xml <- load.profiles(150828, "black")
#' 
#' 
load.profiles <- function(img.date, batch, fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    img.date <- toString(img.date)
    
    fpath <- paste(fpath, img.date, "/", batch, "/", sep = "")
    files <- list.files(fpath, pattern = "\\.xml$")
    
    xml.data <- list()
    
    # read through all xml data and import all xml fields into single df per image
    for (i in 1:length(files)) {
        xml.data[[i]] <- as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(paste(fpath, files[i], sep = ""))),
                                                             list(recursive=TRUE))), FUN = unlist),
                                       stringsAsFactors = F)
    }
    
    # merge all list elements into single df, filling in missing columns with NA
    xml.data <- rbind.fill(xml.data)
    
    # where possible, convert character strings to numeric
    for (i in 1:ncol(xml.data)) {
        if (suppressWarnings(all(!is.na(as.numeric(xml.data[,i]))))) {
            
            xml.data[,i] <- as.numeric(xml.data[,i])
            
            # if column contains only whole numbers, convert to integer rather than numeric
            if (all(xml.data[,i] == floor(xml.data[,i]))) {
                xml.data[,i] <- as.integer(xml.data[,i]) 
            }
        }
    }
    
    # sort data frame in chronological order
    xml.data <- xml.data[order(strptime(xml.data$Time, "%H:%M:%S")),]
    return(xml.data)
}



#' Summarise image profile data
#'
#' Searches through folder structure, extracts all .xml profiles, and summarises into a data frame.
#' @details If NA appears in data frame, the acquisitions for that batch do not all share the same value of that attribute, and should be investigated individually.
#' @param fpath Path to top level of image directory. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing summaries of image profile data.
#' @export
#' @examples
#' df <- summarise.profiles()
#' 
#' 
summarise.profiles <- function(fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    # get list of folders in given path location
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # check that valid path has been given
    if (length(folders) == 0) {
        cat("Image folder",fpath,"not found in",getwd())
    } else {
        
        d <- length(folders)
        xml.summ <- list()
        
        # read each folder in turn
        for (i in 1:d) {
            
            batches <- tolower(list.dirs(paste(fpath, folders[i],"/", sep = ""), full.names = F, recursive = F))
            
            # read each batch in turn
            for (j in 1:length(batches)) {
                
                # extract data frame of xml data (may not be the same format for each day!)
                p <- load.profiles(folders[i], batches[j], fpath)
                
                # for numeric columns, check that min == max; otherwise, replace with NA
                for (k in c(1:ncol(p))[sapply(p, class) == "numeric"]) {
                    if (min(p[,k]) != max(p[,k])) {
                        p[,k] <- NA
                    }
                }
                xml.summ[[length(xml.summ) + 1]] <- cbind(date = folders[i],
                                                          batch = batches[j],
                                                          frames = nrow(p),
                                                          p[1,])
            }
        }
        df <- rbind.fill(xml.summ)
        df$date <- as.character(df$date)
        df$batch <- as.character(df$batch)
        df$frames <- as.integer(df$frames)
        return(df)
    }
}



#' Calculate summary statistics for all images
#'
#' Searches through folder structure, extracts all .tif images, and stores the summary statistics for each in a data frame.
#' @param fpath Path to top level of image directory. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing summaries of pixel values: mean and standard deviation, max, min, quartiles and IQR.
#' @export
#' @examples
#' df <- summarise.images()
#' 
#' 
summarise.images <- function(fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    df <- data.frame(img.date = character(),
                     batch = character(),
                     
                     mean = numeric(),
                     sd = numeric(),
                     
                     min = double(),
                     lq = double(),
                     median = double(),
                     uq = double(),
                     max = double(),
                     iqr = double(),
                     stringsAsFactors = F)                    
    
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # read each folder in turn
    for (i in 1:length(folders)) {
        
        batches <- tolower(list.dirs(paste(fpath, folders[i],"/", sep = ""), full.names = F, recursive = F))
        
        # read each batch in turn
        for (j in 1:length(batches)) {
            
            r <- nrow(df) + 1
            data <- load.images(folders[i], batches[j], fpath)
            summ <- batch.summary(data)
            
            df[r, c(1:2)] <- c(folders[i], batches[j])
            df[r, c(3:6)] <- summ[[1]]
            df[r, c(7:12)] <- summ[[2]]
        }
    }
    return(df)
}



#' Count available images
#'
#' Searches through folder structure to quickly summarise number of available .tif images from each batch on each day. Assumes files are organised as /Image-data/date/batch/image.tif
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing counts of .tif images in batch subfolders for each date, with image acquisition dates as row names.
#' @export
#' @examples
#' image.count <- count.images()
#' 
#' 
count.images <- function(fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    # get list of folders in given path location
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # check that valid path has been given
    if (length(folders) == 0) {
        cat("Image folder",fpath,"not found in",getwd())
    } else {
        
        d <- length(folders)
        df <- data.frame(black = double(), grey = double(), white = double())
        
        for (i in 1:d) {    
            # folder name should be acquisition date (yymmdd)
            img.date <- folders[i]
            rev.date <- paste(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2), sep = "")
            
            # check folder for each batch in turn, get number of .tif files
            b <- length(list.files(paste(fpath, folders[i],"/black/", sep = ""), pattern = "\\.tif$"))
            g <- length(list.files(paste(fpath, folders[i],"/grey/", sep = ""), pattern = "\\.tif$"))
            w <- length(list.files(paste(fpath, folders[i],"/white/", sep = ""), pattern = "\\.tif$"))
            
            # add .tif count to data frame
            df[i,] <- c(b,g,w)
            row.names(df)[i] <- img.date
        }
        return(df)
    }
}


#' Quickly load files
#'
#' Load all pixelwise summaries and bad pixel map in one step.
#' @param fpath Path to top level of stored images. Default is "./Other-data/"
#' @details Imports pixelwise mean and SD for all files, and bad pixel map
#' @export
#' @examples
#' load.pixel.maps()
#' 
load.pixel.maps <- function(fpath = "./Other-data/") {
    pw.w <<- readRDS(paste0(fpath, "Pixelwise-means-white.rds"))
    pw.g <<- readRDS(paste0(fpath, "Pixelwise-means-grey.rds"))
    pw.b <<- readRDS(paste0(fpath, "Pixelwise-means-black.rds"))
    
    pw.sd.w <<- readRDS(paste0(fpath, "Pixelwise-sds-white.rds"))
    pw.sd.g <<- readRDS(paste0(fpath, "Pixelwise-sds-grey.rds"))
    pw.sd.b <<- readRDS(paste0(fpath, "Pixelwise-sds-black.rds"))
    
    bpm <<- read.csv(paste0(fpath, "BadPixelMap-160314.csv"), as.is = T)
}





#' Import a single day's acquisition into a 2048 x 2048 array
#' 
#' Import all images for one acquisition date into a 2048 x 2048 array, with layers corresponding to different power settings. Cropped pixels are padded with NA.
#' @export
#' @examples
#' im.160430 <- import.acq("./Image-data/160430")
#' 
import.acq <- function(acq.folder, subfolders = c("black", "grey", "white"), panel.dim = c(2048, 2048)) {
    
    acq <- array(dim = c(panel.dim, length(subfolders)),
                 dimnames = list(NULL, NULL, subfolders))
    
    for (sf in subfolders) {
        
        # load xml data into data frame (get image list & offsets)
        im.list <- list.files(paste0(acq.folder, "/", sf), pattern = "\\.xml$", full.names = T)
        df <- rbind.fill(lapply(im.list, 
                                function(im.xml) as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(im.xml)), list(recursive = TRUE))), 
                                                                      FUN = unlist), stringsAsFactors = F)))
        # gives ordered list of images - may be needed if investigating as time series
        df <- data.frame(acq.time = df$Time,
                         filenm = df$.attrs.FileName,
                         offset.x = df$CameraProperties.imageOffsetX,
                         offset.y = df$CameraProperties.imageOffsetY,
                         size.x = df$CameraProperties.imageSizeX,
                         size.y = df$CameraProperties.imageSizeY,
                         kV = df$kV,
                         uA = df$uA,
                         exp.time = df$ExposureTime,
                         stringsAsFactors = F)[order(acq.time = df$Time),]
        
        # quick error check - are all images on same area of detector (should always be T)
        if (!identical(sapply(df[3:6], min), sapply(df[3:6], max))) {
            cat("Error - TIF images do not all have same offset & dimension")
            return(df)
        }
        
        # if all images cover same area of detector, continue
        offset.x <- as.integer(df$offset.x[1])
        offset.y <- as.integer(df$offset.y[1])
        size.x <- as.integer(df$size.x[1])
        size.y <- as.integer(df$size.y[1])
        
        ims <- apply(cbind(acq.folder, "/", sf, "/", df$filenm), 1, paste, collapse = "")
        
        # get pixelwise mean of image
        tmp <- apply(abind(lapply(ims, readTIFF, as.is = T), along = 3), 1:2, mean)
        tmp <- t(tmp[nrow(tmp):1, , drop = FALSE])
        
        # position correctly in array
        acq[(offset.x + 1) : (offset.x+size.x),
            (panel.dim[2]-size.y-offset.y + 1) : (panel.dim[2]-offset.y),
            sf] <- tmp
    }
    return(acq)
}



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






#' Apply shading correction
#' 
#' Calculate a shading-corrected image from a set of blanks
#' @param im Set of pixelwise means for black, grey and white acquisitions
#' 
#' @export
#' 
shading.corrected <- function(im, fix.inf = T) {
    sc <- 60000 * (im[,,"grey"] - im[,,"black"]) / (im[,,"white"] - im[,,"black"])
    
    if (fix.inf) {
        sc[is.infinite(sc)] <- NA
    }
    sc
}










####################################################################################################




#' Create sample of 'healthy' pixels
#'
#' Given a matrix of coordinates of bad pixels, create a same-size random sample of healthy pixels that do not appear on the bad pixel list.
#' @param bad.pixels Data frame or matrix whose first two columns are X and Y coordinates of identified bad pixels.
#' @param n Size of sample to create. Default is to create a sample of same size as bad pixel list.
#' @return Matrix containing X and Y coordinates of sampled pixels.
#' @export
#' @examples
#' bp.160314 <- reset.bp(160314)
#' samp.160314 <- sample.healthy(bp.160314)
#' 
#' 
sample.healthy <- function(bad.pixels, n = nrow(bad.pixels)) {
    
    # convert coords to single integer
    # subtract 1997 to because coords start at 1, not 9
    bp.int <- (bad.pixels[, 1] * 1996) + bad.pixels[, 2] - 1997
    
    # sample as integer (much more efficient than permuting)
    s <- sample(c(1:(1996^2))[-c(bp.int)], n)
    
    # decompose into integer division & modulo components to get x & y
    # add 1 to get coordinates to correct range
    cbind(x = s %/% 1996 + 1, y = s %% 1996 + 1)
}


#' Convert bad pixel map to image
#' 
#' Given a bad pixel map with coordinates and types, create a 2d array showing each type of bad pixel identified
#' @param bpx Bad pixel map: data frame with first two columns containing x and y coordinates of bad pixels, and a column 'type' specifying what kind of bad pixel was identified.
#' @param im.dim Vector of dimensions of image array to create. Default is c(1996, 1996).
#' @return 2d image array of the specified dimensions
#' @export
#' 
bpx2im <- function(bpx, im.dim = c(1996, 1996)) {
    im <- array(0, dim = im.dim)
    im[as.matrix(bpx[,1:2])] <- bpx$type
    return(im)
}





#' Combine two bad pixel maps
#'
#' Merge two data frames of bad pixels into one single list of coordinates, with source map and (where possible) error type marked.
#'
#' NB. Currently only set up to compare 'official' (old) pixel map with new one. Should update this.
#' @param new.map Data frame (as created using \link{\code{reset.bp}}) containing coordinates and types of pixels.
#' @param old.map Data frame (eg. official bad pixel map) containing coordinates and types of pixels.
#' @return Data frame containing row and column coordinates of bad pixels, factor giving type of bad pixel, and factor noting which map the pixel was found in (old, new, or both)
#' @export
#' @examples
#' matched <- combine.maps(bp, bpm)
combine.maps <- function(new.map, old.map) {
    
    match <- merge(data.frame(new.map[!duplicated(new.map[,1:2]), c("row", "col", "type")],
                              map1 = "new"),
                   data.frame(old.map, map2 = "old"),
                   by.x = c(1:2), by.y = c(1:2), all = T)
    
    match$map <- apply(match[,c("map1", "map2")], 1, function(x) paste(na.omit(x), collapse = ""))
    match$map[match$map == "newold"] <- "both" 
    match$map <- factor(match$map, levels = c("new", "both", "old"))
    
    match[order(match$map), c("row", "col", "type", "map")]
}



#' Plot a combined pixel map
#'
#' Compare bad pixels identified by two different pixel maps by plotting their values.
#' @param map Data frame (as created using \link{\code{combine.maps}}) containing coordinates and sources of pixels.
#' @param img.date Integer or string giving acquisition date for which bad pixels are to be plotted.
#' @param colour Which image batch should be plotted? Default is "all", which will plot all three sequentially.
#' @export
#' @examples
#' matched <- combine.maps(bp, bpm)
#' plot.map.comparison(matched)
plot.map.comparison <- function(map, img.date, colour = "all") {
    
    # get random sample of 'healthy' pixels for comparison
    pp <- sample.healthy(map)
    
    img.date <- toString(img.date)
    
    if (colour == "all") {
        
        for (colour in c("black", "grey", "white")) {
            plot(pw.m[,, colour, img.date][pp], pch = 20, col = adjustcolor("gold", alpha = 0.2),
                 ylim = c(0, 65535), ylab = "Mean pixel value over 20 images", xlab = "Pixel index",
                 main = paste0("Bad pixels - ", colour, " images, ", img.date))
            
            points(pw.m[,, colour, img.date][as.matrix(map[,1:2])], pch = 20, 
                   col = adjustcolor(c("chartreuse3", "slateblue1", "red")[map$map], alpha = 0.5))
            
            legend("topleft", legend = c(levels(map$map), "'healthy' pixels"), pch = 20,
                   col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
        }
    } else {
        
        plot(pw.m[,, colour, img.date][pp], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), ylab = "Mean pixel value over 20 images", xlab = "Pixel index",
             main = paste0("Bad pixels - ", colour, " images, ", img.date))
        
        points(pw.m[,, colour, img.date][as.matrix(map[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        
        legend("topleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
    }
}


#' Johnson Q-Q plot
#'
#' Produce a Q-Q plot of data quantiles vs Johnson distribution quantiles, with quantiles of interest marked.
#' @param data Image array of 1996 x 1996 pixels
#' @param quantiles Vector of quantiles to plot; default is c(1:999)/1000
#' @param grid.quantiles Vector of quantiles of particular interest, to mark with gridlines
#' @param title String containing plot title. Default is "Johnson Q-Q plot".
#' @param ... Additional plotting arguments
#' @export
#' @examples
#' Johnson.QQ(pw.sd[,,"black", "150828"])
Johnson.QQ <- function(data, quantiles = c(1:999)/1000, grid.quantiles = c(0.01, 0.99), title = "Johnson Q-Q plot", ...) {
    
    data <- data[!is.na(data)]
    jf <- JohnsonFit(data, moment = "quant")
    
    plot(qJohnson(quantiles, jf), quantile(data, quantiles), pch = 20, asp = T,
         ylab = "Observed quantile", xlab = "Johnson quantile", main = title, ...)
    abline(0,1, col = "red", lty = 2)
    
    abline(h = quantile(data, grid.quantiles), col = "skyblue", lty = 2)
    abline(v = qJohnson(grid.quantiles, jf), col = "skyblue", lty = 2)
    
    legend("topleft", cex = 0.7, title = "Gridlines: ", legend = paste0("Q", grid.quantiles * 100), lty = 2, col = "skyblue")
}







#' Get median differences over image
#' 
#' Run a median filter over an image to obtain a smoothed version, then subtract this from the original image to highlight any points that are unusually high in relation to their neighbours
#' @param im Single-layer image array
#' @return Image array containing pixelwise differences
#' @export
#' 
med.diffs <- function(im) {
    im - r2m(focal(m2r(im), matrix(rep(1, 9), ncol = 3), fun = median))
}


#' Summarise bad pixels for each acquisition
#' 
#' @export
#' 
summarise.bpx <- function(bp, by = "type") {
    
    df <- rbind.fill(lapply(lapply(bp, "[", by), function(m) data.frame(t(as.matrix(table(m))))))
    rownames(df) <- names(bp)
    df
}




#' Per-column loess smoothing
#' 
#' Fit a loess spline to each column of the image using \link{\code{lowess}}. The upper and lower detector panels are treated separately.
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @param span Span to be used for loess smoothing. Default is 1/ sqrt(n), with n taken to be the height of the midline.
#' @param midline Midline of subpanels. Default is 992.5.
#' @return 1996x1996 matrix of smoothed values.
#' @export
#' @examples
#' smoothed <- lowess.per.column(panel.res)
#' smoothed.res <- panel.res - smoothed
#' pixel.image(smoothed.res)
#' s.hist(smoothed.res)
lowess.per.column <- function(image, span, midline = 992.5) {
    
    upper.midline <- ceiling(midline)
    lower.midline <- floor(midline)
    upper.edge <- dim(image)[[2]]
    
    # split image into upper & lower panels: smoothing across midline is nonsensical
    upper <- image[, upper.midline:upper.edge]
    lower <- image[, 1:lower.midline]
    
    # default span: 1/ sqrt(n)
    if(missing(span)) span <- 1/ sqrt(midline)
    
    # apply lowess smoothing per column to each panel in turn
    upper.smoo <- do.call("rbind", lapply(apply(upper, 1, lowess, f = span), "[[", 2))
    lower.smoo <- do.call("rbind", lapply(apply(lower, 1, lowess, f = span), "[[", 2))
    
    # join upper & lower panels again 
    smoo <- array(dim = dim(image))
    smoo[, upper.midline:upper.edge] <- upper.smoo
    smoo[, 1:lower.midline] <- lower.smoo
    
    return(smoo)
}


####################################################################################################


#' Linear regression over distance from centre
#' 
#' Fit a circular spot model with linear gradients to the image, including polynomial terms of z as explanatory variables
#' @param image Single-layer (1996x1996) array image to be fitted.
#' @param order Order of polynomial to include. Default is 2
#' @param centre Vector of coordinates of centre of spot. Default is centre of uncropped detector, c(1023.5, 992.5)
#' @param robust Boolean: use robust model fitting (\link{\code{rlm}})? Default is F (use \link{\code{lm}} for regression).
#' @return Linear model fitted to given image
#' @export
#' @examples
#' circ.lm <- spot.lm(pw.m)
#' circ.res <- matrix(circ.lm$residuals, ncol = 1996)
#' pixel.image(circ.res)
#' 
spot.lm <- function(image, order = 2, centre = c(1023.5, 992.5), robust = F) {
    
    # get distances from spot centre 
    z.dist <- data.frame(x =  c(1:1996^2)%%1996, y = c(1:1996^2)%/%1996 + 1)
    z.dist$x[z.dist$x == 0] <- 1996                                 # correct x-coord
    z.dist$y[z.dist$x == 1996] <- z.dist$y[z.dist$x == 1996] - 1    # correct y-coord
    z.dist$z <- sqrt((z.dist$x - centre[1])^2 + (z.dist$y - centre[2])^2)
    
    # distance from centre vs pixel values
    zz <- cbind(melt(image), z.dist)[,c(4,5,6,3)]
    
    # fit & return linear model
    if (order > 1) {
        if (robust) {
            mod <- rlm(value ~ poly(z, order), zz)
        } else {
            mod <- lm(value ~ poly(z, order), zz)
        }
    } else {
        if (robust) {
            mod <- rlm(value ~ z, zz)
        } else {
            mod <- lm(value ~ z, zz)
        }
    }
    return(mod)
}


#' Fit hyper- or hypo-elliptical spot model
#' 
#' Uses linear regression over x- and y-distance from centre to the power n to fit a hyper- or hypo-elliptical spot
#' @param image Single-layer (1996x1996) array image to be fitted.
#' @param n Power to which x- and y- distances are to be raised to fit model.
#' @param order Order of polynomial to include. Default is 2
#' @param centre Vector of coordinates of centre of spot. Default is centre of uncropped detector, c(1023.5, 992.5)
#' @param robust Boolean: use robust model fitting (\link{\code{rlm}})? Default is F (use \link{\code{lm}} for regression).
#' @return Linear model fitted to given image
#' @details 1 < n < 2: rhombus with convex edges; n = 2: an ellipse; n > 2: rectangle with rounded corners.
#' @export
#' @examples
#' spot <- spot.lm(pw.m[,,"white", "150828"], n = 2, robust = T)
#' spot.res <- matrix(spot$residuals, ncol = 1996)
#' pixel.image(spot.res)
he.spot.lm <- function(image, n = 2, order = 2, centre = c(1023.5, 992.5), robust = F) {
    
    # get x & y coords from spot centre
    z.dist <- data.frame(x =  c(1:1996^2)%%1996, y = c(1:1996^2)%/%1996 + 1)
    z.dist$x[z.dist$x == 0] <- 1996
    z.dist$y[z.dist$x == 1996] <- z.dist$y[z.dist$x == 1996] - 1
    
    zz <- cbind(melt(image), z.dist,
                xx = abs(z.dist$x - centre[1])^n, 
                yy = abs(z.dist$y - centre[2])^n)[, c(4,5,6,7, 3)]
    
    # fit & return linear model
    if (order > 1) {
        if (robust) {
            mod <- rlm(value ~ poly(xx, order) + poly(yy, order), zz)
        } else {
            mod <- lm(value ~ poly(xx, order) + poly(yy, order), zz)
        }
    } else {
        if (robust) {
            mod <- rlm(value ~ xx + yy, zz)
        } else {
            mod <- lm(value ~ xx + yy, zz)
        }
    }
    return(mod)
}


#' Per-panel linear regression
#' 
#' Fit a linear gradient over each of the 32 subpanels within the image
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @param terms String specifying terms to be fitted to each subpanel. Default is "x + y", fitting a linear model to the x and y coordinates of each subpanel without considering interactions.
#' @param robust Boolean: use robust regression (\link{\code{rlm}})? Default is F (use \link{\code{lm}} for regression).
#' @param left.pad Integer value: how many pixels are cropped from left image edge? Default is 2.
#' @param upper.pad Integer value: how many pixels are cropped from upper image edge? Default is 20.
#' @param rotate.lower boolean: fit a model to the lower panels as they appear in the image, or rotate by 180 degrees to align amplifier coordinates with those of top panels? Default is F (fit models to image as is)
#' @return List containing the terms applied, a matrix of fitted values, and a matrix of model coefficients for each panel.
#' @export
#' @examples
#' panel.lm <- panel.lm(circ.res, "poly(x, 2) + y")    # finds coefficients of intercept, x, x^2, y.
#' panel.lm <- panel.lm(circ.res, "y")                 # finds coefficients of intercept and y.
#' panel.res <- circ.res - panel.lm$fitted.values
#' pixel.image(panel.res)
panel.lm <- function (image, terms = "x + y", robust = F, left.pad = 2, upper.pad = 20, rotate.lower = F) {
    
    # convert formula string to lower case to ensure match
    terms <- tolower(terms)
    
    # create empty arrays to hold output
    coeffs <- c()
    smoothed.panels <- array(dim = c(128, 1024, 32))
    
    # iterate over subpanels, fit specified linear model to each
    sp <- subpanels(image)
    
    for (s in 1:32) {
        if (s > 16 & rotate.lower) {
            sp[,,s] <- array(rev(sp[,,s]), dim = dim(sp[,,s]))
        }
        df <- melt(sp[,,s])
        colnames(df) <- c("x", "y", "value")
        
        if (robust) {
            lm <- rlm(as.formula(paste0("value ~ ", terms)), data = df)
        } else {
            lm <- lm(as.formula(paste0("value ~ ", terms)), data = df)
        }
        
        coeffs <- rbind(coeffs, coef(lm))
        smoothed.panels[, , s] <- predict(lm, df)
    }
    list(formula = paste0("value ~ ", terms), fitted.values = join.panels(smoothed.panels, left.pad, upper.pad), models = coeffs)
}


#' Mini-panel linear regression
#' 
#' Fit a linear gradient over each of the 4 minipanels within each of the 32 subpanels within the image
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @param terms String specifying terms to be fitted to each subpanel. Default is "x + y", fitting a linear model to the x and y coordinates of each subpanel without considering interactions.
#' @param robust Boolean: use robust regression (\link{\code{rlm}})? Default is F (use \link{\code{lm}} for regression).
#' @param left.pad Integer value: how many pixels are cropped from left image edge? Default is 2.
#' @param upper.pad Integer value: how many pixels are cropped from upper image edge? Default is 20.
#' @return List containing the terms applied, a matrix of fitted values, and a matrix of model coefficients for each panel.
#' @export
#' @examples
#' panel.lm <- panel.lm(circ.res, "y")                 # finds coefficients of intercept and y.
#' panel.res <- circ.res - panel.lm$fitted.values
#' minipanel <- minipanel.lm(panel.res, terms = "x + y", robust = T)
#' res <- panel.res - minipanel$fitted.values
minipanel.lm <- function (image, terms = "x + y", robust = F, left.pad = 2, upper.pad = 20) {
    terms <- tolower(terms)
    coeffs <- c()
    smoothed.panels <- array(dim = c(128, 1024, 32))
    sp <- subpanels(image)
    for (s in 1:32) {
        for(i in 1:4) {
            df <- melt(sp[, (1024 - (i * 256) + 1): (1024 - ((i-1) * 256)), s])
            colnames(df) <- c("x", "y", "value")
            
            if (robust) {
                lm <- rlm(as.formula(paste0("value ~ ", terms)), data = df)
            } else {
                lm <- lm(as.formula(paste0("value ~ ", terms)), data = df)
            }
            coeffs <- rbind(coeffs, coef(lm))
            smoothed.panels[, (1024 - (i * 256) + 1): (1024 - ((i-1) * 256)), s] <- predict(lm, df)
        }
    }
    list(formula = paste0("value ~ ", terms),
         fitted.values = join.panels(smoothed.panels, left.pad, upper.pad), 
         models = coeffs)
}






#' Return fitted values from panelwise coefficients
#' 
#' Given a matrix of coefficients such as from \link{\code{fit.panel.lm}}, return a list containing values fitted to the panel coordinates.
#' @param coeffs 32 x 3 matrix of coefficients (offset, x, y) for each panel
#' @return List containing 1996x1996 matrix of per-panel offset values and 1996x1996 matrix of per-panel gradient values.
#' @export
#' @examples
#' panel.lm <- fit.panel.lm(circ.res)
#' zz <- fit.panels(panel.lm$models)
#' pixel.image(zz$grad)
fit.panels <- function(coeffs) {
    
    offset.p <- array(dim = c(128, 1024, 32))
    grad.p <- array(dim = c(128, 1024, 32))
    
    for (i in 1:32) {
        offset.p[,,i] <- coeffs[i,"offset"]
        tmp <- melt(grad.p[,,i])
        grad.p[,,i] <- (coeffs[i,"x"] * tmp$X1) + (coeffs[i,"y"] * tmp$X2)
    }
    
    list(offset = join.panels(offset.p),
         grad = join.panels(grad.p))
}



#' Split image array into subpanels
#'
#' Take a single image layer and return 32-layer labelled subpanel array to allow easier computation over panels.
#' @param data Single-layer 1996x1996 image array
#' @param left.pad Integer numer of pixels cropped from left-hand edge of image. Default is 2
#' @param upper.pad Integer numer of pixels cropped from upper edge of image. Default is 20
#' @export
#' @examples
#' load.images(150828, "black")
#' pw.m <- pixelwise.mean(b.150828)
#' 
#' zz <- subpanels(pw.m)
#' pixel.image(zz[,,"U8"], break.levels = sd.levels(b.150828))
#' 
subpanels <- function(data, left.pad = 2, upper.pad = 20) {
    
    # only works on single layer - unwieldy when larger
    
    panel.names <- apply(cbind(c(rep("U", 16), rep("L", 16)),
                               rep(c(1:16), 2)), 
                         1, paste, collapse = "")
    
    # pad original image with NA to bring to 2048 x 2048 square
    n <- array(NA, dim = c(2048, 2048))
    n[left.pad + c(1:1996), 2048-1996-upper.pad + c(1:1996)] <- data
    
    # create array to hold split data
    m <- array(dim = c(128, 1024, 32), 
               dimnames = list(dimnames(data)[[1]], dimnames(data)[[2]], panel.names))
    
    # extract panels
    for (i in 1:16) {
        m[ , , i] <- n[(128 * (i-1)) + c(1: 128), 1025:2048]        # upper row
        m[ , , i + 16] <- n[(128 * (i-1)) + c(1: 128), 1:1024]      # lower row
    }
    return(m)
}


#' Merge subpanels into image array
#' 
#' Take a subpanel image array and recombine into single panel image
#' @param data 128x1024x32 array of subpanel images
#' @param left.pad Integer numer of pixels cropped from left-hand edge of image. Default is 2
#' @param upper.pad Integer numer of pixels cropped from upper edge of image. Default is 20
#' @export
#' @examples
#' 
#' zz <- split.panels(pw.m)
#' qq <- join.panels(zz)
#' all(pw.m == qq)               # TRUE
join.panels <- function(data, left.pad = 2, upper.pad = 20) {
    m <- array(dim = c(2048, 2048), dimnames = list(dimnames(data)[[1]], dimnames(data)[[2]]))
    
    # extract panels
    for (i in 1:16) {
        m[(128 * (i-1)) + c(1: 128), 1025:2048] <- data[ , , i]         # upper row
        m[(128 * (i-1)) + c(1: 128), 1:1024] <- data [ , , i + 16]      # lower row
    }
    
    # return valid 1996x1996 region
    # (can't simply remove NA - after smoothing, NA may exist within panel)
    return(m[left.pad + c(1: 1996), 2048-1996-upper.pad + c(1:1996)])
}



#' Calculate bivariate normal density of image surface
#' 
#' Calculate the bivariate Gaussian density of an image surface with given parameters (rho fixed at 0)
#' @param par Vector of named parameters A (constant offset), x0 (x-origin), y0 (y-origin), sig.x (SD over x), sig.y (SD over y)
#' @param obs Data frame with columns x and y containing coordinates at which to evaluate the function. 
#' @return Vector containing calculated density
#' @export
#' 
#' @examples
#' gv <- setNames(melt(pw.m[,,"grey", "160430"]), nm = c("x", "y", "z"))
#' tmp <- bvn(c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500), gv)
#' pixel.image(array(tmp, dim = c(2048, 2048)))
#' 
bvn <- function(par, obs) {
    A <- par["A"]; x0 <- par["x0"]; y0 <- par["y0"]
    sig.x <- par["sig.x"]; sig.y <- par["sig.y"]
    
    A * exp(- 0.5 * ((((obs$x - x0) / sig.x)^2) + ((obs$y - y0) / sig.y)^2))
}



#' Support function: get sum of squared errors of bivariate normal density
#' 
#' Function to calculate bivariate Gaussian density (with rho = 0) with given parameters for image surface, and return sum of squared residuals. Primarily a support function for \link{\code{optim}}.
#' @param par Vector of named parameters, to be passed to function to be evaluated
#' @param obs Data frame containing columns named x, y and z, containing coordinates at which function should be evaluated, and observed values at those points.
#' @param fn Function to be evaluated
#' @export
#' 
bvn.ss <- function(par, obs) {
    
    est <- bvn(par, obs)
    sum((est - obs$z)^2, na.rm = T)
}



#' Least Squares fit of bivariate Gaussian surface to image
#' 
#' Use \link{\code{optim}} to fit an elliptical bivariate Gaussian surface to a two-dimensional image by minimizing sum of squared errors.
#' @param im 2d image array to which surface is to be fitted.
#' @param par Starting values for the parameters to be fitted.
#' @param x0.target Range of values within which the x-origin will be constrained to lie. Default is image area.
#' @param y0.target Range of values within which the y-origin will be constrained to lie. Default is image area.
#' @return List of fitted parameters and convergence statistics produced by \link{\code{optim}}.
#' @export
#' 
#' @examples
#' spot <- gaussian.spot.ls(pw.m[,,"grey", "160430"], c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
#'
gaussian.spot.ls <- function(im, par, x0.target = c(0, ncol(im)), y0.target = c(0, nrow(im))) {
    
    # convert image to data.frame of coordinates & observations
    gv <- setNames(melt(im), nm = c("x", "y", "z"))
    
    # optimise by minimising sum of squared errors of bivariate normal distribution
    optim(par = par, bvn.ss, obs = gv, method = "L-BFGS-B",
          lower = c(-Inf, x0.target[1], y0.target[1], 0, 0), 
          upper = c(Inf, x0.target[2], y0.target[2], Inf, Inf))
}



#' Produce bivariate density surface
#' 
#' Produce bivariate density surface from given parameters and return as array
#' @param par Vector of named parameters
#' @param arr.dim Dimensions of array to be produced. Default is c(2048, 2048).
#' @return Array of given dimensions containing evaluated surface at each point.
#' @export
#' 
#' @examples
#' fitted.spot <- gaussian.spot.ls(pw.m[,,"grey", "160430"], c(A = 15000, x0 = 1024.5, y0 = 1024.5, sig.x = 500, sig.y = 500))
#' pixel.image(gaussian.spot.mat(fitted.spot$par, arr.dim = c(2048, 2048)))
#' 
gaussian.spot.mat <- function(par, arr.dim = c(2048, 2048)) {
    
    obs <- setNames(melt(array(dim = arr.dim)), nm = c("x", "y", "z"))
    est <- bvn(par, obs)
    
    array(est, dim = arr.dim)
}



#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at mean and +- 0.5, 1, 2, 3 sd.
#' @param data Matrix of data to plot
#' @param midpoint String: use mean or median as midpoint? Default is "mean".
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'                    data,
#'                    levels = sd.levels(data))
#' 
#' 
sd.levels <- function(data, midpoint = "mode") {
    
    data <- data[!is.na(data)]
    data <- data[!is.infinite(data)]
    
    mfn <- switch(midpoint,
                  mean = "mean",
                  median = "median",
                  mode = "modal.density")
    
    m <- eval(parse(text = paste0(mfn, "(data)")))
    
    sort(c(min(data), m + (c(-6, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 6) * sd(data)), max(data)))
}



#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at mean and +- 0.5, 1, 2, 3 sd.
#' @param data Matrix of data to plot
#' @param midpoint String: use mean or median as midpoint? Default is "mean".
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'                    data,
#'                    levels = sd.levels(data))
#' 
#' 
th.levels <- function(data, fix.centre = NA) {
    
    data <- data[!is.na(data)]
    data <- data[!is.infinite(data)]
    
    if(is.na(fix.centre)) {
        m <- modal.density(data)
    } else {
        m <- fix.centre
    }
    
    sort(c(sapply(c(1, 2, 3, 4.5, 6), function(nn) asymmetric.mad(data, nn, fix.centre = fix.centre)),
           m,
           range(data)))
}



#' Colour scheme for contour plot
#'
#' Support function: returns a vector of colours to be used in contour plotting.
#' @return Vector of colours. When used in conjunction with  \link{\code{sd.levels()}}, 
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'     data,
#'     levels = sd.levels(mean(data), sd(data)),
#'     col = sd.colours())
#' 
#' 
sd.colours <- function() {
    c("black",             # min to -6sd        # min to -6mad
      "darkblue",     # -6 to -3sd              # -6 to -5mad
      "blue",             # -3 to -2sd          # -5 to -3mad
      
      "cyan3",      # -2 to -1sd                # -3 to -2mad
      
      "green1",       # -1 to -0.5sd            # -2 to -1mad
      "greenyellow",    # 0.5 to 0sd            # -1 to 0mad
      "yellow",             # 0 to 0.5sd        # 0 to 1mad
      "gold",           # 0.5 to 1sd            # 1 to 2mad
      
      "orange",              # 1 to 2sd         # 2 to 3mad
      
      "red3",        # 2 to 3sd                 # 3 to 5mad
      "purple",           # 3 to 6sd            # 5 to 6mad
      "deeppink")         # 6sd to max          # 6mad to max  
}



#' Display pixel values as image
#'
#' Create an image with pixels shaded according to their distance from the mean value. Displays at lower resolution than contour plot (plots are 20% of the size), but much quicker (10% of time taken).
#' @details Pixels more than 2sd above the mean are red, pink or purple. Pixels within 2sd of the mean are orange, yellow or green. Pixels more than 2sd below the mean are blue.
#' @param data 2-dimensional matrix containing values to be plotted
#' @param x.range Vector range showing x-range of cells to be included in plot
#' @param y.range Vector range showing y-range of cells to be included in plot
#' @param title String containing title to be printed with plot
#' @param midpoint Either 'mean' or 'median'. Default is 'mean'
#' @param break.levels Vector of values used as breakpoints when binning values
#' @export
#' @examples
#' pixel.image(pw.mean)
#' 
#' 
pixel.image <- function(data, title = "", x.range = c(1:nrow(data)), y.range = c(1:ncol(data)), 
                        break.levels = th.levels(data), 
                        panels = F, x.lab = "", y.lab = "", ...) {
    
    image(x.range, y.range, 
          data[x.range, y.range], 
          col = sd.colours(),
          breaks = break.levels,
          main = title, xlab = x.lab, ylab = y.lab,
          asp = T, 
          ...)
    
    if (panels) draw.panels()
}


#' Draw outlines of convex hulls
#' 
#' Given a set of coordinates, identify each cluster of adjacent pixels and draw the convex hull.
#' @param px Matrix or data frame of x and y coordinates to be clumped and outlined
#' 
#' @export
#' 
draw.outlines <- function(px, im.dim = c(2048, 2048), ...) {
    
    # clump adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = im.dim)), dir = 4)
    
    # coordinates of each clump, with clump id
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    invisible(lapply(unique(xy$id),
                     function(i) {
                         ch <- chull(xy[xy$id == i,1:2])
                         lines(xy[xy$id == i,1:2][c(ch, ch[1]),], ...)
                     }))
}



#' Plot defective pixels
#'
#' Plot bad pixel map on 2048 x 2048 grid
#' @param px Coordinates of pixels to be plotted
#' @export
#' @examples
#' pixel.plot(which(pw.mean > 25000, arr.ind = T))
#' 
#' 
pixel.plot <- function(px, xlim = c(0,2048), ylim = c(0,2048), pch = 15, panels = F, cex = 0.4,
                       main = "", xlab = "", ylab = "", ...) {
    
    plot(px[,1:2], asp = T, 
         xlim = xlim, ylim = ylim, pch = pch, main = main, xlab = xlab, ylab = ylab, cex = cex, ...)
    
    if (panels) draw.panels()
}


####################################################################################################


#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at outlier (lq - 1.5*IQR), lq, 3 even bins between lq and median, median, 3 even bins between median and uq, uq, and outlier (1.5*iqr + uq).
#' @param data Matrix of data to plot
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'                    data,
#'                    levels = iqr.levels(data))
#' 
#' 
iqr.levels <- function(data) {
    lq <- quantile(data, 0.25)
    med <- median(data)
    uq <- quantile(data, 0.75)
    iqr <- IQR(data)
    
    c(0,
      max(lq - 1.5 * iqr,         # low outliers
          lq/2),
      lq,
      lq + (med - lq) / 3,
      med - (med - lq) / 3,
      med,
      med + (uq - med) / 3,
      uq - (uq - med) / 3,
      uq,
      min(uq + 1.5 * iqr,    # high outliers
          uq + ((65535-uq) / 2)),
      65535)
}





#' Colour scheme for bad pixel plot
#'
#' Support function: returns a vector of colours to be used in bad pixel plotting.
#' @param block Vector of categories to omit from the plot, which will be set to colour NA. Default is c("edge", "s.bright").
#' @return Vector of colours.
#' @export
#' 
bp.colours <- function(block = c("edge", "s.bright")) {
    
    cats <- c("no response", "dead", "hot", "v.bright","bright", "s.bright", "screen spot", "edge", "v.dim", "dim", "s.dim")
    bp.cols <- c("purple", "black", "red", "orange", "gold", "yellow", "grey", "darkgrey", "green3", "green", "lightskyblue")
    bp.cols[which(cats %in% block)] <- NA
    return(bp.cols)
}


#' Plot bad pixels
#' 
#' Plot the coordinates and types of a bad pixel map
#' @param bpm Data frame containing xy coordinates of pixels to be plotted, and a field "type" specifying the category of each bad pixel identified.
#' @param block Vector of categories to omit from the plot, which will be set to colour NA. Default is c("edge", "s.bright").
#' @param ... Additional optional graphical parameters to pass to plotting function
#' @export
#' 
plot.bad.px <- function(bpm, pch = 15, xlab = "", ylab = "", block = c("edge", "s.bright"), ...) {
    plot(bpm[,1:2], pch = pch, col = bp.colours(block = block)[bpm$type], asp = T, xlab = xlab, ylab = ylab, ...)
}

#' Contour plot of pixel values
#'
#' Create a filled contour plot with pixels shaded according to their distance from the mean value. Takes approx. 10x as long as \link{pixel.image}.
#' @details Pixels more than 2sd above the mean are red, pink or purple. Pixels within 2sd of the mean are orange, yellow or green. Pixels more than 2sd below the mean are blue.
#' @param data 2-dimensional matrix containing values to be plotted
#' @param title String containing title to be printed with plot
#' @param midpoint Either 'mean' or 'median'. Default is 'mean'
#' @export
#' @examples
#' pixel.contour(pw.mean)
#' 
#' 
pixel.contour <- function(data, title = "", midpoint = "mean") {
    
    # if data has more than 2 dimensions, truncate to take only first 'layer'
    d <- length(dim(data))
    if (d > 2) {
        data <- data[,, rep(1, d-2)]
    }
    
    filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
                   data,
                   levels = sd.levels(data, midpoint),
                   col = sd.colours(),
                   asp = T,
                   main = title
    )
}








#' Get coordinates of panel edges
#'
#' Identify coordinates of individual detector panels, given their dimensions and those of the whole sensor.
#' @details Each detector panel outputs 128 x 1024 pixels, with some cropping at each edge of the combined image.
#' @param left.crop Number of pixels cropped from left-hand edge of image. Default is 2.
#' @param upper.crop Number of pixels cropped from lower edge of image. Default is 20.
#' @param width Width of each panel, in pixels. Default is 128.
#' @param height Height of each panel, in pixels. Default is 1024.
#' @param x.dim Width of image returned, in pixels. Default is 1996.
#' @param y.dim Height of image returned, in pixels. Default is 1996.
#' @export
#' @examples
#' panels <- panel.edges()
#' 
#' 
panel.edges <- function(left.crop = 0, upper.crop = 0, width = 128, height = 1024, x.dim = 2048, y.dim = 2048) {
    
    #    seq(from = 0, to = y.dim, by = height) + y.dim - height + upper.crop + 1
    
    list(y = c(1, y.dim - height + upper.crop + 1, y.dim + 1),
         x = c(1, c(1:15)* width - left.crop + 1, x.dim + 1))
}


#' Add panel edges to plot
#'
#' Mark edges of individual panels in array
#' @details Each detector panel outputs 128 x 1024 pixels, with some cropping at each edge of the combined image.
#' @param p List containing y, a named vector of y-coordinates of panel corners, and x, a named vector of x-coordinates of panel corners.
#' @export
#' @examples
#' pixel.image(b.150828)
#' draw.panels()            # marks specified panel boundaries
#' 
#' 
draw.panels <- function(p = panel.edges(), ...) {
    
    # horizontal
    for (yl in p$y[2:(length(p$y)-1)]-0.5) {
        lines(range(p$x)-0.5, c(yl, yl), ...)
    }
    
    # vertical
    for (xl in p$x[2:(length(p$x)-1)]-0.5) {
        lines(c(xl, xl), range(p$y)-0.5, ...)
    }
}


#' Vector overlay plot
#'
#' Base plotting function, specifying type "o" (points plotted over lines), pch = 20, cex = 0.7. Other plotting parameters can be passed as arguments as usual.
#' @param data Vector of values to plot (eg. column of pixel values)
#' @param add Boolean: add values to current plot, or start a new plot? Default is F.
#' @export
#' @examples
#' o.plot(b.150828[4,,1])
#' o.plot(b.150828[5,,1], col = adjustcolor("red", alpha = 0.2), add = T)
#' 
#' 
o.plot <- function(data, add = F, ...) {
    if (add) {
        points(data, type = "o", pch = 20, cex = 0.7, ...)
    } else {
        plot(data, type = "o", pch = 20, cex = 0.7, ...)
    }
}


#' Get coordinates of points to plot zoom
#' 
#' Given a set of points and the surrounding area to plot, obtain coordinates to create a plot of a small subset of an image
#' @param points Matrix or data frame, the first two columns of which are the x and y coordinates of the points to be plotted.
#' @param surround Number of additional pixels to plot in every direction around the selected points
#' @return matrix of coordinates to plot
#' @export
#'  
get.focus <- function(points, surround = 3) {
    n <- 2 * surround + 1
    focus <- matrix(c(round(mean(points$x), 0) + rep(c(-surround: surround), n),
                      round(mean(points$y), 0) + sort(rep(c(-surround: surround), n))), ncol = 2)
    focus[focus <= 0] <- 0
    focus[focus >= 2048] <- 2048
    return(focus)
}


#' Histogram plot with SD levels shown
#' 
#' Produce a histogram with coloured bars at the bottom showing scale used. Defaults can be used to show legend for plots produced using \link{\code{pixel.image}}.
#' @param data Array of numbers to be plotting in histogram
#' @param scale Vector of scale cutpoints. Default is \link{\code{sd.levels}}.
#' @param scale.colours Vector of scale colours. Default is \link{\code{sd.colours}}
#' @param xlim Vector of lower and upper x-limits. If not provided, will use central 0.95 of normal distribution with mean & sd of observed data.
#' @export
#' @examples
#' pixel.image(b.150828)
#' scale.hist(b.150828)
s.hist <- function(data, scale = sd.levels(data), scale.colours = sd.colours(), xlim, scale.height = 600, ...) {
    
    data <- data[!is.na(data)]
    
    if (missing(xlim)) {
        xlim <- c(floor(qnorm(0.025, mean(data), sd(data))/10)*10,
                  ceiling(qnorm(0.975, mean(data), sd(data))/10)*10)
    }
    
    hist(data, breaks = "fd", xlim = xlim, ...)
    
    # add colours to indicate scale of pixel map
    cl <- cut(xlim[1]:xlim[2], scale)
    points(xlim[1]:xlim[2], rep(-scale.height, length(xlim[1]:xlim[2])), pch = 15, col = scale.colours[cl])
}


#' Plot focal area
#' 
#' Create pixel image of area surrounding a point of interest. Pixel values are labelled and bad pixels identified highlighted
#' @param im 2d image array to plot
#' @param centre Vector of x and y coordinates of pixel of particular interest
#' @param surround Integer: set size of area surrounding pixel of interest to be displayed. Default is 5, which displays an 11x11 square centred on the pixel of interest.
#' @param dp Integer: display pixel values to how many decimal places? Default is 1.
#' @param scale.by Integer: divide pixel values by this number for easier display and comparison. Default is 1000, so a value of 65535 will be displayed as 65.5
#' @param lbl.cex Magnification factor to apply to labels in each cell. Default is 0.7
#' @param bad.px Optional bad pixel map to be used to highlight any pixels in the display area that have been identified as defective.
#' @param bpx.cex Magnification factor to apply to symbol used to highlight bad pixels. Default is 2.5
#' @param ... Additional optional graphical arguments
#' @export
#' 
focal.plot <- function(im, centre, surround = 5, dp = 1, scale.by = 1000, lbl.cex = 0.7, bad.px, bpx.cex = 2.5, ...) {
    
    ff <- get.focus(data.frame(x = centre[1], y = centre[2]), surround = surround)
    
    pixel.image(im, xlim = range(ff[,1]), ylim = range(ff[,2]), ...)
    text(ff, labels = round(im[ff]/scale.by, dp), cex = lbl.cex)
    
    if (!missing(bad.px)) {
        points(bad.px[,1:2], pch = 0, cex = bpx.cex)
    }
}


#' Draw subpanels on 2048 x 2048 array
#' 
#' Marks locations of 32 subpanels of 128 x 1024 pixels on plot
#' @param ... Optional graphical arguments
#' @export
#' 
draw.panels.2048 <- function(...) {
    
    # horizontal midline
    lines (c(0, 2048), c(1024.5, 1024.5), ...)
    
    for (i in 1:15) {
        lines(c(128 * i + 0.5, 128 * i + 0.5), c(0, 2048), ...)
    }
}


#' Colour ramp for pixel images
#' 
#' @export
sd.ramp <- function() {
    colorRampPalette(c("darkblue", "cyan3", "green2", "lemonchiffon", "gold", "red3", "purple"), 
                     space = "Lab")
}







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
screen.spots.old <- function(im, min.diam = 5, smooth.span = 1/15, midline = 1024.5, edge.crop = 10, coords = T) {
    
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


#' convert coefficients to matrix
#'
#' Create a matrix for use in contour plotting, based on a PPM (point process model)
#' @export
#' 
ppmfit.matrix <- function(ppm) {
    
    zz <- setNames(melt(array(dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))), nm = c("x", "y", "z"))
    cc <- coef(ppm)
    
    zz$z <- cc[1] + (zz$x * cc["x"]) + (zz$y * cc["y"]) + 
        (zz$x^2 * cc["I(x^2)"]) + (zz$y^2 * cc["I(y^2)"]) + (zz$x * zz$y * cc["I(x * y)"]) 
    
    array(zz$z, dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))
}



#' Get transition matrices for bad pixel maps
#' 
#' Given a list of bad pixel maps in coordinate form, return a list of transition matrices showing movements between states observed between one acquisition and the next.
#' @param bpx List of dataframes containing bad pixel coordinates and types
#' @return List of tables showing movement between states from one acquisition to the next
#' @export
#' 
get.transitions <- function(bp) {
    tr <- list()
    
    lvls <- c("normal", levels(bp[[1]]$type))
    
    for (i in 1:(length(bp) - 1)) {
        tr[[i]] <- table("From" = ordered(lvls[c(bpx2im(bp[[i]])) + 1], levels = lvls),
                         "To" = ordered(lvls[c(bpx2im(bp[[i+1]])) + 1], levels = lvls))
        names(tr)[[i]] <- paste(names(bp)[c(i,(i+1))], collapse = "-")
    }
    tr
}



#' Flat-field/offset correction of image
#' 
#' Use full set of black, white and grey images to obtain offset-corrected image. Requires array \code{pw.m}, as generated by \link{\code{load.pixel.means()}}
#' @param dt Integer/string: date to create flat-field image
#' @return 1996x1996 array containing corrected image
#' @export
#' @examples
#' corr <- flat.field.corrected(160314)
flat.field.corrected <- function(dt) {
    dt <- toString(dt)
    
    corr <- (pw.m[,,"grey", dt] - pw.m[,,"black", dt]) /
        (pw.m[,,"white", dt] - pw.m[,,"black", dt]) * 60000
    
    # replace NA (where white == black)
    corr[is.na(corr)] <- pw.m[,,"black", dt][is.na(corr)]
    corr
}



#' Get pixelwise mean of daily acquisitions
#'
#' Return pixelwise mean value across an array
#' @param data Array with 3 or more dimensions
#' @return Array containing pixelwise mean of all acquired values
#' @export
#' @examples
#' pw.mean <- pixelwise.mean(daily)
#' 
#' 
pixelwise.mean <- function(data) {
    # alter this to allow pixelwise mean over 4d array?
    apply(data, c(1,2), mean)
}


#' Get pixelwise standard deviation of daily acquisitions
#'
#' Return pixelwise SD across an array
#' @param data Array with 3 or more dimensions
#' @return Array containing pixelwise standard deviation of all acquired values
#' @export
#' @examples
#' pw.sd <- pixelwise.sd(daily)
#' 
#' 
pixelwise.sd <- function(data) {
    # alter this to allow pixelwise sd over 4d array?
    apply(data, c(1,2), sd)
}



#' Get summary statistics for an image batch
#'
#' Return a vector of summary statistics for a batch of acquisitions: mean, standard deviation, median, upper and lower quartiles, interquartile range, and max and min pixelwise mean value.
#' @param data Array with 3 or more dimensions
#' @return Vector containing named summary statistics
#' @details Pixelwise max and min value are used to check if image is scaled from absolute 0:1 or from an arbitrary minimum point
#' @export
#' @examples
#' b.150828.summ <- batch.summary(b.150828)
#' 
#' 
batch.summary <- function(data) {
    
    pw.sd <- pixelwise.sd(data)
    pw.rng <- apply(data, c(1, 2), max) - apply(data, c(1, 2), min)
    
    # output vector of summary statistics
    list(c(mean = mean(data),
           sd = sd(data),
           pw.sd = mean(pw.sd),
           pw.rng = mean(pw.rng)),
         c(min = min(data),
           lq = quantile(data,0.25),
           median = median(data),
           uq = quantile(data, 0.75),
           max = max(data),
           iqr = IQR(data)))
}


#' Store summary statistics for an image batch in a data frame
#'
#' Populate a data frame with the summary statistics provided by \link{batch.summary}
#' @param data Array with 3 or more dimensions
#' @return Updates global data frame containing summary statistics for multiple images.
#' @export
#' @examples
#' save.summary(b.150828)
#' 
#' 
save.summary <- function(data) {
    
    obj.nm <- toString(as.list(sys.call())[2])
    img.date <- substring(obj.nm,3,8)
    
    
    if (!exists("img.summs")) {
        img.summs <<- data.frame(obj.nm = character(),
                                 img.date = character(),
                                 batch = character(),
                                 g.mean = numeric(),
                                 g.sd = numeric(),
                                 g.lq = numeric(),
                                 g.median = numeric(),
                                 g.uq = numeric(),
                                 g.iqr = numeric(),
                                 pw.min = numeric(),
                                 pw.max = numeric(),
                                 stringsAsFactors = F)
        r <- 1
    } else {
        if (obj.nm %in% img.summs$obj.nm) {
            cat("Row already exists. Delete it manually if you want to overwrite.")
            return()
        } else {
            r <- nrow(img.summs) + 1
        }
    }
    
    img.date <- substring(obj.nm,3,8)  
    
    batches <- matrix(c("b", "g", "w", "black", "grey", "white"), ncol = 2, nrow = 3)
    batch <- batches[which(batches[,1] == substring(obj.nm,1,1)),2]
    
    s <- batch.summary(data)
    
    img.summs[r, c(1:3)] <<- c(obj.nm, img.date, batch)
    img.summs[r,c(4:11)] <<- s
}




#' Create map of 'absolute' bad pixels
#'
#' Create a data frame containing coordinates of 'hot' and 'dead' pixels in each image batch for a given date.
#' 
#' 'Hot' pixels are always on, with a pixelwise mean value of 65535 across all 20 images. 
#' 'Dead' pixels are always off, with a pixelwise mean value of 0 across all 20 images. 
#' 
#' @details Requires objects \code{pw.b}, \code{pw.g} and \code{pw.g}, as created by \link{\code{load.pixel.maps}}.
#' @param img.date Integer or string giving date to be 
#' @return Data frame containing coordinates and classifications of all absolutely defective pixels.
#' @export
#' @examples
#' bp.160314 <- reset.bp(160314)
#' bp.160314 <- reset.bp("160314")          # gives same result
#' 
reset.bp <- function(img.dt) {
    img.dt <- toString(img.dt)
    
    if (exists("pw.m")) {
        bp <- rbind(data.frame(which(pw.m[ , , "black", img.dt] == 0, arr.ind = T), src = "black", type = "dead"),
                    data.frame(which(pw.m[ , , "black", img.dt] == 65535, arr.ind = T), src = "black", type = "hot"),
                    data.frame(which(pw.m[ , , "grey", img.dt] == 0, arr.ind = T), src = "grey", type = "dead"),
                    data.frame(which(pw.m[ , , "grey", img.dt] == 65535, arr.ind = T), src = "grey", type = "hot"),
                    data.frame(which(pw.m[ , , "white", img.dt] == 0, arr.ind = T), src = "white", type = "dead"),
                    data.frame(which(pw.m[ , , "white", img.dt] == 65535, arr.ind = T), src = "white", type = "hot"))
    } else {
        bp <- rbind(data.frame(which(pw.b[ , , img.dt] == 0, arr.ind = T), src = "black", type = "dead"),
                    data.frame(which(pw.b[ , , img.dt] == 65535, arr.ind = T), src = "black", type = "hot"),
                    data.frame(which(pw.g[ , , img.dt] == 0, arr.ind = T), src = "grey", type = "dead"),
                    data.frame(which(pw.g[ , , img.dt] == 65535, arr.ind = T), src = "grey", type = "hot"),
                    data.frame(which(pw.w[ , , img.dt] == 0, arr.ind = T), src = "white", type = "dead"),
                    data.frame(which(pw.w[ , , img.dt] == 65535, arr.ind = T), src = "white", type = "hot"))
    }
    bp
}

#' Identify 'bad' pixels
#'
#' Use thresholds set out in technical manual to identify 'bad' pixels
#' @param data Array containing a sequence of images from a single batch
#' @param too.bright Numeric value to use as threshold for underperforming bright pixels: the multiple of median value above which a pixel is deemed to be unacceptably bright. Default is 1.5 x median.
#' @param too.dim Numeric value to use as threshold for underperforming dark pixels: the multiple of median value below which a pixel is deemed to be unacceptably bright. Default is 0.45 x median.
#' @param noisy Numeric value to use as threshold for noisy pixels: the number of median pixel SDs above which a pixel is deemed to be noisy. Default is 6.
#' @return Data frame containing coordinates and classifications of all underperforming pixels in the given array.
#' @export
#' @examples
#' bp.b.150828 <- threshold.pixels(b.150828)
#' 
#' 
threshold.pixels <- function(data, too.bright = 1.5, too.dim = 0.45, noisy = 6) {
    
    # create temporary list to hold output
    bp <- list()
    
    pw.mean <- pixelwise.mean(data)
    pw.sd <- pixelwise.sd(data)
    median.value <- median(data)
    median.sd <- median(pw.sd)
    
    
    # always use tryNULL when creating data frames, in case no rows are returned
    
    # order in which each type is added to df is order of priority when classifying
    # (eg. duplicates in later categories will be removed)
    
    # dead pixels: daily mean value is exactly 0
    bp$dead <- tryNULL(cbind(data.frame(which(pw.mean == 0.0, arr.ind = T)), 
                             type = "Dead"))
    
    # hot pixels: daily mean value is exactly 65535
    bp$hot <- tryNULL(cbind(data.frame(which(pw.mean == 65535.0, arr.ind = T)),
                            type = "Hot"))
    
    # underperforming bright pixels: value > 1.5x median
    bp$bright <- tryNULL(cbind(data.frame(unique(which(data > (too.bright * median.value), arr.ind = T)[,c(1:2)])),
                               type = "Too bright"))
    if (!is.null(bp$bright)) {colnames(bp$bright) <- c("row", "col", "type")}
    
    # underperforming dark pixels: value < 0.45x median
    bp$dim <- tryNULL(cbind(data.frame(unique(which(data < (too.dim * median.value), arr.ind = T)[,c(1:2)])),
                            type = "Too dark"))
    if (!is.null(bp$dim)) {colnames(bp$dim) <- c("row", "col", "type")}
    
    # noise: pixel sigma > 6x median sigma
    bp$noisy <- tryNULL(cbind(data.frame(which(pw.sd > (noisy * median.sd), arr.ind = T)),
                              type = "Noisy"))
    
    # merge separate data frames into one
    df <- data.frame(do.call("rbind", bp))
    
    # remove duplicates and return df
    return(df[!duplicated(df[,c(1:2)]),])
}


#' Identify extreme values along transect
#' 
#' Find extreme pixel values in single vector
#' @details Uses Loess smoothing to find 'spine' of vector and identifies any points that fall more than a certain constant multiple of the SD from the central spine.
#' @param transect Vector of pixel values to be assessed for 'bad' pixels
#' @param loess.span Filter span to be used in smoothing the data. Default is 1/15.
#' @param keep Central proportion of points to keep. Default is 0.95, which results in a cutoff for 'bad' pixels of 1.96sd.
#' @param display Boolean: plot the resulting transect? Default is T.
#' @param legend Boolean: include legend in plot? Default is F.
#' @param legend.pos String: where to display legend? default is "topright".
#' @param ... Additional parameters to be passed to plotting function
#' @return A list containing named vectors of Loess-smoothed values, residuals, and indices of pixels flagged.
#' @export
#' @examples 
#' zz <- bp.by.res.sd(pw.m[992, 1:992], display = T)
bp.by.res.sd <- function(transect, loess.span = 1/15, keep = 0.95, display = T, legend = F, legend.pos = "topright",  ...) {
    
    k <- qnorm(1-(1-keep)/2, 0, 1)
    
    smoo <- lowess(transect, f = loess.span)$y
    res <- transect - smoo
    
    abs.mean <- mean(abs(res))
    abs.sd <- sd(abs(res))
    
    bp <- which(abs(res) > (abs.mean + k * abs.sd))
    
    if (display) {
        plot(transect, type = "l", col = adjustcolor("grey", alpha = 0.5), ...)
        points(transect, pch = 20, cex = 0.7)
        
        points(smoo, type = "l", col = adjustcolor("green3", alpha = 0.5), lwd = 3)
        points(smoo + abs.mean, type = "l", col = "blue", lwd = 2)
        points(smoo - abs.mean, type = "l", col = "blue", lwd = 2)
        
        points(smoo - abs.mean - (k * abs.sd), type = "l", col = "purple", lwd = 2)
        points(smoo + abs.mean + (k * abs.sd), type = "l", col = "purple", lwd = 2)
        
        points(cbind(bp, transect[bp], ncol = 2), col = "red")
        
        if (legend) {
            legend(legend.pos, bty = "n", lwd = 2,
                   col = c("green3", "blue", "purple", "red"), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 1), 
                   legend = c("Smoothed value", "Smoothed + mean amplitude", 
                              paste0("Smoothed + mean amplitude + ", round(k, 2), " SD"),
                              "Bad pixels"))
        }
    }
    
    list(smoothed = smoo, residual = res, flagged = bp)
} 



#' Identify extreme values along transect
#' 
#' Find extreme pixel values in single vector
#' @details Uses Loess smoothing to find 'spine' of vector and identifies any points that fall more than a certain constant multiple of the NMAD from the central spine.
#' @param transect Vector of pixel values to be assessed for 'bad' pixels
#' @param loess.span Filter span to be used in smoothing the data. Default is 1/15.
#' @param keep Central proportion of points to keep. Default is 0.95, which results in a cutoff for 'bad' pixels of 1.96sd.
#' @param display Boolean: plot the resulting transect? Default is T.
#' @param legend Boolean: include legend in plot? Default is F.
#' @param legend.pos String: where to display legend? default is "topright".
#' @param ... Additional parameters to be passed to plotting function
#' @return A list containing named vectors of Loess-smoothed values, residuals, and indices of pixels flagged.
#' @export
#' @examples 
#' zz <- bp.by.res.sd(pw.m[992, 1:992], display = T)
bp.by.res.nmad <- function(transect, loess.span = 1/15, keep = 0.95, display = T, legend = F, legend.pos = "topright",  ...) {
    
    k <- qnorm(1-(1-keep)/2, 0, 1)
    
    smoo <- lowess(transect, f = loess.span)$y
    res <- transect - smoo
    
    abs.mean <- mean(abs(res))
    abs.mad <- mad(abs(res))
    
    bp <- which(abs(res) > (abs.mean + k * abs.mad))
    
    if (display) {
        plot(transect, type = "l", col = adjustcolor("grey", alpha = 0.5), ...)
        points(transect, pch = 20, cex = 0.7)
        
        points(smoo, type = "l", col = adjustcolor("green3", alpha = 0.5), lwd = 3)
        
        points(smoo + abs.mean, type = "l", col = "blue", lwd = 2)
        points(smoo - abs.mean, type = "l", col = "blue", lwd = 2)
        
        points(smoo - abs.mean - (k * abs.mad), type = "l", col = "purple", lwd = 2)
        points(smoo + abs.mean + (k * abs.mad), type = "l", col = "purple", lwd = 2)
        
        points(cbind(bp, transect[bp], ncol = 2), col = "red")
        
        if (legend) {
            legend(legend.pos, bty = "n", lwd = 2,
                   col = c("green3", "blue", "purple", "red"), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 1), 
                   legend = c("Smoothed value", "Smoothed + mean amplitude", 
                              paste0("Smoothed + mean amplitude + ", round(k, 2), " NMAD"),
                              "Bad pixels"))
        }
    }
    
    list(smoothed = smoo, residual = res, flagged = bp)
} 


#' Identify extreme values in image
#' 
#' Find columnwise extreme values in image array
#' @details For each column of each panel, use Loess smoothing to find central spline and calculate residuals. Any points with a residual above a certain threshold are flagged as bad pixels.
#' @param image Single-layer 1996x1996 image array
#' @param method String: either "SD" or "NMAD". Threshold will be set according to multiples of SD or NMAD of columnwise residuals.
#' @param n.sigma Numeric scalar: how many multiples of the SD/NMAD will be used as threshold for 'extreme' values.
#' @param span Numeric value: smoothing span to be used in Loess spline
#' @param display Boolean: plot the results?
#' @param panels Boolean: if plotting the results, show subpanel edges?
#' @export
#' 
get.bad.pixels <- function(image, method = "sd", n.sigma = 1.96, span = 1/15, display = T, panels = T) {
    
    # subfunction to return indices of extreme values 
    bad.pixels <- function(transect, n.sigma, span, method) {
        smoo <- lowess(transect, f = span)$y
        res <- transect - smoo
        
        abs.mean <- mean(abs(res))
        if (method %in% c("mad", "nmad", "MAD", "NMAD")) {
            bound <- mad(abs(res))
        } else {
            bound <- sd(abs(res))
        }
        
        which(abs(res) > (abs.mean + n.sigma * bound))
    }
    
    # get coordinates of selected pixels
    lower.out <- apply(image[, 1:992], 1, bad.pixels, n.sigma, span, method)
    
    # label each list of pixel indices with column ref
    names(lower.out) <- c(1:1996)
    
    # remove any empty elements from list
    lower.out <- Filter(length, lower.out)
    
    # convert per-column indices to coordinates
    lower.coords <- do.call("rbind", Map(function(lower.out, i) cbind(col = as.integer(i), row = lower.out), lower.out, names(lower.out)))
    
    # repeat for upper panel
    upper.out <- apply(image[, 993:1996], 1, bad.pixels, n.sigma, span, method)
    names(upper.out) <- c(1:1996)
    upper.out <- Filter(length, upper.out)
    upper.coords <- do.call("rbind", Map(function(upper.out, i) cbind(col = as.integer(i), row = upper.out + 992), upper.out, names(upper.out)))
    
    coords <- rbind(lower.coords, upper.coords)
    
    if (display) {
        plot(coords, pch = 20, col = adjustcolor("blue", alpha = 0.3), 
             main = paste0("Bad pixels: ", n.sigma, method, " from smoothed column spline"))
        if (panels) draw.panels()
    }
    
    list(coords = coords, method = method, n.sigma = n.sigma, span = span)
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



#' Linear regression of observed pixel values
#' 
#' Perform linear regression of specified variables
#' @param im Acquisition array on which regression is to be performed. Must contain black, white and grey images.
#' @param terms Text string providing formula to be fitted. Default is "g ~ b * w", predicting grey value based on black and white.
#' @param res.only Boolean: return only a matrix of residuals (T) or (F) a list containing the fitted data, residuals and original data, along with the RMSE and adjusted RMSE.
#' @return If res.only == T, returns a matrix of residuals. If F, returns a list containing a data.frame and the RMSE of the fitted model.
#' @export
#' 
fit.w.lm <- function(im, terms = "g ~ b * w", midline = 1024.5, res.only = T) {
    
    df <- setNames(data.frame(melt(im[, , "black"]), 
                              melt(im[, , "grey"]), 
                              melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    if (!is.na(midline)) {
        df$upper <- df$y > midline
        terms <- paste0(terms, " + upper")
    }
    
    # fit linear model to central part of image only (excludes edge effects)
    w.lm <- lm(as.formula(terms), 
               data = df[findInterval(df$x, c(40.5, 2008.5)) == 1 & 
                             findInterval(df$y, c(40.5, 2008.5)) == 1, ])
    
    df$fv <- predict(w.lm, df)
    
    target <- gsub(" ~.*$", "", terms)
    
    df$res <- eval(parse(text = paste0("df$", target))) - df$fv
    
    if (res.only) {
        return(array(df$res, dim = dim(im[,,"black"])))
    } else {
        return(list(df = df, r2 = round(summary(w.lm)$adj.r.squared, 3), rmse = round(summary(w.lm)$sigma, 2)))
    }
}


#' Scatter plot of linear regression
#' 
#' @export
#' 
wlm.plot <- function(dt) {
    
    dt <- toString(dt)
    
    wlm <- eval(parse(text = paste0("lm.", dt)))
    
    pdf(paste0(fpath, "white-fit-", dt, ".pdf"))
    par(mar = c(4, 4, 1, 1))
    smoothScatter(wlm$df$w, wlm$df$fv, nrpoints = 0, xlim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed", ylab = "Fitted value")
    abline(line(wlm$df$w, wlm$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
    dev.off()
    
    write(paste0("Adj $r^2$ ", round(wlm$r2, 3), "; ",
                 "RMSE ", round(wlm$rmse, 2)),
          paste0(fpath, "fitted-wv-all-", dt, ".txt"))
}
