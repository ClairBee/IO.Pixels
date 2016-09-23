
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