

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

