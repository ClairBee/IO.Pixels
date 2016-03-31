

#' Identify 'bad' pixels using quantiles of residuals after parametric fitting
#'
#' Load a specified image and fit parametric circular spot and panelwise gradients, then cut residuals at specified quantiles to classify 'bad' pixels
#' @details Requires that pixelwise mean files are already loaded, as by \link{\code{load.pixel.means}}
#' @param img.date Integer/string specifying image date to be imported
#' @param batch String specifying colour batch to be imported (black, grey or white)
#' @param circ.o Integer (or NA if no spot model is required) specifying order of polynomial model to fit to circular spot.
#' @param panel.o Integer (or NA if no per-panel model is required) specifying order of polynomial model to fit to panel X and Y coordinates.
#' @param lq Lower quartile at which to start looking for 'outliers'. Default is 0.001
#' @param uq Upper quartile at which to start looking for 'outliers'. Default is 0.999
#' @param dist Multiples of IQR away from lq/uq at which a point is classified as 'bad'.
#' @param details Boolean: return lists of identified 'bad pixels' (T) or only summary (F)? Default is F.
#' @export
#' @examples
#' load.pixel.means()
#' zz <- rbind(bp.parametric.quantiles(160314, "white"),
#'             bp.lowess.quantiles(160314, "white"))
#' 
bp.parametric.quantiles <- function(img.date, batch, circ.o = 2, panel.o = 1, 
                                    lq = 0.001, uq = 0.999, 
                                    dist = 1.5, details = F) {
    img.date <- toString(img.date)
    im <- pw.m[,, batch, img.date]
    
    
    # fit circular model (if specified)
    if (!is.na(circ.o)) {
        if (circ.o == 1) {
            circ.lm <- fit.circular.lm(im)
            c.par <- "linear spot"
        } else {
            circ.lm <- fit.circular.lm.poly(im, o = circ.o)
            c.par <- paste0("poly (", circ.o, ") spot")
        }
        circ.res <- matrix(circ.lm$residuals, ncol = 1996)
    } else {
        circ.res <- im
        c.par = "no spot"
    }
    
    # fit per-panel model (if specified)
    if (!is.na(panel.o)) {
        if (panel.o == 1) {
            panel.lm <- fit.panel.lm(circ.res)
            p.par <- "linear panels"
        } else {
            panel.lm <- fit.panel.lm.poly(circ.res, o = panel.o)
            p.par <- paste0("poly (", panel.o, ") panels")
        }
        res <- circ.res - panel.lm$fitted.values
    } else {
        res <- circ.res
        p.par <- "no panels"
    }
    
    # cutoff points for 'bad' pixels, based on residual quantiles
    ul <- quantile(res, uq) + (dist * IQR(res))
    ll <- quantile(res, lq) - (dist * IQR(res))
    
    summ <- data.frame(method = paste(c.par, p.par, "quantiles", sep = "; "),
                       img.date, batch,
                       params = paste(lq, uq, dist, sep = ", "), 
                       lower = ll, upper = ul, 
                       n.low = length(which(res < ll)),
                       n.high = length(which(res > ul)),
                       res.mean = mean(res), 
                       res.sd = sd(res), res.mad = mad(res), res.med = median(res), 
                       row.names = NULL)
    
    if (details) {
        return(list(details = summ, 
                    low = which(res < ll, arr.ind = T), 
                    high = which(res > ul, arr.ind = T)))
    } else {
        return(summ)
    }
}

