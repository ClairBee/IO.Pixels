
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
