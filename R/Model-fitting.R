
#' Linear regression over distance from centre
#' 
#' Fit a circular spot model with linear gradients to the image.
#' @param image Single-layer (1996x1996) array image to be fitted.
#' @param centre Vector of coordinates of centre of spot. Default is centre of uncropped detector, c(1023.5, 992.5)
#' @return Linear model fitted to given image
#' @export
#' @examples
#' circ.lm <- fit.circular.lm(pw.m)
#' circ.res <- matrix(circ.lm$residuals, ncol = 1996)
#' pixel.image(circ.res)
fit.circular.lm <- function(image, centre = c(1023.5, 992.5)) {
    
    # get distances from spot centre
    z.dist <- merge(x = c(1:dim(image)[1]), y = c(1:dim(image)[2]))
    z.dist$z <- sqrt((z.dist$x - centre[1])^2 + (z.dist$y - centre[2])^2)
    
    # distance from centre vs pixel values
    zz <- cbind(melt(image), z.dist)[,c(4,5,6,3)]
    
    # fit & return linear model
    return(lm(value ~ z, zz))
}


#' Per-panel linear regression
#' 
#' Fit a linear gradient over each of the 32 subpanels within the image
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @return List containing a matrix of fitted values, and a matrix of model coefficients for each panel.
#' @export
#' @examples
#' panel.lm <- fit.panel.lm(circ.res)
#' panel.res <- circ.res - panel.lm$fitted.values
#' pixel.image(panel.res)
fit.panel.lm <- function(image) {
    sp <- subpanels(image)
    
    coeffs <- array(dim = c(32, 3),
                    dimnames = list(dimnames(sp)[[3]], c("offset", "x", "y")))
    smoothed.panels <- array(dim = c(128, 1024, 32))
    
    for (s in 1:32) {
        lm <- lm(value ~ X1 + X2, data = melt(sp[ , , s]))
        coeffs[s,] <- coef(lm)
        smoothed.panels[,,s] <- predict(lm, melt(sp[ , , s]))
    }
    
    list(fitted.values = join.panels(smoothed.panels), models = coeffs)
}


#' Per-column loess smoothing
#' 
#' Fit a loess spline to each column of the image using \link{\code{lowess}}. The upper and lower detector panels are treated separately.
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @param span Span to be used for loess smoothing. Default is 1/15.
#' @return 1996x1996 matrix of smoothed values.
#' @export
#' @examples
#' smoothed <- lowess.per.column(panel.res)
#' smoothed.res <- panel.res - smoothed
#' pixel.image(smoothed.res)
#' s.hist(smoothed.res)
lowess.per.column <- function(image, span = 1/15) {
    
    # split image into upper & lower panels: smoothing across midline is nonsensical
    upper <- image[, 993:1996]
    lower <- image[, 1:992]
    
    # apply lowess smoothing per column to each panel in turn
    upper.smoo <- do.call("rbind", lapply(apply(upper, 1, lowess, f = span), "[[", 2))
    lower.smoo <- do.call("rbind", lapply(apply(lower, 1, lowess, f = span), "[[", 2))
    
    # join upper & lower panels again 
    smoo <- array(dim = c(1996, 1996))
    smoo[, 993:1996] <- upper.smoo
    smoo[, 1:992] <- lower.smoo
    
    return(smoo)
}