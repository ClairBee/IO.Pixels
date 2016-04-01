
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
    z.dist <- cbind(x = c(1:1996^2) %/% 1996 + 1, y = c(1:1996^2) %% 1996)
    z.dist$z <- sqrt((z.dist$x - centre[1])^2 + (z.dist$y - centre[2])^2)
    
    # distance from centre vs pixel values
    zz <- cbind(melt(image), z.dist)[,c(4,5,6,3)]
    
    # fit & return linear model
    return(lm(value ~ z, zz))
}


#' Linear regression over distance from centre, including polynomial terms of z
#' 
#' Fit a circular spot model with linear gradients to the image, including polynomial terms of z as explanatory variables
#' @param image Single-layer (1996x1996) array image to be fitted.
#' @param o Order of polynomial to include. Default is 2
#' @param centre Vector of coordinates of centre of spot. Default is centre of uncropped detector, c(1023.5, 992.5)
#' @return Linear model fitted to given image
#' @export
#' @examples
#' circ.poly.lm <- fit.circular.lm.poly(pw.m)
#' circ.poly.res <- matrix(circ.poly.lm$residuals, ncol = 1996)
#' pixel.image(circ.poly.res)
#' 
fit.circular.lm.poly <- function(image, o = 2, centre = c(1023.5, 992.5)) {
    
    # get distances from spot centre
    z.dist <- cbind(x = c(1:1996^2) %/% 1996 + 1, y = c(1:1996^2) %% 1996)
    z.dist$z <- sqrt((z.dist$x - centre[1])^2 + (z.dist$y - centre[2])^2)
    
    # distance from centre vs pixel values
    zz <- cbind(melt(image), z.dist)[,c(4,5,6,3)]
    
    # fit & return linear model
    return(lm(value ~ poly(z, o), zz))
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


#' Per-panel linear regression, with polynomial terms
#' 
#' Fit a linear gradient over each of the 32 subpanels within the image, including polynomial terms of x and y.
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @param o Order of polynomial to include. Default is 2
#' @return List containing a matrix of fitted values, and a matrix of model coefficients for each panel.
#' @export
#' @examples
#' panel.lm.2 <- fit.panel.lm.poly(circ.res.2, o = 2)
#' panel.res.2 <- circ.res.2 - panel.lm.2$fitted.values
#' pixel.image(panel.res.2)
fit.panel.lm.poly <- function(image, o = 2) {
    sp <- subpanels(image)
    
    coeffs <- array(dim = c(32, 5),
                    dimnames = list(dimnames(sp)[[3]], c("offset", "x","x^2", "y", "y^2")))
    smoothed.panels <- array(dim = c(128, 1024, 32))
    
    for (s in 1:32) {
        lm <- lm(value ~ poly(X1,o) + poly(X2,o), data = melt(sp[ , , s]))
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