
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
    z.dist <- data.frame(x = c(1:1996^2) %/% 1996 + 1, y = c(1:1996^2) %% 1996)
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


#' Per-panel linear regression
#' 
#' Fit a linear gradient over each of the 32 subpanels within the image
#' @param image Single-layer array image (1996x1996) to be fitted.
#' @param terms String specifying terms to be fitted to each subpanel. Default is "x + y", fitting a linear model to the x and y coordinates of each subpanel without considering interactions.
#' @return List containing the terms applied, a matrix of fitted values, and a matrix of model coefficients for each panel.
#' @param robust Boolean: use robust model fitting (\link{\code{rlm}})? Default is F (use \link{\code{lm}} for regression).
#' @export
#' @examples
#' panel.lm <- panel.lm(circ.res, "poly(x, 2) + y")    # finds coefficients of intercept, x, x^2, y.
#' panel.lm <- panel.lm(circ.res, "y")                 # finds coefficients of intercept and y.
#' panel.res <- circ.res - panel.lm$fitted.values
#' pixel.image(panel.res)
panel.lm <- function (image, terms = "x + y", robust = F) {
    
    # convert formula string to lower case to ensure match
    terms <- tolower(terms)
    
    # create empty arrays to hold output
    coeffs <- c()
    smoothed.panels <- array(dim = c(128, 1024, 32))
    
    # iterate over subpanels, fit specified linear model to each
    sp <- subpanels(image)
    
    for (s in 1:32) {
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
    list(formula = paste0("value ~ ", terms), fitted.values = join.panels(smoothed.panels), models = coeffs)
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