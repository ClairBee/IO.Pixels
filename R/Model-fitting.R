
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
#' @return List containing the terms applied, a matrix of fitted values, and a matrix of model coefficients for each panel.
#' @export
#' @examples
#' panel.lm <- panel.lm(circ.res, "poly(x, 2) + y")    # finds coefficients of intercept, x, x^2, y.
#' panel.lm <- panel.lm(circ.res, "y")                 # finds coefficients of intercept and y.
#' panel.res <- circ.res - panel.lm$fitted.values
#' pixel.image(panel.res)
panel.lm <- function (image, terms = "x + y", robust = F, left.pad = 2, upper.pad = 20) {
    
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