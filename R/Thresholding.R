
# THRESHOLDING & CLASSIFICATION OF VALUES                                                       ####

#' Modal density
#' 
#' Calculate location of modal density of a distribution
#' @param dat Data for which bounds are to be calculated
#' @return Value at which modal density occurs
#' @export
#' 
modal.density <- function(dat) {
    zz <- density(dat, n = 65536, na.rm = T)
    zz$x[which.max(zz$y)]
}



#' Asymmetric NMAD thresholds
#' 
#' Calculate modal density mu of data and find boundaries at distance n * MAD either side of mu.
#' @param dat Data for which bounds are to be calculated
#' @param n Number of multiples of MAD to be applied either side of mu
#' @return Vector containing upper and lower bounds
#' @export
#' 
asymmetric.mad <- function(dat, n = 6) {
    
    mu <- modal.density(dat)

    # calculate MAD on either side of maximum density, centred at maximum density
    c(mu - n * mad(dat[which(dat <= mu)], center = mu),
      mu + n * mad(dat[which(dat > mu)], center = mu))
}


#' Calculate class boundaries
#' 
#' For a pixelwise mean image, calculate thresholds for each classification
#' @param dat Data for which boundaries are to be calculated
#' @param b Boundaries to calculate (as well as \code{\link{asymmetric.mad}}): default is c(0.25, 0.5)
#' @return Vector of cutpoints
#' @export
#' 
class.boundaries <- function(dat, b = c(0.25, 0.5)) {
    
    rng <- c(floor(min(dat, na.rm = T)), ceiling(max(dat, na.rm = T)) + 1)
    inner.bounds <- asymmetric.mad(dat)
    
    vb <- inner.bounds[2] + ((rng[2] - inner.bounds[2]) * b)
    
    sort(c(rng, inner.bounds, vb))
}


#' Label extreme-valued pixels
#'
#' Full explanation still required
#' @export
#' 
classify.px <- function(im, levels = c("dim", NA, "s.bright", "bright", "v.bright")) {
    
        class <- array(findInterval(im, class.boundaries(im)), dim = dim(im))
        class[class == 2] <- NA
        
        if(length(levels) < 5) {
            levels <- c(levels, rep(levels[length(levels)], 5))[1:5]
        }

        data.frame(which(!is.na(class), arr.ind = T),
                   type = ordered(levels[class[which(!is.na(class))]],
                                  levels = rev(unique(levels))))
}

####################################################################################################

# FITTING A LINEAR REGRESSION                                                                   ####

#' Linear regression of observed pixel values
#' 
#' Support function: perform linear regression of specified variables
#' @param im Acquisition array on which regression is to be performed. Must contain black, white and grey images.
#' @param terms Text string providing formula to be fitted. Default is "g ~ b * w", predicting grey value based on black and white.
#' @param res.only Boolean: return only a matrix of residuals (T) or (F) a list containing the fitted data, residuals and original data, along with the RMSE and adjusted RMSE.
#' @return If res.only == T, returns a matrix of residuals. If F, returns a list containing a data.frame and the RMSE of the fitted model.
#' @export
#' 
fit.gv.lm <- function(im, terms = "g ~ b * w", midline = 1024.5, res.only = T) {
    
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

####################################################################################################

# FITTING A BIVARIATE GAUSSIAN SPOT                                                             ####

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

####################################################################################################

