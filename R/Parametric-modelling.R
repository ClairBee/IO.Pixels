
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