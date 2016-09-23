
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