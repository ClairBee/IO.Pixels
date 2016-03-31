
#' Create sample of 'healthy' pixels
#'
#' Given a matrix of coordinates of bad pixels, create a same-size random sample of healthy pixels that do not appear on the bad pixel list.
#' @param bad.pixels Data frame or matrix whose first two columns are X and Y coordinates of identified bad pixels.
#' @param n Size of sample to create. Default is to create a sample of same size as bad pixel list.
#' @return Matrix containing X and Y coordinates of sampled pixels.
#' @export
#' @examples
#' bp.160314 <- reset.bp(160314)
#' samp.160314 <- sample.healthy(bp.160314)
#' 
#' 
sample.healthy <- function(bad.pixels, n = nrow(bad.pixels)) {
    
    # convert coords to single integer
    # subtract 1997 to because coords start at 1, not 9
    bp.int <- (bad.pixels[, 1] * 1996) + bad.pixels[, 2] - 1997

    # sample as integer (much more efficient than permuting)
    s <- sample(c(1:(1996^2))[-c(bp.int)], n)
    
    # decompose into integer division & modulo components to get x & y
    # add 1 to get coordinates to correct range
    cbind(x = s %/% 1996 + 1, y = s %% 1996 + 1)
}

