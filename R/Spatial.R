
#' convert coefficients to matrix
#'
#' Create a matrix for use in contour plotting, based on a PPM (point process model)
#' @export
#' 
ppmfit.matrix <- function(ppm) {
    
    zz <- setNames(melt(array(dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))), nm = c("x", "y", "z"))
    cc <- coef(ppm)
    
    zz$z <- cc[1] + (zz$x * cc["x"]) + (zz$y * cc["y"]) + 
        (zz$x^2 * cc["I(x^2)"]) + (zz$y^2 * cc["I(y^2)"]) + (zz$x * zz$y * cc["I(x * y)"]) 
    
    array(zz$z, dim = c(ppm$Q$data$window$xrange[2], ppm$Q$data$window$yrange[2]))
}