
#' Convert points to raster
#' 
#' Support function: convert matrix of coordinates into a raster, so that convolution, clumping etc can be applied
#' @param px data.frame or matrix, the first two columns of which contain the X and Y coordinates of points to be marked on the raster
#' @param im.dim Vector of dimensions for the final raster. Default is c(2048, 2048).
#' @return Raster object with the coordinates marked as 1, everything else as NA.
#' @export
#' 
px2r <- function(px, im.dim = c(2048, 2048)) {
    
    mat <- array(0, dim = im.dim)
    mat[as.matrix(px[,1:2])] <- 1
    
    # rotate matrix by reversing and transposing, then convert to raster
    mat <- t(mat)[im.dim[2]:1,]
    raster(mat, xmn = 0.5, xmx = im.dim[1] + .5, ymn = 0.5, ymx = im.dim[2] + .5)
}


#' Convert matrix to raster
#' 
#' Convert an image array to a raster file, maintaining the correct orientation for plotting
#' @param im Single-layer image array to be converted to a raster
#' @return RasterLayer object containing the image values
#' @export
#' 
m2r <- function(im) {
    # rotate matrix by reversing and transposing, then convert to raster
    mat <- t(im)[dim(im)[2]:1,]
    raster(mat, xmn = 0.5, xmx = dim(im)[1] + .5, ymn = 0.5, ymx = dim(im)[2] + .5)
}


#' Convert raster to matrix
#' 
#' Convert a raster to an image array, maintaining the correct orientation for plotting
#' @param rast Single-layer RasterLayer object to be converted to a 2d array
#' @return Single-layer image array
#' @export
#' 
r2m <- function(rast) {
    mat <- matrix(getValues(rast), ncol = dim(rast)[2], byrow = T)
    t(mat[dim(mat)[1]:1,,drop = F])
}