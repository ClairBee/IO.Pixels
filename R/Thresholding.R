

#' Load one day's acquisitions from one channel
#'
#' Import one day's images from a single channel into an array. If the function is called by the user, the data is loaded into an object named automatically by the function. If called by another function, will return an array.
#' @param data Array containing a sequence of images from a single channel
#' @param too.bright Numeric value to use as threshold for underperforming bright pixels: the multiple of median value above which a pixel is deemed to be unacceptably bright. Default is 1.5 x median.
#' @param too.dim Numeric value to use as threshold for underperforming dark pixels: the multiple of median value below which a pixel is deemed to be unacceptably bright. Default is 0.45 x median.
#' @param noisy Numeric value to use as threshold for noisy pixels: the number of median pixel SDs above which a pixel is deemed to be noisy. Default is 6.
#' @return Data frame containing coordinates and classifications of all underperforming pixels in the given array.
#' @export
#' @examples
#' bp.b.150828 <- threshold.pixels(b.150828)
#' 
#' 
threshold.pixels <- function(data, too.bright = 1.5, too.dim = 0.45, noisy = 6) {
    
    # create temporary list to hold output
    bp <- list()
    
    pw.mean <- pixelwise.mean(data)
    pw.sd <- pixelwise.sd(data)
    median.value <- median(data)
    median.sd <- median(pw.sd)
    
    
    # always use tryNULL when creating data frames, in case no rows are returned
    
    # order in which each type is added to df is order of priority when classifying
    # (eg. duplicates in later categories will be removed)
    
    # dead pixels: daily mean value is exactly 0
    bp$dead <- tryNULL(cbind(data.frame(which(pw.mean == 0.0, arr.ind = T)), 
                             type = "Dead"))
    
    # hot pixels: daily mean value is exactly 65535
    bp$hot <- tryNULL(cbind(data.frame(which(pw.mean == 65535.0, arr.ind = T)),
                            type = "Hot"))
    
    # underperforming bright pixels: value > 1.5x median
    bp$bright <- tryNULL(cbind(data.frame(unique(which(data > (too.bright * median.value), arr.ind = T)[,c(1:2)])),
                               type = "Too bright"))
    if (!is.null(bp$bright)) {colnames(bp$bright) <- c("row", "col", "type")}
    
    # underperforming dark pixels: value < 0.45x median
    bp$dim <- tryNULL(cbind(data.frame(unique(which(data < (too.dim * median.value), arr.ind = T)[,c(1:2)])),
                            type = "Too dark"))
    if (!is.null(bp$dim)) {colnames(bp$dim) <- c("row", "col", "type")}
    
    # noise: pixel sigma > 6x median sigma
    bp$noisy <- tryNULL(cbind(data.frame(which(pw.sd > (noisy * median.sd), arr.ind = T)),
                              type = "Noisy"))
    
    # merge separate data frames into one
    df <- data.frame(do.call("rbind", bp))
    
    # remove duplicates and return df
    return(df[!duplicated(df[,c(1:2)]),])
}