
#' Load Daily Images
#'
#' Import one day's images into a single array
#' @param img.date Date of images to import, in format yymmdd.
#' @param channel Specify channel to import: black, grey or white.
#' @param x Width of image. Default is 1996.
#' @param y Height of image. Default is 1996.
#' @param z Number of images to import. Default is 20.
#' @param progress Display progress bar? T or F (default is F)
#' @param path Path to top level of stored images. Default is "./Image-data/"
#' @export
#' @examples
#' daily <- load.daily(150828, "black")
#' 
#' 
load.daily <- function(img.date, channel, x = 1996, y = 1996, z = 20, progress = F, path = "./Image-data/") {
    
    # convert date to filenm format, set up array
    img.date <- toString(img.date)
    rev.date <- paste(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2), sep = "")
    m <- array(dim = c(x, y, z))
    
    # show progress bar or not?
    if (progress) {pb <- txtProgressBar(max = z, style = 3)}
    
    # read each file in turn, assign values into array
    for (i in 1:z) {
        filenm <- paste(path,img.date,"/",channel,"/",channel,"_",rev.date,"_",i,".tif", sep = "")
    
        m[,,i] <- readTIFF(filenm)
        if (progress) {setTxtProgressBar(pb,i)}
    } 
    return(m)
}
