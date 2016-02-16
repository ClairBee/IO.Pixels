#' Load one day's acquisitions from one channel
#'
#' Import one day's images from a single channel into an array
#' @param img.date Date of images to import, in format yymmdd.
#' @param channel Specify channel to import: black, grey or white.
#' @param x Width of image. Default is 1996.
#' @param y Height of image. Default is 1996.
#' @param z Number of images to import. Default is 20.
#' @param progress Display progress bar? T or F (default is F)
#' @param fpath Path to top level of stored images. Default is "./Image-data/"
#' @export
#' @examples
#' acq <- load.acq(150828, "black")
#' 
#' 
load.acq <- function(img.date, channel, x = 1996, y = 1996, z = 20, progress = F, fpath = "./Image-data/") {
    
    # convert date to filenm format, set up array
    img.date <- toString(img.date)
    rev.date <- paste(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2), sep = "")
    m <- array(dim = c(x, y, z))
    
    # show progress bar or not?
    if (progress) {pb <- txtProgressBar(max = z, style = 3)}
    
    # read each file in turn, assign values into array
    for (i in 1:z) {
        filenm <- paste(fpath,img.date,"/",channel,"/",channel,"_",rev.date,"_",i,".tif", sep = "")
        tmp <- readTIFF(filenm)
        if (length(dim(tmp)) == 2) {
            tmp <- tmp
        } else {
            tmp <- tmp[,,1]
        }
        # flip & transpose matrix to match orientation of original image
        m[,,i] <- t(tmp[nrow(tmp):1,,drop=FALSE])
        
        if (progress) {setTxtProgressBar(pb,i)}
    } 
    return(m)
}


#' Count available images
#'
#' Searches through folder structure to count number of available .tif images. Assumes files are organised as /Image-data/date/channel/image.tif
#' @param fpath Path to top level of stored images. Default is "./Image-data/"
#' @return Data frame containing counts of .tif images in channel subfolders for each date, with image acquisition dates as row names.
#' @export
#' @examples
#' image.count <- count.images()
#' 
#' 
count.images <- function(fpath = "./Image-data/") {
    
    # get list of folders in given path location
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # check that valid path has been given
    if (length(folders) == 0) {
        cat("Image folder",fpath,"not found in",getwd())
    } else {
        
        d <- length(folders)
        df <- data.frame(black = double(), grey = double(), white = double())
        
        for (i in 1:d) {    
            # folder name should be acquisition date (yymmdd)
            img.date <- folders[i]
            rev.date <- paste(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2), sep = "")
            
            # check folder for each channel in turn, get number of .tif files
            b <- length(list.files(paste(fpath, folders[i],"/black/", sep = ""), pattern = "\\.tif$"))
            g <- length(list.files(paste(fpath, folders[i],"/grey/", sep = ""), pattern = "\\.tif$"))
            w <- length(list.files(paste(fpath, folders[i],"/white/", sep = ""), pattern = "\\.tif$"))
            
            # add .tif count to data frame
            df[i,] <- c(b,g,w)
            row.names(df)[i] <- img.date
        }
        return(df)
    }
}


#' Get pixelwise mean of daily acquisitions
#'
#' Return pixelwise mean value across an array
#' @param data Array with 3 or more dimensions
#' @return Array containing pixelwise mean of all acquired values
#' @export
#' @examples
#' pw.mean <- pixelwise.mean(daily)
#' 
#' 
pixelwise.mean <- function(data) {
    # alter this to allow pixelwise mean over 4d array
    apply(data, c(1,2), mean)
}


#' Get pixelwise standard deviation of daily acquisitions
#'
#' Return pixelwise SD across an array
#' @param data Array with 3 or more dimensions
#' @return Array containing pixelwise standard deviation of all acquired values
#' @export
#' @examples
#' pw.sd <- pixelwise.sd(daily)
#' 
#' 
pixelwise.sd <- function(data) {
    # alter this to allow pixelwise sd over 4d array
    apply(data, c(1,2), sd)
}



