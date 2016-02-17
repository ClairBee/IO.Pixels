#' Load one day's acquisitions from one channel
#'
#' Import one day's images into a single array
#' @param img.date Date of images to import, in format yymmdd.
#' @param channel Specify channel to import: black, grey or white.
#' @param x Width of image. Default is 1996.
#' @param y Height of image. Default is 1996.
#' @param z Number of images to import. Default is 20.
#' @param progress Display progress bar? T or F (default is F)
#' @param fpath Path to top level of stored images. Default is "./Image-data/"
#' @export
#' @examples
#' b.150828 <- load.image.batch(150828, "black")
#' 
#' 
load.image.batch <- function(img.date, channel, fpath = "./Image-data/") {
    
    # convert date to filenm format, set up array
    img.date <- toString(img.date)
    rev.date <- paste(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2), sep = "")
    
    # check if xml summ already in df? Currently always overwrites
    xml.data <- load.xml.batch(img.date, channel, fpath)
    x <- min(xml.data$x.dim)
    y <- min(xml.data$y.dim)
    z <- nrow(xml.data)
            
    fpath <- paste(fpath, img.date, "/", channel, "/", sep = "")
    files <- list.files(fpath, pattern = "\\.tif$")
    
    m <- array(dim = c(x, y, z))
    
    # read each file in turn, assign values into array
    for (i in 1:length(files)) {
        filenm <- paste(fpath,channel,"_",rev.date,"_",i,".tif", sep = "")
        tmp <- readTIFF(filenm)
        if (length(dim(tmp)) > 2) {
            tmp <- tmp[,,rep(1,length(dim(tmp))-2)]
        }

        # flip & transpose matrix to match orientation of original image
        m[,,i] <- t(tmp[nrow(tmp):1,,drop=FALSE])
    }
    
    # overwrite df containing file summary
    store.batch.data(img.date, channel, xml.data, m)
    return(m)
}


#' Load xml profiles of one day's acquisitions from one channel
#'
#' Support function: import xml profiles of one day's images into a single data frame
#' @param img.date Date of images to import, in format yymmdd.
#' @param channel Specify channel to import: black, grey or white.
#' @param fpath Path to top level of stored images. Default is "./Image-data/"
#' @export
#' @examples
#' b.150828.xml <- load.xml.batch(150828, "black")
#' 
#' 
load.xml.batch <- function(img.date, channel, fpath = "./Image-data/") {
    
    img.date <- toString(img.date)
    
    fpath <- paste(fpath, img.date, "/", channel, "/", sep = "")
    files <- list.files(fpath, pattern = "\\.xml$")
    
    xml.data <- data.frame(acq.time = character(),
                           x.offset = double(),
                           y.offset = double(),
                           x.dim = double(),
                           y.dim = double(),
                           stringsAsFactors = F)
    
    for (i in 1:length(files)) {
        tif.profile <- xmlToList(xmlParse(paste(fpath, files[i], sep = "")))
        xml.data[i,1] <- tif.profile$Time
        xml.data[i,c(2:5)] <- c(as.integer(tif.profile$CameraProperties$imageOffsetX),
                                as.integer(tif.profile$CameraProperties$imageOffsetY),
                                as.integer(tif.profile$CameraProperties$imageSizeX),
                                as.integer(tif.profile$CameraProperties$imageSizeY))
    }
    xml.data <- xml.data[order(xml.data$acq.time),]
    row.names(xml.data) <- c(1:length(files))
    return(xml.data)
}


#' Store summary of image batch
#'
#' Support function: store key summaries of each day's images in a data frame
#' @param img.date Date of images to import, in format yymmdd.
#' @param channel Specify channel to import: black, grey or white.
#' @param xml.data Output from \link{load.xml.batch} containing tiff profile summaries
#' @param image.data Matrix of imported tiff images
#' @export
#' @examples
#' store.batch.data(150828, "black", b.150828.xml, b.150828)
#' 
#' 
store.batch.data <- function(img.date, channel, xml.data, image.data) {
    
    img.date <- toString(img.date)
    obj.nm <- paste(substr(channel,1,1),".", img.date, sep = "") 
    
    # create df to hold summary, if not already existing
    if (!exists("img.summ")) {
        img.summ <<- data.frame(obj.nm = character(),
                                acq.date = character(),
                                acq.time = character(),
                                channel = character(),
                                mean = numeric(),
                                sd = numeric(),
                                x.offset = double(),
                                y.offset = double(),
                                x.dim = double(),
                                y.dim = double(),
                                frames = double(),
                                stringsAsFactors = F)
        r <- 1                                                     # first row
    } else {
        if (obj.nm %in% row.names(img.summ)) {
            r <- which(row.names(img.summ) == obj.nm, arr.ind = T)     # overwrite existing
        } else {
            r <- nrow(img.summ) + 1                                # next available row
        }
    }
    
    # update row r of df with summary information
    
    img.summ[r, 1:4] <<- c(obj.nm, img.date, xml.data$acq.time[1], channel)
    
    img.summ$mean[r] <<- mean(image.data)
    img.summ$sd[r] <<- sd(image.data)
    
    img.summ$frames[r] <<- nrow(xml.data)
    
    # check that offset & dimensions are same in all images
    if (max(xml.data$x.offset) == min(xml.data$x.offset)) {
        img.summ$x.offset[r] <<- min(xml.data$x.offset)}
    
    if (max(xml.data$y.offset) == min(xml.data$y.offset)) {
        img.summ$y.offset[r] <<- min(xml.data$y.offset)}
    
    if (max(xml.data$x.dim) == min(xml.data$x.dim)) {
        img.summ$x.dim[r] <<- min(xml.data$x.dim)}
    
    if (max(xml.data$y.dim) == min(xml.data$y.dim)) {
        img.summ$y.dim[r] <<- min(xml.data$y.dim)}
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



