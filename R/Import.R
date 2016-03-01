
#' Load one day's acquisitions from one channel
#'
#' Import one day's images from a single channel into an array. If the function is called by the user, the data is loaded into an object named automatically by the function. If called by another function, will return an array.
#' @param img.date Date of images to import, in format yymmdd.
#' @param channel Specify channel to import: black, grey or white.
#' @param fpath Path to top level of stored images. Default is "./Image-data/"
#' @param x Width of image. Default is 1996.
#' @param y Height of image. Default is 1996.
#' @param z Number of images to import. Default is 20.
#' @export
#' @examples
#' load.images(150828, "black")   # will create array b.150828
#' 
#' 
load.images <- function(img.date, channel, fpath = "./Image-data/", x = 1996, y = 1996, z = 20) {
    
    img.date <- toString(img.date)
    obj.nm <- paste0(substring(channel,1,1), ".",img.date)
    global = F
    
    # was function called directly or within another function?
    if (environmentName(parent.frame()) == "R_GlobalEnv") {
        global = T
        
        # check if object already exists: don't overwrite
        if (exists(obj.nm)) {
            err <- paste("An object named",obj.nm,"already exists. If you want to reload the data, please delete it first.")
            stop(err)
        }            
    }
    
    # convert date to filenm format, set up array
    rev.date <- paste0(substr(img.date,5,6),substr(img.date,3,4),substr(img.date,1,2))
    m <- array(dim = c(x, y, z))
    
    # read each file in turn, assign values into array
    # DON'T iterate over file names - images will not be in the correct order
    for (i in 1:z) {
        filenm <- paste0(fpath,img.date,"/",channel,"/",channel,"_",rev.date,"_",i,".tif")
        tmp <- readTIFF(filenm, as.is = T)
        
        # if imported array has more than 2 dimensions, keep only the first
        if (length(dim(tmp)) == 2) {
            tmp <- tmp
        } else {
            tmp <- tmp[,,rep(1, length(dim(tmp)) - 2)]
        }
        
        # flip & transpose matrix to match orientation of original image
        m[,,i] <- t(tmp[nrow(tmp):1,,drop=FALSE])
    } 
    
    # if called from global evironment, automatically create & name object
    # otherwise, return object so that it can be assigned by calling function
    if (global) {
        assign(obj.nm, m, envir = .GlobalEnv)
    } else {
        return(m)
    }
}


#' Load one day's acquisitions from all channels
#'
#' Import one day's images from into  1996x1996x20x3 array.
#' @param img.date Date of images to import, in format yymmdd.
#' @param fpath Path to top level of stored images. Default is "./Image-data/"
#' @export
#' @examples
#' img.150828 <- load.daily(150828)
#' 
#' 
load.daily <- function(img.date, fpath = "./Image-data/") {
    m <- array(dim = c(1996, 1996, 20, 3))
    m[,,,1] <- load.images(img.date, "black")
    m[,,,2] <- load.images(img.date, "grey")
    m[,,,3] <- load.images(img.date, "white")
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
#' b.150828.xml <- load.profiles(150828, "black")
#' 
#' 
load.profiles <- function(img.date, channel, fpath = "./Image-data/") {
    
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
    
    return(xml.data)
}



#' Summarise image profile data
#'
#' Searches through folder structure, extracts all .xml profiles, and summarises into a data frame.
#' @details If NA appears in data frame, the acquisitions for that channel do not all share the same value of that attribute, and should be investigated individually.
#' @param fpath Path to top level of image directory. Default is "./Image-data/"
#' @return Data frame containing summaries of image profile data.
#' @export
#' @examples
#' df <- summarise.profiles()
#' 
#' 
summarise.profiles <- function(fpath = "./Image-data/") {
    
    # get list of folders in given path location
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # check that valid path has been given
    if (length(folders) == 0) {
        cat("Image folder",fpath,"not found in",getwd())
    } else {
        
        # create df to hold summary
        df <- data.frame(acq.date = character(),
                         acq.time = character(),
                         channel = character(),
                         x.offset = double(),
                         y.offset = double(),
                         x.dim = double(),
                         y.dim = double(),
                         frames = double(),
                         stringsAsFactors = F)
        
        d <- length(folders)
        
        # read each folder in turn
        for (i in 1:d) {
            
            channels <- tolower(list.dirs(paste(fpath, folders[i],"/", sep = ""), full.names = F, recursive = F))
            
            # read each channel in turn
            for (j in 1:length(channels)) {
                
                r <- nrow(df) + 1
                profile <- load.profiles(folders[i], channels[j], fpath)
                
                df[r, c(1:3)] <- c(folders[i], profile$acq.time[1], channels[j])
                
                df$frames[r] <- nrow(profile)
                
                # check that offset & dimensions are same in all images
                if (max(profile$x.offset) == min(profile$x.offset)) {
                    df$x.offset[r] <- min(profile$x.offset)}
                
                if (max(profile$y.offset) == min(profile$y.offset)) {
                    df$y.offset[r] <- min(profile$y.offset)}
                
                if (max(profile$x.dim) == min(profile$x.dim)) {
                    df$x.dim[r] <- min(profile$x.dim)}
                
                if (max(profile$y.dim) == min(profile$y.dim)) {
                    df$y.dim[r] <- min(profile$y.dim)}
            }
        }
        return(df)
    }
}



#' Calculate summary statistics for all images
#'
#' Searches through folder structure, extracts all .tif images, and stores the summary statistics for each in a data frame.
#' @param fpath Path to top level of image directory. Default is "./Image-data/"
#' @return Data frame containing summaries of pixel values: mean and standard deviation, max, min, quartiles and IQR.
#' @export
#' @examples
#' df <- summarise.images()
#' 
#' 
summarise.images <- function(fpath = "./Image-data/") {
    
    df <- data.frame(img.date = character(),
                     channel = character(),
                     
                     mean = numeric(),
                     sd = numeric(),
                     
                     min = double(),
                     lq = double(),
                     median = double(),
                     uq = double(),
                     max = double(),
                     iqr = double(),
                     stringsAsFactors = F)                    
    
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # read each folder in turn
    for (i in 1:length(folders)) {
        
        channels <- tolower(list.dirs(paste(fpath, folders[i],"/", sep = ""), full.names = F, recursive = F))
        
        # read each channel in turn
        for (j in 1:length(channels)) {
            
            r <- nrow(df) + 1
            data <- load.images(folders[i], channels[j], fpath)
            summ <- batch.summary(data)
            
            df[r,c(1:2)] <- c(folders[i], channels[j])
            df[r,c(3:4)] <- summ[[1]]
            df[r,c(5:10)] <- summ[[2]]
        }
    }
    return(df)
}



#' Count available images
#'
#' Searches through folder structure to quickly summarise number of available .tif images from each channel on each day. Assumes files are organised as /Image-data/date/channel/image.tif
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



