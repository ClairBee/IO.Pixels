

#' Create array of all pixelwise means
#'
#' Load all pixelwise mean summaries (already created and stored as .rds objects) into a single labelled array.
#' @details RDS objects for import can be downloaded from https://github.com/ClairBee/Pixels/tree/master/02_Objects/images
#' @param fpath Path to top level of stored .rds objects. Default is "./02_Objects/images"
#' @param acq.list Vector containing acquisitions to import. If left blank, will import every file in the target folder named "pwm-xxxxxx.rds".
#' @export
#' 
load.pixel.means <- function(fpath = "./02_Objects/images", acq.list) {
    
    # remove trailing blackslash, if necessary
    fpath <- gsub("/$", "", fpath)
    
    if (missing (acq.list)) {
        acq.list <- gsub(".rds", "", gsub("pwm-", "", list.files(fpath, pattern = "pwm-[a-z, A-Z, 0-9]+\\.rds$", full.names = F)))
    }
    abind(sapply(acq.list, function(nm) readRDS(paste0(fpath, "/pwm-", nm, ".rds")), 
                 simplify = F),
          along = 4)
}

#' Create array of all pixelwise SDs
#'
#' Load all pixelwise standard deviation matrices into a single labelled array.
#' @param acq.list Vector containing acquisitions to import. If left blank, will import every file in the target folder named "pwm-xxxxxx.rds".
#' @param fpath Path to top level of stored .rds objects. Default is "./02_Objects/images"
#' @export
#' 
load.pixel.sds <- function(acq.list, fpath = "./02_Objects/stdevs") {
    if (missing (acq.list)) {
        acq.list <- gsub(".rds", "", gsub("pwsd-", "", list.files(fpath, pattern = "pwsd-[a-z, A-Z, 0-9]+\\.rds$", full.names = F)))
    }
    abind(sapply(acq.list, function(nm) readRDS(paste0(fpath, "/pwsd-", nm, ".rds")), 
                 simplify = F),
          along = 4)
}


#' Create array from stored .rds objects
#'
#' Load all .rds ojects of a specified type into a single labelled array.
#' @details RDS objects for import can be downloaded from https://github.com/ClairBee/Pixels/tree/master/02_Objects
#' @param fpath Path to top level of stored .rds objects. Default is "./02_Objects/images"
#' @param otype String specifying which type of image should be imported. Default is "pwm" (pixelwise mean images).
#' @param acq.list Vector containing acquisitions to import. If left blank, will import every file in the target folder named "<otype>-xxxxxx.rds".
#' @export
#' 
load.objects <- function(fpath = "./02_Objects/images", otype = "pwm", acq.list) {
    
    # remove trailing blackslash, if necessary
    fpath <- gsub("/$", "", fpath)
    otype <- paste0(gsub("-$", "", otype), "-")
    
    if (missing (acq.list)) {
        acq.list <- gsub(".rds", "", gsub(otype, "", list.files(fpath, pattern = paste0(otype, "[a-z, A-Z, 0-9]+\\.rds$"), full.names = F)))
    }
    abind(sapply(acq.list, function(nm) readRDS(paste0(fpath, "/", otype, nm, ".rds")), 
                 simplify = F),
          rev.along = 0)
}

####################################################################################################

#' Load one day's acquisitions from one batch
#'
#' Import one day's images from a single batch into an array. If the function is called by the user, the data is loaded into an object named automatically by the function. If called by another function, will return an array.
#' @param img.date Date of images to import, in format yymmdd.
#' @param batch Specify batch to import: black, grey or white.
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @param x Width of image. Default is 1996.
#' @param y Height of image. Default is 1996.
#' @param z Number of images to import. Default is 20.
#' @export
#' @examples
#' load.images(150828, "black")   # will create array b.150828
#' 
#' 
load.images <- function(img.date, batch, fpath = "/home/clair/Documents/Pixels/Image-data/", x = 1996, y = 1996, z = 20) {
    
    img.date <- toString(img.date)
    obj.nm <- paste0(substring(batch,1,1), ".",img.date)
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
        filenm <- paste0(fpath,img.date,"/",batch,"/",batch,"_",rev.date,"_",i,".tif")
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


#' Load one day's acquisitions from all batches
#'
#' Import one day's images from into  1996x1996x20x3 array.
#' @param img.date Date of images to import, in format yymmdd.
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @export
#' @examples
#' img.150828 <- load.daily(150828)
#' 
#' 
load.daily <- function(img.date, fpath = "/home/clair/Documents/Pixels/Image-data/") {
    m <- array(dim = c(1996, 1996, 20, 3),
               dimnames = list(NULL, NULL, NULL, c("black", "grey", "white")))
    
    m[,,,1] <- load.images(img.date, "black")
    m[,,,2] <- load.images(img.date, "grey")
    m[,,,3] <- load.images(img.date, "white")
    return(m)
}



#' Load xml profiles of one day's acquisitions from one batch
#'
#' Support function: import xml profiles of one day's images into a single data frame
#' @param img.date Date of images to import, in format yymmdd.
#' @param batch Specify batch to import: black, grey or white.
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing all xml profile data for the specified image set.
#' @export
#' @examples
#' b.150828.xml <- load.profiles(150828, "black")
#' 
#' 
load.profiles <- function(img.date, batch, fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    img.date <- toString(img.date)
    
    fpath <- paste(fpath, img.date, "/", batch, "/", sep = "")
    files <- list.files(fpath, pattern = "\\.xml$")
    
    xml.data <- list()
    
    # read through all xml data and import all xml fields into single df per image
    for (i in 1:length(files)) {
        xml.data[[i]] <- as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(paste(fpath, files[i], sep = ""))),
                                                             list(recursive=TRUE))), FUN = unlist),
                                       stringsAsFactors = F)
    }
    
    # merge all list elements into single df, filling in missing columns with NA
    xml.data <- rbind.fill(xml.data)
    
    # where possible, convert character strings to numeric
    for (i in 1:ncol(xml.data)) {
        if (suppressWarnings(all(!is.na(as.numeric(xml.data[,i]))))) {
            
            xml.data[,i] <- as.numeric(xml.data[,i])
            
            # if column contains only whole numbers, convert to integer rather than numeric
            if (all(xml.data[,i] == floor(xml.data[,i]))) {
                xml.data[,i] <- as.integer(xml.data[,i]) 
            }
        }
    }
         
    # sort data frame in chronological order
    xml.data <- xml.data[order(strptime(xml.data$Time, "%H:%M:%S")),]
    return(xml.data)
}



#' Summarise image profile data
#'
#' Searches through folder structure, extracts all .xml profiles, and summarises into a data frame.
#' @details If NA appears in data frame, the acquisitions for that batch do not all share the same value of that attribute, and should be investigated individually.
#' @param fpath Path to top level of image directory. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing summaries of image profile data.
#' @export
#' @examples
#' df <- summarise.profiles()
#' 
#' 
summarise.profiles <- function(fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    # get list of folders in given path location
    folders <- list.dirs(fpath, full.names = F, recursive = F)
    
    # check that valid path has been given
    if (length(folders) == 0) {
        cat("Image folder",fpath,"not found in",getwd())
    } else {
        
        d <- length(folders)
        xml.summ <- list()
        
        # read each folder in turn
        for (i in 1:d) {
            
            batches <- tolower(list.dirs(paste(fpath, folders[i],"/", sep = ""), full.names = F, recursive = F))
            
            # read each batch in turn
            for (j in 1:length(batches)) {
                
                # extract data frame of xml data (may not be the same format for each day!)
                p <- load.profiles(folders[i], batches[j], fpath)
                
                # for numeric columns, check that min == max; otherwise, replace with NA
                for (k in c(1:ncol(p))[sapply(p, class) == "numeric"]) {
                    if (min(p[,k]) != max(p[,k])) {
                        p[,k] <- NA
                    }
                }
                xml.summ[[length(xml.summ) + 1]] <- cbind(date = folders[i],
                                                          batch = batches[j],
                                                          frames = nrow(p),
                                                          p[1,])
            }
        }
        df <- rbind.fill(xml.summ)
        df$date <- as.character(df$date)
        df$batch <- as.character(df$batch)
        df$frames <- as.integer(df$frames)
        return(df)
    }
}



#' Calculate summary statistics for all images
#'
#' Searches through folder structure, extracts all .tif images, and stores the summary statistics for each in a data frame.
#' @param fpath Path to top level of image directory. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing summaries of pixel values: mean and standard deviation, max, min, quartiles and IQR.
#' @export
#' @examples
#' df <- summarise.images()
#' 
#' 
summarise.images <- function(fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
    df <- data.frame(img.date = character(),
                     batch = character(),
                     
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
        
        batches <- tolower(list.dirs(paste(fpath, folders[i],"/", sep = ""), full.names = F, recursive = F))
        
        # read each batch in turn
        for (j in 1:length(batches)) {
            
            r <- nrow(df) + 1
            data <- load.images(folders[i], batches[j], fpath)
            summ <- batch.summary(data)
            
            df[r, c(1:2)] <- c(folders[i], batches[j])
            df[r, c(3:6)] <- summ[[1]]
            df[r, c(7:12)] <- summ[[2]]
        }
    }
    return(df)
}



#' Count available images
#'
#' Searches through folder structure to quickly summarise number of available .tif images from each batch on each day. Assumes files are organised as /Image-data/date/batch/image.tif
#' @param fpath Path to top level of stored images. Default is "/home/clair/Documents/Pixels/Image-data/"
#' @return Data frame containing counts of .tif images in batch subfolders for each date, with image acquisition dates as row names.
#' @export
#' @examples
#' image.count <- count.images()
#' 
#' 
count.images <- function(fpath = "/home/clair/Documents/Pixels/Image-data/") {
    
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
            
            # check folder for each batch in turn, get number of .tif files
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


#' Quickly load files
#'
#' Load all pixelwise summaries and bad pixel map in one step.
#' @param fpath Path to top level of stored images. Default is "./Other-data/"
#' @details Imports pixelwise mean and SD for all files, and bad pixel map
#' @export
#' @examples
#' load.pixel.maps()
#' 
load.pixel.maps <- function(fpath = "./Other-data/") {
    pw.w <<- readRDS(paste0(fpath, "Pixelwise-means-white.rds"))
    pw.g <<- readRDS(paste0(fpath, "Pixelwise-means-grey.rds"))
    pw.b <<- readRDS(paste0(fpath, "Pixelwise-means-black.rds"))
    
    pw.sd.w <<- readRDS(paste0(fpath, "Pixelwise-sds-white.rds"))
    pw.sd.g <<- readRDS(paste0(fpath, "Pixelwise-sds-grey.rds"))
    pw.sd.b <<- readRDS(paste0(fpath, "Pixelwise-sds-black.rds"))
    
    bpm <<- read.csv(paste0(fpath, "BadPixelMap-160314.csv"), as.is = T)
}





#' Import a single day's acquisition into a 2048 x 2048 array
#' 
#' Import all images for one acquisition date into a 2048 x 2048 array, with layers corresponding to different power settings. Cropped pixels are padded with NA.
#' @export
#' @examples
#' im.160430 <- import.acq("./Image-data/160430")
#' 
import.acq <- function(acq.folder, subfolders = c("black", "grey", "white"), panel.dim = c(2048, 2048)) {
    
    acq <- array(dim = c(panel.dim, length(subfolders)),
                 dimnames = list(NULL, NULL, subfolders))
    
    for (sf in subfolders) {
        
        # load xml data into data frame (get image list & offsets)
        im.list <- list.files(paste0(acq.folder, "/", sf), pattern = "\\.xml$", full.names = T)
        df <- rbind.fill(lapply(im.list, 
                                function(im.xml) as.data.frame(lapply(do.call("c", c(xmlToList(xmlParse(im.xml)), list(recursive = TRUE))), 
                                                                      FUN = unlist), stringsAsFactors = F)))
        # gives ordered list of images - may be needed if investigating as time series
        df <- data.frame(acq.time = df$Time,
                         filenm = df$.attrs.FileName,
                         offset.x = df$CameraProperties.imageOffsetX,
                         offset.y = df$CameraProperties.imageOffsetY,
                         size.x = df$CameraProperties.imageSizeX,
                         size.y = df$CameraProperties.imageSizeY,
                         kV = df$kV,
                         uA = df$uA,
                         exp.time = df$ExposureTime,
                         stringsAsFactors = F)[order(acq.time = df$Time),]
        
        # quick error check - are all images on same area of detector (should always be T)
        if (!identical(sapply(df[3:6], min), sapply(df[3:6], max))) {
            cat("Error - TIF images do not all have same offset & dimension")
            return(df)
        }
        
        # if all images cover same area of detector, continue
        offset.x <- as.integer(df$offset.x[1])
        offset.y <- as.integer(df$offset.y[1])
        size.x <- as.integer(df$size.x[1])
        size.y <- as.integer(df$size.y[1])
        
        ims <- apply(cbind(acq.folder, "/", sf, "/", df$filenm), 1, paste, collapse = "")
        
        # get pixelwise mean of image
        tmp <- apply(abind(lapply(ims, readTIFF, as.is = T), along = 3), 1:2, mean)
        tmp <- t(tmp[nrow(tmp):1, , drop = FALSE])
        
        # position correctly in array
        acq[(offset.x + 1) : (offset.x+size.x),
            (panel.dim[2]-size.y-offset.y + 1) : (panel.dim[2]-offset.y),
            sf] <- tmp
    }
    return(acq)
}