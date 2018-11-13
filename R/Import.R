
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
    

    
    if (length(acq.list) == 1) {
        imp <- readRDS(paste0(fpath, "/", otype, acq.list, ".rds"))
    } else {
        imp <- abind(sapply(acq.list, function(nm) readRDS(paste0(fpath, "/", otype, nm, ".rds")), 
                            simplify = F),
                     rev.along = 0)
    }
    return(imp)
}



#' Import system-generated bad pixel map
#'
#' Extract coordinates of bad pixels from system-generated bad pixel map
#' @param dt Date of image for which pixel map is to be imported.
#' @param fpath Path to image folder
#' @return Data.frame containing x and y coordinates 
system.map <- function(dt, fpath = "./Image-data/", img = paste0("./02_Objects/images/pwm-", dt, ".rds")) {
    map <- setNames(as.data.frame(do.call("rbind", 
                                   lapply(lapply(xmlToList(xmlParse(paste0(fpath, dt, "/BadPixelMap.bpm.xml")))$BadPixels, 
                                                 "[", 1:2), 
                                          as.numeric))), 
             nm = c("x", "y"))
    
    # use processed image to identify coordinate offsets
    os <- apply(which(!is.na(readRDS(img)[,,1]), arr.ind = T), 2, min)
    
    map$y <- dim(im)[2] - map$y - (os[2] - 1)
    map$x <- map$x + os[1]
    map
}
