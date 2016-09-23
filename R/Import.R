
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