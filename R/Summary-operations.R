

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
    # alter this to allow pixelwise mean over 4d array?
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
    # alter this to allow pixelwise sd over 4d array?
    apply(data, c(1,2), sd)
}



#' Get summary statistics for an image batch
#'
#' Return a vector of summary statistics for a batch of acquisitions: mean, standard deviation, median, upper and lower quartiles, interquartile range, and max and min pixelwise mean value.
#' @param data Array with 3 or more dimensions
#' @return Vector containing named summary statistics
#' @details Pixelwise max and min value are used to check if image is scaled from absolute 0:1 or from an arbitrary minimum point
#' @export
#' @examples
#' b.150828.summ <- batch.summary(b.150828)
#' 
#' 
batch.summary <- function(data) {
    
    # output vector of summary statistics
    list(c(mean = mean(data),
           sd = sd(data)),
         c(min = min(data),
           lq = quantile(data,0.25),
           median = median(data),
           uq = quantile(data, 0.75),
           max = max(data),
           iqr = IQR(data)))
}


#' Store summary statistics for an image batch in a data frame
#'
#' Populate a data frame with the summary statistics provided by \link{batch.summary}
#' @param data Array with 3 or more dimensions
#' @return Updates global data frame containing summary statistics for multiple images.
#' @export
#' @examples
#' save.summary(b.150828)
#' 
#' 
save.summary <- function(data) {
    
    obj.nm <- toString(as.list(sys.call())[2])
    img.date <- substring(obj.nm,3,8)

    
    if (!exists("img.summs")) {
        img.summs <<- data.frame(obj.nm = character(),
                                 img.date = character(),
                                 channel = character(),
                                 g.mean = numeric(),
                                 g.sd = numeric(),
                                 g.lq = numeric(),
                                 g.median = numeric(),
                                 g.uq = numeric(),
                                 g.iqr = numeric(),
                                 pw.min = numeric(),
                                 pw.max = numeric(),
                                 stringsAsFactors = F)
        r <- 1
    } else {
        if (obj.nm %in% img.summs$obj.nm) {
            cat("Row already exists. Delete it manually if you want to overwrite.")
            return()
        } else {
            r <- nrow(img.summs) + 1
        }
    }
    
    img.date <- substring(obj.nm,3,8)  
    
    channels <- matrix(c("b", "g", "w", "black", "grey", "white"), ncol = 2, nrow = 3)
    channel <- channels[which(channels[,1] == substring(obj.nm,1,1)),2]
    
    s <- batch.summary(data)
    
    img.summs[r, c(1:3)] <<- c(obj.nm, img.date, channel)
    img.summs[r,c(4:11)] <<- s
}