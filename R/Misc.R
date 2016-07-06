
#' Create sample of 'healthy' pixels
#'
#' Given a matrix of coordinates of bad pixels, create a same-size random sample of healthy pixels that do not appear on the bad pixel list.
#' @param bad.pixels Data frame or matrix whose first two columns are X and Y coordinates of identified bad pixels.
#' @param n Size of sample to create. Default is to create a sample of same size as bad pixel list.
#' @return Matrix containing X and Y coordinates of sampled pixels.
#' @export
#' @examples
#' bp.160314 <- reset.bp(160314)
#' samp.160314 <- sample.healthy(bp.160314)
#' 
#' 
sample.healthy <- function(bad.pixels, n = nrow(bad.pixels)) {
    
    # convert coords to single integer
    # subtract 1997 to because coords start at 1, not 9
    bp.int <- (bad.pixels[, 1] * 1996) + bad.pixels[, 2] - 1997

    # sample as integer (much more efficient than permuting)
    s <- sample(c(1:(1996^2))[-c(bp.int)], n)
    
    # decompose into integer division & modulo components to get x & y
    # add 1 to get coordinates to correct range
    cbind(x = s %/% 1996 + 1, y = s %% 1996 + 1)
}


#' Convert bad pixel map to image
#' 
#' Given a bad pixel map with coordinates and types, create a 2d array showing each type of bad pixel identified
#' @param bpx Bad pixel map: data frame with first two columns containing x and y coordinates of bad pixels, and a column 'type' specifying what kind of bad pixel was identified.
#' @param im.dim Vector of dimensions of image array to create. Default is c(1996, 1996).
#' @return 2d image array of the specified dimensions
#' @export
#' 
bpx2im <- function(bpx, im.dim = c(1996, 1996)) {
    im <- array(0, dim = im.dim)
    im[as.matrix(bpx[,1:2])] <- bpx$type
    return(im)
}





#' Combine two bad pixel maps
#'
#' Merge two data frames of bad pixels into one single list of coordinates, with source map and (where possible) error type marked.
#'
#' NB. Currently only set up to compare 'official' (old) pixel map with new one. Should update this.
#' @param new.map Data frame (as created using \link{\code{reset.bp}}) containing coordinates and types of pixels.
#' @param old.map Data frame (eg. official bad pixel map) containing coordinates and types of pixels.
#' @return Data frame containing row and column coordinates of bad pixels, factor giving type of bad pixel, and factor noting which map the pixel was found in (old, new, or both)
#' @export
#' @examples
#' matched <- combine.maps(bp, bpm)
combine.maps <- function(new.map, old.map) {
    
    match <- merge(data.frame(new.map[!duplicated(new.map[,1:2]), c("row", "col", "type")],
                      map1 = "new"),
                 data.frame(old.map, map2 = "old"),
                 by.x = c(1:2), by.y = c(1:2), all = T)
    
    match$map <- apply(match[,c("map1", "map2")], 1, function(x) paste(na.omit(x), collapse = ""))
    match$map[match$map == "newold"] <- "both" 
    match$map <- factor(match$map, levels = c("new", "both", "old"))
    
    match[order(match$map), c("row", "col", "type", "map")]
}



#' Plot a combined pixel map
#'
#' Compare bad pixels identified by two different pixel maps by plotting their values.
#' @param map Data frame (as created using \link{\code{combine.maps}}) containing coordinates and sources of pixels.
#' @param img.date Integer or string giving acquisition date for which bad pixels are to be plotted.
#' @param colour Which image batch should be plotted? Default is "all", which will plot all three sequentially.
#' @export
#' @examples
#' matched <- combine.maps(bp, bpm)
#' plot.map.comparison(matched)
plot.map.comparison <- function(map, img.date, colour = "all") {
    
    # get random sample of 'healthy' pixels for comparison
    pp <- sample.healthy(map)
    
    img.date <- toString(img.date)
    
    if (colour == "all") {
        
        for (colour in c("black", "grey", "white")) {
            plot(pw.m[,, colour, img.date][pp], pch = 20, col = adjustcolor("gold", alpha = 0.2),
                 ylim = c(0, 65535), ylab = "Mean pixel value over 20 images", xlab = "Pixel index",
                 main = paste0("Bad pixels - ", colour, " images, ", img.date))
            
            points(pw.m[,, colour, img.date][as.matrix(map[,1:2])], pch = 20, 
                   col = adjustcolor(c("chartreuse3", "slateblue1", "red")[map$map], alpha = 0.5))
            
            legend("topleft", legend = c(levels(map$map), "'healthy' pixels"), pch = 20,
                   col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
        }
    } else {
        
        plot(pw.m[,, colour, img.date][pp], pch = 20, col = adjustcolor("gold", alpha = 0.2),
             ylim = c(0, 65535), ylab = "Mean pixel value over 20 images", xlab = "Pixel index",
             main = paste0("Bad pixels - ", colour, " images, ", img.date))
        
        points(pw.m[,, colour, img.date][as.matrix(map[,1:2])], pch = 20, 
               col = adjustcolor(c("chartreuse3", "slateblue1", "red")[bp.matched$map], alpha = 0.5))
        
        legend("topleft", legend = c(levels(bp.matched$map), "'healthy' pixels"), pch = 20,
               col =  adjustcolor(c("chartreuse3", "slateblue1", "red", "gold"), alpha = 0.5), bty = "n")
    }
}


#' Johnson Q-Q plot
#'
#' Produce a Q-Q plot of data quantiles vs Johnson distribution quantiles, with quantiles of interest marked.
#' @param data Image array of 1996 x 1996 pixels
#' @param quantiles Vector of quantiles to plot; default is c(1:999)/1000
#' @param grid.quantiles Vector of quantiles of particular interest, to mark with gridlines
#' @param title String containing plot title. Default is "Johnson Q-Q plot".
#' @param ... Additional plotting arguments
#' @export
#' @examples
#' Johnson.QQ(pw.sd[,,"black", "150828"])
Johnson.QQ <- function(data, quantiles = c(1:999)/1000, grid.quantiles = c(0.01, 0.99), title = "Johnson Q-Q plot", ...) {
    
    jf <- JohnsonFit(data, moment = "quant")
    
    plot(qJohnson(quantiles, jf), quantile(data, quantiles), pch = 20, asp = T,
         ylab = "Observed quantile", xlab = "Johnson quantile", main = title, ...)
    abline(0,1, col = "red", lty = 2)
    
    abline(h = quantile(data, grid.quantiles), col = "skyblue", lty = 2)
    abline(v = qJohnson(grid.quantiles, jf), col = "skyblue", lty = 2)
    
    legend("topleft", cex = 0.7, title = "Gridlines: ", legend = paste0("Q", grid.quantiles * 100), lty = 2, col = "skyblue")
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



#' Get median differences over image
#' 
#' Run a median filter over an image to obtain a smoothed version, then subtract this from the original image to highlight any points that are unusually high in relation to their neighbours
#' @param im Single-layer image array
#' @return Image array containing pixelwise differences
#' @export
#' 
med.diffs <- function(im) {
    im - r2m(focal(m2r(im), matrix(rep(1, 9), ncol = 3), fun = median))
}


#' Summarise bad pixels for each acquisition
#' 
#' @export
#' 
summarise.bpx <- function(bp, by = "type") {
    
    df <- rbind.fill(lapply(lapply(bp, "[", by), function(m) data.frame(t(as.matrix(table(m))))))
    rownames(df) <- names(bp)
    df
}



#' Apply shading correction
#' 
#' Calculate a shading-corrected image from a set of blanks
#' @param im Set of pixelwise means for black, grey and white acquisitions
#' 
#' @export
#' 
shading.corrected <- function(im) {
    60000 * (im[,,"grey"] - im[,,"black"]) / (im[,,"white"] - im[,,"black"])
}