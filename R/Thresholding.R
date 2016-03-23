

#' Identify 'bad' pixels
#'
#' Use thresholds set out in technical manual to identify 'bad' pixels
#' @param data Array containing a sequence of images from a single batch
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


#' Identify extreme values along transect
#' 
#' Find extreme pixel values in single vector
#' @details Uses Loess smoothing to find 'spine' of vector and identifies any points that fall more than a certain constant multiple of the SD from the central spine.
#' @param transect Vector of pixel values to be assessed for 'bad' pixels
#' @param loess.span Filter span to be used in smoothing the data. Default is 1/15.
#' @param keep Central proportion of points to keep. Default is 0.95, which results in a cutoff for 'bad' pixels of 1.96sd.
#' @param display Boolean: plot the resulting transect? Default is T.
#' @param legend Boolean: include legend in plot? Default is F.
#' @param legend.pos String: where to display legend? default is "topright".
#' @param ... Additional parameters to be passed to plotting function
#' @return A list containing named vectors of Loess-smoothed values, residuals, and indices of pixels flagged.
#' @export
#' @examples 
#' zz <- bp.by.res.sd(pw.m[992, 1:992], display = T)
bp.by.res.sd <- function(transect, loess.span = 1/15, keep = 0.95, display = T, legend = F, legend.pos = "topright",  ...) {
    
    k <- qnorm(1-(1-keep)/2, 0, 1)
    
    smoo <- lowess(transect, f = loess.span)$y
    res <- transect - smoo
    
    abs.mean <- mean(abs(res))
    abs.sd <- sd(abs(res))
    
    bp <- which(abs(res) > (abs.mean + k * abs.sd))
    
    if (display) {
        plot(transect, type = "l", col = adjustcolor("grey", alpha = 0.5), ...)
        points(transect, pch = 20, cex = 0.7)
        
        points(smoo, type = "l", col = adjustcolor("green3", alpha = 0.5), lwd = 3)
        points(smoo + abs.mean, type = "l", col = "blue", lwd = 2)
        points(smoo - abs.mean, type = "l", col = "blue", lwd = 2)
        
        points(smoo - abs.mean - (k * abs.sd), type = "l", col = "purple", lwd = 2)
        points(smoo + abs.mean + (k * abs.sd), type = "l", col = "purple", lwd = 2)
        
        points(cbind(bp, transect[bp], ncol = 2), col = "red")
        
        if (legend) {
            legend(legend.pos, bty = "n", lwd = 2,
                   col = c("green3", "blue", "purple", "red"), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 1), 
                   legend = c("Smoothed value", "Smoothed + mean amplitude", 
                              paste0("Smoothed + mean amplitude + ", round(k, 2), " SD"),
                              "Bad pixels"))
        }
    }
    
    list(smoothed = smoo, residual = res, flagged = bp)
} 



#' Identify extreme values along transect
#' 
#' Find extreme pixel values in single vector
#' @details Uses Loess smoothing to find 'spine' of vector and identifies any points that fall more than a certain constant multiple of the NMAD from the central spine.
#' @param transect Vector of pixel values to be assessed for 'bad' pixels
#' @param loess.span Filter span to be used in smoothing the data. Default is 1/15.
#' @param keep Central proportion of points to keep. Default is 0.95, which results in a cutoff for 'bad' pixels of 1.96sd.
#' @param display Boolean: plot the resulting transect? Default is T.
#' @param legend Boolean: include legend in plot? Default is F.
#' @param legend.pos String: where to display legend? default is "topright".
#' @param ... Additional parameters to be passed to plotting function
#' @return A list containing named vectors of Loess-smoothed values, residuals, and indices of pixels flagged.
#' @export
#' @examples 
#' zz <- bp.by.res.sd(pw.m[992, 1:992], display = T)
bp.by.res.nmad <- function(transect, loess.span = 1/15, keep = 0.95, display = T, legend = F, legend.pos = "topright",  ...) {
    
    k <- qnorm(1-(1-keep)/2, 0, 1)
    
    smoo <- lowess(transect, f = loess.span)$y
    res <- transect - smoo
    
    abs.mean <- mean(abs(res))
    abs.mad <- mad(abs(res))
    
    bp <- which(abs(res) > (abs.mean + k * abs.mad))
    
    if (display) {
        plot(transect, type = "l", col = adjustcolor("grey", alpha = 0.5), ...)
        points(transect, pch = 20, cex = 0.7)
        
        points(smoo, type = "l", col = adjustcolor("green3", alpha = 0.5), lwd = 3)
        
        points(smoo + abs.mean, type = "l", col = "blue", lwd = 2)
        points(smoo - abs.mean, type = "l", col = "blue", lwd = 2)
        
        points(smoo - abs.mean - (k * abs.mad), type = "l", col = "purple", lwd = 2)
        points(smoo + abs.mean + (k * abs.mad), type = "l", col = "purple", lwd = 2)
        
        points(cbind(bp, transect[bp], ncol = 2), col = "red")
        
        if (legend) {
            legend(legend.pos, bty = "n", lwd = 2,
                   col = c("green3", "blue", "purple", "red"), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 1), 
                   legend = c("Smoothed value", "Smoothed + mean amplitude", 
                              paste0("Smoothed + mean amplitude + ", round(k, 2), " NMAD"),
                              "Bad pixels"))
        }
    }
    
    list(smoothed = smoo, residual = res, flagged = bp)
} 


#' Identify extreme values in image
#' 
#' Find columnwise extreme values in image array
#' @details For each column of each panel, use Loess smoothing to find central spline and calculate residuals. Any points with a residual above a certain threshold are flagged as bad pixels.
#' @param image Single-layer 1996x1996 image array
#' @param method String: either "SD" or "NMAD". Threshold will be set according to multiples of SD or NMAD of columnwise residuals.
#' @param n.sigma Numeric scalar: how many multiples of the SD/NMAD will be used as threshold for 'extreme' values.
#' @param span Numeric value: smoothing span to be used in Loess spline
#' @param display Boolean: plot the results?
#' @param panels Boolean: if plotting the results, show subpanel edges?
#' @export
#' 
get.bad.pixels <- function(image, method = "sd", n.sigma = 1.96, span = 1/15, display = T, panels = T) {
        
    # subfunction to return indices of extreme values 
    bad.pixels <- function(transect, n.sigma, span, method) {
        smoo <- lowess(transect, f = span)$y
        res <- transect - smoo
        
        abs.mean <- mean(abs(res))
        if (method %in% c("mad", "nmad", "MAD", "NMAD")) {
            bound <- mad(abs(res))
        } else {
            bound <- sd(abs(res))
        }
        
        which(abs(res) > (abs.mean + n.sigma * bound))
    }
    
    # get coordinates of selected pixels
    lower.out <- apply(image[, 1:992], 1, bad.pixels, n.sigma, span, method)
    
    # label each list of pixel indices with column ref
    names(lower.out) <- c(1:1996)
    
    # remove any empty elements from list
    lower.out <- Filter(length, lower.out)
    
    # convert per-column indices to coordinates
    lower.coords <- do.call("rbind", Map(function(lower.out, i) cbind(col = as.integer(i), row = lower.out), lower.out, names(lower.out)))
    
    # repeat for upper panel
    upper.out <- apply(image[, 993:1996], 1, bad.pixels, n.sigma, span, method)
    names(upper.out) <- c(1:1996)
    upper.out <- Filter(length, upper.out)
    upper.coords <- do.call("rbind", Map(function(upper.out, i) cbind(col = as.integer(i), row = upper.out + 992), upper.out, names(upper.out)))
    
    coords <- rbind(lower.coords, upper.coords)
    
    if (display) {
        plot(coords, pch = 20, col = adjustcolor("blue", alpha = 0.3), 
             main = paste0("Bad pixels: ", n.sigma, method, " from smoothed column spline"))
        if (panels) draw.panels()
    }
    
    list(coords = coords, method = method, n.sigma = n.sigma, span = span)
}