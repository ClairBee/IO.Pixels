
#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at mean and +- 0.5, 1, 2, 3 sd.
#' @param mean Mean value to use as midpoint of levels
#' @param sd Standard deviation value to determine distance of cutpoints from mean
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'     data,
#'     levels = sd.levels(mean(data), sd(data)))
#' 
#' 
sd.levels <- function(mean, sd) {
    c(0, mean + (c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3) * sd), 1)
}



#' Colour scheme for contour plot
#'
#' Support function: returns a vector of colours to be used in contour plotting.
#' @return Vector of colours. When used in conjunction with  \link{\code{sd.levels()}}, 
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'     data,
#'     levels = sd.levels(mean(data), sd(data)),
#'     col = sd.colours())
#' 
#' 
sd.colours <- function() {
    c("midnightblue",     # 0 to -3sd
      "blue",             # -3 to -2sd
      "dodgerblue3",      # -2 to -1sd
      
      "greenyellow",       # -1 to -0.5sd
      "yellow",           # 0.5 to 0sd
      "gold",             # 0 to 0.5sd
      "orange",           # 0.5 to 1sd
      
      "red3",              # 1 to 2sd
      "violetred",        # 2 to 3sd
      "purple")           # 3sd to 1
}



#' Contour plot of pixel values
#'
#' Create a filled contour plot with pixels shaded according to their distance from the mean value. Takes approx. 10x as long as \link{pixel.image}.
#' @details Pixels more than 2sd above the mean are red, pink or purple. Pixels within 2sd of the mean are orange, yellow or green. Pixels more than 2sd below the mean are blue.
#' @param data 2-dimensional matrix containing values to be plotted
#' @param title String containing title to be printed with plot
#' @param midpoint Either 'mean' or 'median'. Default is 'mean'
#' @export
#' @examples
#' pixel.contour(pw.mean)
#' 
#' 
pixel.contour <- function(data, title = "", midpoint = "mean") {
    
    # if data has more than 2 dimensions, truncate to take only first 'layer'
    d <- length(dim(data))
    if (d > 2) {
        data <- data[,, rep(1, d-2)]
    }
    
    if (midpoint == "median") {
        m <- median(data)
    } else {
        m <- mean(data)
    }
    s <- sd(data)
    
    filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
                   data,
                   levels = sd.levels(m, s),
                   col = sd.colours(),
                   asp = T,
                   main = title
    )
}


#' Display pixel values as image
#'
#' Create an image with pixels shaded according to their distance from the mean value. Displays at lower resolution than contour plot (plots are 20% of the size), but much quicker (10% of time taken).
#' @details Pixels more than 2sd above the mean are red, pink or purple. Pixels within 2sd of the mean are orange, yellow or green. Pixels more than 2sd below the mean are blue.
#' @param data 2-dimensional matrix containing values to be plotted
#' @param x.range Vector range showing x-range of cells to be included in plot
#' @param y.range Vector range showing y-range of cells to be included in plot
#' @param title String containing title to be printed with plot
#' @param midpoint Either 'mean' or 'median'. Default is 'mean'
#' @export
#' @examples
#' pixel.image(pw.mean)
#' 
#' 
pixel.image <- function(data, title = "", x.range = c(1:1996), y.range = c(1:1996), midpoint = "mean") {
    
    if (midpoint == "median") {
        m <- median(data)
    } else {
        m <- mean(data)
    }
    s <- sd(data)
    
    image(x.range, y.range, 
          data[x.range, y.range], 
          col = sd.colours(),
          breaks = sd.levels(m, s),
          main = title,
          asp = T)
}



#' Add panel edges to pixel image
#'
#' Add lines showing the edges of individual detector panels to a pixel value image (eg. created using \link{pixel.image})
#' @details Each detector panel outputs 128 x 1024 pixels, with some cropping at each edge of the combined image.
#' @param width Width of each panel, in pixels. Default is 128.
#' @param height Height of each panel, in pixels. Default is 1024.
#' @param left.crop Number of pixels cropped from left-hand edge of image. Default is 2.
#' @param left.crop Number of pixels cropped from lower edge of image. Default is 32.
#' @export
#' @examples
#' pixel.image(pw.mean)
#' show.panels()
#' 
#' 
show.panels <- function(left.crop = 2, lower.crop = 32, width = 128, height = 1024) {
    abline(h = height - lower.crop + 0.5)
    abline(v = (c(0:15)*128 - left.crop + 0.5))
}