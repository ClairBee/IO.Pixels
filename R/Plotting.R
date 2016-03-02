
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
sd.levels <- function(data, midpoint = "mean") {
    
    if (midpoint == "median") {
        m <- median(data)
    } else {
        m <- mean(data)
    }
    
    c(min(data), m + (c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3) * sd(data)), max(data))
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
        
    filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
                   data,
                   levels = sd.levels(data, midpoint),
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
#' @param break.levels Vector of values used as breakpoints when binning values
#' @export
#' @examples
#' pixel.image(pw.mean)
#' 
#' 
pixel.image <- function(data, title = "", x.range, y.range, midpoint = "mean", break.levels) {
        
    if (missing("break.levels")) {break.levels <- sd.levels(data, midpoint)}
    
    if (missing("x.range")) {x.range <- c(1:nrow(data))}
        
    if (missing("y.range")) {y.range <- c(1:ncol(data))}
    
    image(x.range, y.range, 
          data[x.range, y.range], 
          col = sd.colours(),
          breaks = break.levels,
          main = title,
          asp = T)
}



#' Get coordinates of panel edges
#'
#' Identify coordinates of individual detector panels, given their dimensions and those of the whole sensor.
#' @details Each detector panel outputs 128 x 1024 pixels, with some cropping at each edge of the combined image.
#' @param left.crop Number of pixels cropped from left-hand edge of image. Default is 2.
#' @param upper.crop Number of pixels cropped from lower edge of image. Default is 20.
#' @param width Width of each panel, in pixels. Default is 128.
#' @param height Height of each panel, in pixels. Default is 1024.
#' @param x.dim Width of image returned, in pixels. Default is 1996.
#' @param y.dim Height of image returned, in pixels. Default is 1996.
#' @export
#' @examples
#' panels <- panel.edges()
#' 
#' 
panel.edges <- function(left.crop = 2, upper.crop = 20, width = 128, height = 1024, x.dim = 1996, y.dim = 1996) {
    
    list(y = c(1, y.dim - height + upper.crop + 1, y.dim + 1),
         x = c(1, c(1:15)* width - left.crop + 1, x.dim + 1))
}


#' Add panel edges to plot
#'
#' Mark edges of individual panels in array
#' @details Each detector panel outputs 128 x 1024 pixels, with some cropping at each edge of the combined image.
#' @param p List containing y, a named vector of y-coordinates of panel corners, and x, a named vector of x-coordinates of panel corners.
#' @export
#' @examples
#' pixel.image(b.150828)
#' draw.panels()            # marks specified panel boundaries
#' 
#' 
draw.panels <- function(p = panel.edges()) {
    abline(h = p$y - 0.5)
    abline(v = p$x - 0.5)
}
