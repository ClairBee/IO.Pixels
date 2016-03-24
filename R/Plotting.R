
#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at mean and +- 0.5, 1, 2, 3 sd.
#' @param data Matrix of data to plot
#' @param midpoint String: use mean or median as midpoint? Default is "mean".
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'                    data,
#'                    levels = sd.levels(data))
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


#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at outlier (lq - 1.5*IQR), lq, 3 even bins between lq and median, median, 3 even bins between median and uq, uq, and outlier (1.5*iqr + uq).
#' @param data Matrix of data to plot
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'                    data,
#'                    levels = iqr.levels(data))
#' 
#' 
iqr.levels <- function(data) {
    lq <- quantile(data, 0.25)
    med <- median(data)
    uq <- quantile(data, 0.75)
    iqr <- IQR(data)
    
    c(0,
      max(lq - 1.5 * iqr,         # low outliers
          lq/2),
      lq,
      lq + (med - lq) / 3,
      med - (med - lq) / 3,
      med,
      med + (uq - med) / 3,
      uq - (uq - med) / 3,
      uq,
      min(uq + 1.5 * iqr,    # high outliers
          uq + ((65535-uq) / 2)),
      65535)
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
pixel.image <- function(data, title = "", x.range = c(1:nrow(data)), y.range = c(1:ncol(data)), 
                        midpoint = "mean", break.levels = sd.levels(data, "mean"), 
                        panels = F, x.lab = "", y.lab = "", ...) {
    
    image(x.range, y.range, 
          data[x.range, y.range], 
          col = sd.colours(),
          breaks = break.levels,
          main = title, xlab = x.lab, ylab = y.lab,
          asp = T, 
          ...)
    
    if (panels) draw.panels()
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


#' Vector overlay plot
#'
#' Base plotting function, specifying type "o" (points plotted over lines), pch = 20, cex = 0.7. Other plotting parameters can be passed as arguments as usual.
#' @param data Vector of values to plot (eg. column of pixel values)
#' @param add Boolean: add values to current plot, or start a new plot? Default is F.
#' @export
#' @examples
#' o.plot(b.150828[4,,1])
#' o.plot(b.150828[5,,1], col = adjustcolor("red", alpha = 0.2), add = T)
#' 
#' 
o.plot <- function(data, add = F, ...) {
    if (add) {
        points(data, type = "o", pch = 20, cex = 0.7, ...)
    } else {
        plot(data, type = "o", pch = 20, cex = 0.7, ...)
    }
}



#' Histogram plot with SD levels shown
#' 
#' Producs a histogram with coloured bars at the bottom showing scale used. Defaults can be used to show legend for plots produced using \link{\code{pixel.image}}.
#' @param data Array of numbers to be plotting in histogram
#' @param scale Vector of scale cutpoints. Default is \link{\code{sd.levels()}}.
#' @param scale.colours Vector of scale colours. Default it \link{\code{sd.lcolours()}}
#' @param xlim Vector of lower and upper x-limits. If not provided, will use central 95% of normal distribution with mean & sd of observed data.
#' @export
#' @examples
#' pixel.image(b.150828)
#' scale.hist(b.150828)
s.hist <- function(data, scale = sd.levels(data), scale.colours = sd.colours(), xlim, ...) {
    if (missing(xlim)) {
        xlim <- c(floor(qnorm(0.025, mean(data), sd(data))/10)*10,
                  ceiling(qnorm(0.975, mean(data), sd(data))/10)*10)
    }
    
    hist(data, breaks = "fd", xlim = xlim, ...)
    
    # add colours to indicate scale of pixel map
    cl <- cut(xlim[1]:xlim[2], scale)
    points(xlim[1]:xlim[2], rep(-600, length(xlim[1]:xlim[2])), pch = 15, col = scale.colours[cl])
}