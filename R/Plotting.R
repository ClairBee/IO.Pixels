
#' Plot given coordinates
#'
#' Plot pixel map
#' @param px Coordinates of pixels to be plotted
#' @export
#' @examples
#' pixel.plot(which(pw.m[,,"white"] < 15000, arr.ind = T))
#' 
pixel.plot <- function(px, xlim = c(0,2048), ylim = c(0,2048), pch = 15, panels = F, cex = 0.4,
                       main = "", xlab = "", ylab = "", ...) {
    
    plot(px[,1:2], asp = T, 
         xlim = xlim, ylim = ylim, pch = pch, main = main, xlab = xlab, ylab = ylab, cex = cex, ...)
}


#' Colours for plotting bad pixels
#' 
#' Support function: colour scheme for given pixel types
#' @return List of colours
#' @export
#' 
px.cols <- function() {
    c("magenta3", "black", "purple", "red", "orange", "gold", 
       "blue", "skyblue", "green1", "green3", "cyan1", "cyan3")
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
#' 
pixel.image <- function(data, main = "", x.range = c(1:nrow(data)), y.range = c(1:ncol(data)), 
                        break.levels = th.levels(data), panels = F, x.lab = "", y.lab = "", ...) {
    
    image(x.range, y.range, 
          data[x.range, y.range], 
          col = sd.colours(),
          breaks = break.levels,
          main = main, xlab = x.lab, ylab = y.lab,
          asp = T, 
          ...)
    
    if (panels) draw.panels()
}


#' Set contour plot levels
#'
#' Support function: returns a vector of levels to be used in contour plotting, cutting data at mean and +- 1,2,3,4.5,6 mad.
#' @param data Matrix of data to plot
#' @return Vector of calculated levels
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'                    data,
#'                    levels = sd.levels(data))
#' 
th.levels <- function(data) {
    
    data <- data[!is.na(data)]
    data <- data[!is.infinite(data)]
       
    m <- modal.density(data)

    sort(c(sapply(c(1, 2, 3, 4.5, 6), function(nn) asymmetric.mad(data, nn)),
           m,
           range(data)))
}


#' Colour scheme for contour plot
#'
#' Support function: returns a vector of colours to be used in contour plotting.
#' @return Vector of colours. When used in conjunction with  \link{\code{th.levels()}}, shades thresholded values.
#' @export
#' @examples
#'     filled.contour(x = c(1:dim(data)[1]), y = c(1:dim(data)[2]),
#'     data,
#'     levels = sd.levels(mean(data), sd(data)),
#'     col = sd.colours())
#' 
#' 
sd.colours <- function() {
    c("black",             # min to -6sd        # min to -6mad
      "darkblue",     # -6 to -3sd              # -6 to -5mad
      "blue",             # -3 to -2sd          # -5 to -3mad
      
      "cyan3",      # -2 to -1sd                # -3 to -2mad
      
      "green1",       # -1 to -0.5sd            # -2 to -1mad
      "greenyellow",    # 0.5 to 0sd            # -1 to 0mad
      "yellow",             # 0 to 0.5sd        # 0 to 1mad
      "gold",           # 0.5 to 1sd            # 1 to 2mad
      
      "orange",              # 1 to 2sd         # 2 to 3mad
      
      "red3",        # 2 to 3sd                 # 3 to 5mad
      "purple",           # 3 to 6sd            # 5 to 6mad
      "deeppink")         # 6sd to max          # 6mad to max  
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
draw.panels <- function(p = panel.edges(), ...) {
    
    # horizontal
    for (yl in p$y[2:(length(p$y)-1)]-0.5) {
        lines(range(p$x)-0.5, c(yl, yl), ...)
    }
    
    # vertical
    for (xl in p$x[2:(length(p$x)-1)]-0.5) {
        lines(c(xl, xl), range(p$y)-0.5, ...)
    }
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
panel.edges <- function(left.crop = 0, upper.crop = 0, width = 128, height = 1024, x.dim = 2048, y.dim = 2048) {
    
    #    seq(from = 0, to = y.dim, by = height) + y.dim - height + upper.crop + 1
    
    list(y = c(1, y.dim - height + upper.crop + 1, y.dim + 1),
         x = c(1, c(1:15)* width - left.crop + 1, x.dim + 1))
}


#' Draw outlines of convex hulls
#' 
#' Given a set of coordinates, identify each cluster of adjacent pixels and draw the convex hull.
#' @param px Matrix or data frame of x and y coordinates to be clumped and outlined
#' 
#' @export
#' 
draw.outlines <- function(px, im.dim = c(2048, 2048), ...) {
    
    # clump adjacent pixels
    cc <- clump(px2r(px, im.dim = im.dim), dir = 4)
    
    # coordinates of each clump, with clump id
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    invisible(lapply(unique(xy$id),
                     function(i) {
                         ch <- chull(xy[xy$id == i,1:2])
                         lines(xy[xy$id == i,1:2][c(ch, ch[1]),], ...)
                     }))
}
