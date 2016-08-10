
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
    
    data <- data[!is.na(data)]
    data <- data[!is.infinite(data)]
    
    if (midpoint == "median") {
        m <- median(data)
    } else {
        m <- mean(data)
    }
    
    sort(c(min(data), m + (c(-6, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 6) * sd(data)), max(data)))
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
    c("black",             # min to -6sd
      "darkblue",     # -6 to -3sd
      "blue",             # -3 to -2sd
      
      "cyan3",      # -2 to -1sd
      
      "green1",       # -1 to -0.5sd
      "greenyellow",    # 0.5 to 0sd
      "yellow",             # 0 to 0.5sd
      "gold",           # 0.5 to 1sd
      
      "orange",              # 1 to 2sd
      
      "red3",        # 2 to 3sd
      "purple",           # 3 to 6sd
      "deeppink")         # 6sd to max    
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


#' Draw outlines of convex hulls
#' 
#' Given a set of coordinates, identify each cluster of adjacent pixels and draw the convex hull.
#' @param px Matrix or data frame of x and y coordinates to be clumped and outlined
#' 
#' @export
#' 
draw.outlines <- function(px, im.dim = c(2048, 2048), ...) {
    
    # clump adjacent pixels
    cc <- clump(m2r(bpx2im(data.frame(px, type = 1), im.dim = im.dim)), dir = 4)
    
    # coordinates of each clump, with clump id
    xy <- data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                     id = getValues(cc)[!is.na(getValues(cc))])
    
    invisible(lapply(unique(xy$id),
           function(i) {
               ch <- chull(xy[xy$id == i,1:2])
               lines(xy[xy$id == i,1:2][c(ch, ch[1]),], ...)
           }))
}



#' Plot defective pixels
#'
#' Plot bad pixel map on 2048 x 2048 grid
#' @param px Coordinates of pixels to be plotted
#' @export
#' @examples
#' pixel.plot(which(pw.mean > 25000, arr.ind = T))
#' 
#' 
pixel.plot <- function(px, xlim = c(0,2048), ylim = c(0,2048), pch = 15, panels = F, cex = 0.4,
                       main = "", xlab = "", ylab = "", ...) {
    
    plot(px[,1:2], asp = T, 
         xlim = xlim, ylim = ylim, pch = pch, main = main, xlab = xlab, ylab = ylab, cex = cex, ...)
    
    if (panels) draw.panels()
}


####################################################################################################


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





#' Colour scheme for bad pixel plot
#'
#' Support function: returns a vector of colours to be used in bad pixel plotting.
#' @param block Vector of categories to omit from the plot, which will be set to colour NA. Default is c("edge", "s.bright").
#' @return Vector of colours.
#' @export
#' 
bp.colours <- function(block = c("edge", "s.bright")) {
    
    cats <- c("no response", "dead", "hot", "v.bright","bright", "s.bright", "screen spot", "edge", "v.dim", "dim", "s.dim")
    bp.cols <- c("purple", "black", "red", "orange", "gold", "yellow", "grey", "darkgrey", "green3", "green", "lightskyblue")
    bp.cols[which(cats %in% block)] <- NA
    return(bp.cols)
}


#' Plot bad pixels
#' 
#' Plot the coordinates and types of a bad pixel map
#' @param bpm Data frame containing xy coordinates of pixels to be plotted, and a field "type" specifying the category of each bad pixel identified.
#' @param block Vector of categories to omit from the plot, which will be set to colour NA. Default is c("edge", "s.bright").
#' @param ... Additional optional graphical parameters to pass to plotting function
#' @export
#' 
plot.bad.px <- function(bpm, pch = 15, xlab = "", ylab = "", block = c("edge", "s.bright"), ...) {
    plot(bpm[,1:2], pch = pch, col = bp.colours(block = block)[bpm$type], asp = T, xlab = xlab, ylab = ylab, ...)
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


#' Get coordinates of points to plot zoom
#' 
#' Given a set of points and the surrounding area to plot, obtain coordinates to create a plot of a small subset of an image
#' @param points Matrix or data frame, the first two columns of which are the x and y coordinates of the points to be plotted.
#' @param surround Number of additional pixels to plot in every direction around the selected points
#' @return matrix of coordinates to plot
#' @export
#'  
get.focus <- function(points, surround = 3) {
    n <- 2 * surround + 1
    focus <- matrix(c(round(mean(points$x), 0) + rep(c(-surround: surround), n),
                      round(mean(points$y), 0) + sort(rep(c(-surround: surround), n))), ncol = 2)
    focus[focus <= 0] <- 0
    focus[focus >= 2048] <- 2048
    return(focus)
}


#' Histogram plot with SD levels shown
#' 
#' Produce a histogram with coloured bars at the bottom showing scale used. Defaults can be used to show legend for plots produced using \link{\code{pixel.image}}.
#' @param data Array of numbers to be plotting in histogram
#' @param scale Vector of scale cutpoints. Default is \link{\code{sd.levels}}.
#' @param scale.colours Vector of scale colours. Default is \link{\code{sd.colours}}
#' @param xlim Vector of lower and upper x-limits. If not provided, will use central 0.95 of normal distribution with mean & sd of observed data.
#' @export
#' @examples
#' pixel.image(b.150828)
#' scale.hist(b.150828)
s.hist <- function(data, scale = sd.levels(data), scale.colours = sd.colours(), xlim, scale.height = 600, ...) {
    
    data <- data[!is.na(data)]
    
    if (missing(xlim)) {
        xlim <- c(floor(qnorm(0.025, mean(data), sd(data))/10)*10,
                  ceiling(qnorm(0.975, mean(data), sd(data))/10)*10)
    }
    
    hist(data, breaks = "fd", xlim = xlim, ...)
    
    # add colours to indicate scale of pixel map
    cl <- cut(xlim[1]:xlim[2], scale)
    points(xlim[1]:xlim[2], rep(-scale.height, length(xlim[1]:xlim[2])), pch = 15, col = scale.colours[cl])
}


#' Plot focal area
#' 
#' Create pixel image of area surrounding a point of interest. Pixel values are labelled and bad pixels identified highlighted
#' @param im 2d image array to plot
#' @param centre Vector of x and y coordinates of pixel of particular interest
#' @param surround Integer: set size of area surrounding pixel of interest to be displayed. Default is 5, which displays an 11x11 square centred on the pixel of interest.
#' @param dp Integer: display pixel values to how many decimal places? Default is 1.
#' @param scale.by Integer: divide pixel values by this number for easier display and comparison. Default is 1000, so a value of 65535 will be displayed as 65.5
#' @param lbl.cex Magnification factor to apply to labels in each cell. Default is 0.7
#' @param bad.px Optional bad pixel map to be used to highlight any pixels in the display area that have been identified as defective.
#' @param bpx.cex Magnification factor to apply to symbol used to highlight bad pixels. Default is 2.5
#' @param ... Additional optional graphical arguments
#' @export
#' 
focal.plot <- function(im, centre, surround = 5, dp = 1, scale.by = 1000, lbl.cex = 0.7, bad.px, bpx.cex = 2.5, ...) {
    
    ff <- get.focus(data.frame(x = centre[1], y = centre[2]), surround = surround)
    
    pixel.image(im, xlim = range(ff[,1]), ylim = range(ff[,2]), ...)
    text(ff, labels = round(im[ff]/scale.by, dp), cex = lbl.cex)
    
    if (!missing(bad.px)) {
        points(bad.px[,1:2], pch = 0, cex = bpx.cex)
    }
}


#' Draw subpanels on 2048 x 2048 array
#' 
#' Marks locations of 32 subpanels of 128 x 1024 pixels on plot
#' @param ... Optional graphical arguments
#' @export
#' 
draw.panels.2048 <- function(...) {
    
    # horizontal midline
    lines (c(0, 2048), c(1024.5, 1024.5), ...)
    
    for (i in 1:15) {
        lines(c(128 * i + 0.5, 128 * i + 0.5), c(0, 2048), ...)
    }
}


#' Colour ramp for pixel images
#' 
#' @export
sd.ramp <- function() {
    colorRampPalette(c("darkblue", "cyan3", "green2", "lemonchiffon", "gold", "red3", "purple"), 
                     space = "Lab")
}