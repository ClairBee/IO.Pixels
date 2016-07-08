
#' Basic bad pixel map
#' 
#' @export
#' 
basic.bpx <- function(acq, screen.spots = T) {

    # calculate thresholds
    th <- apply(acq[,,c("black", "grey")], 3, function(im) {
        med <- median(im, na.rm = T)
        c(v.dim = med * 0.5, dim = med * 0.75,
          bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
    })
    
    px <- list(edge = edge.px(acq, edge.width = 40),
               no.resp = no.response(bright.image = acq[,,"grey"], dark.image = acq[,,"black"]),
               hot = which(acq[,,"black"] == 65535, arr.ind = T),
               dead = which(acq[,,"grey"] == 0, arr.ind = T),
               line = which(find.lines(acq[, , "black"], midline = 1024.5) > 0, arr.ind = T),
               v.bright = which(acq[, , "black"] > th["v.bright", "black"] | 
                                    acq[, , "grey"] > th["v.bright", "grey"], arr.ind = T),
               bright = which(acq[, , "black"] > th["bright", "black"] |
                                  acq[, , "grey"] > th["bright", "grey"], arr.ind = T),
               v.dim = which(acq[, , "black"] < th["v.dim", "black"] | 
                                 acq[, , "grey"] < th["v.dim", "grey"], arr.ind = T), 
               dim = which(acq[, , "black"] < th["dim", "black"] | 
                               acq[, , "grey"] < th["dim", "grey"], arr.ind = T))
    
    if (screen.spots) {
        px$screen.spot <- screen.spots(acq[,,"white"], enlarge = T, ignore.edges = 40)
    }
    
    # remove empty types & combine into data frame
    px <- px[sapply(px, length) > 0]
    px <- rbind.fill(lapply(names(px), function(df) data.frame(px[[df]], type = toString(df))))
    
    Cat <- c("no.resp", "dead", "hot", "v.bright", "bright", "line.b", "edge", "line.d", "screen.spot", "v.dim", "dim")
    px$type <- ordered(px$type, levels = Cat)
    
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    return(px)
}



#' Linearity-based pixel map
#' 
#' @export
#' 
nl.bpx <- function(acq, w.th = 1000) {
 
       # calculate thresholds
    th <- apply(acq[,,c("black", "grey")], 3, function(im) {
        med <- median(im, na.rm = T)
        c(v.dim = med * 0.5, dim = med * 0.75,
          bright = med + (65535 - med) / 4, v.bright = med + (65535 - med) / 2)
    })
    
    wlm <- fit.w.lm(acq)$df
    
    px <- rbind(data.frame(which(acq[,,"black"] > th["bright", "black"], arr.ind = T), type = "warm"),
                data.frame(which(acq[,,"grey"] > th["bright", "grey"], arr.ind = T), type = "bright"),
                data.frame(setNames(wlm[abs(wlm$res) > w.th, 1:2], nm = c("row", "col")), type = "nonlinear"),
                data.frame(which(acq[,,"grey"] - acq[,,"black"] < 500 & acq[,,"grey"] > acq[,,"black"], arr.ind = T), type = "stuck"),
                data.frame(which(acq[,,"black"] < th["dim", "black"], arr.ind = T), type = "cool"),
                data.frame(which(acq[,,"grey"] < th["dim", "grey"], arr.ind = T), type = "dim"))
    
    px$type <- ordered(px$type, levels = c("stuck", "warm", "bright", "nonlinear", "cool", "dim"))
    px <- px[order(px$type),]
    px <- px[!duplicated(px[,1:2]),]
    
    return(px)
}


#' To add: local pixel map
#' 