
#' Get bad pixel map for a given date
#' 
#' Use thresholding and absolute values to classify bad pixels in a single acquisition, using the thresholding limits set in \link{\code{bad.bx.limits}}.
#' @details Requires that the pixelwise mean array \code{pw.m} and a residual array \code{res} be available in the global environment.
#' @param dt Integer or string giving acquisition date (format yymmdd)
#' @export
#' @examples
#' bp <- bad.pixels("160314")
#' 
bad.px <- function(dt) {
    dt <- toString(dt)
    
    qq <- screen.spots.xy(dt)
    
    bp <- rbind(setNames(data.frame(matrix(c(sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996)), rep(1:1996, 20),
                                             rep(1:1996, 20), sort(rep(1:10, 1996)), sort(rep(1987:1996, 1996))), ncol = 2), 
                                    ordered("edge", levels = c("no response", "dead", "hot", "v.bright", "bright", "s.bright", "screen spot", "edge", "v.dim", "dim", "s.dim"))),
                         c("row", "col", "type")),
                data.frame(no.response(dt), type = "no response"),
                data.frame(which(pw.m[, , "black", dt] == 65535, arr.ind = T), type = "hot"),
                data.frame(which(pw.m[, , "white", dt] == 0, arr.ind = T), type = "dead"),
                screen.spots.xy(dt),
                get.dim.bright.px(res[, , "grey", dt]),
                get.dim.bright.px(res[, , "white", dt]))
    
    bp <- bp[order(bp$type),]
    bp[!duplicated(bp[,1:2]),]
}