
#' Get transition matrices for bad pixel maps
#' 
#' Given a list of bad pixel maps in coordinate form, return a list of transition matrices showing movements between states observed between one acquisition and the next.
#' @param bpx List of dataframes containing bad pixel coordinates and types
#' @return List of tables showing movement between states from one acquisition to the next
#' @export
#' 
get.transitions <- function(bp) {
    tr <- list()
    
    lvls <- c("normal", levels(bp[[1]]$type))
    
    for (i in 1:(length(bp) - 1)) {
        tr[[i]] <- table("From" = ordered(lvls[c(bpx2im(bp[[i]])) + 1], levels = lvls),
                         "To" = ordered(lvls[c(bpx2im(bp[[i+1]])) + 1], levels = lvls))
        names(tr)[[i]] <- paste(names(bp)[c(i,(i+1))], collapse = "-")
    }
    tr
}


