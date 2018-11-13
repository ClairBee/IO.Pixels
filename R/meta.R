
#' @export
refresh <- function() {
    require(devtools)
    require(roxygen2)
    org.dir <- getwd()
    target.dir <- paste0("/home/clair/R/My-packages/IO.Pixels")
    
    if (org.dir != target.dir) {setwd(target.dir)}
    document()
    setwd("..")
    install("IO.Pixels")
    setwd(org.dir)
}