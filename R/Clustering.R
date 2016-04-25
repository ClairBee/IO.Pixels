
#' Identify clusters and superclusters of pixels
#' 
#' Identify clusters of adjacent pixels of each type in a bad pixel map, and also superclusters of adjacent hot and bright pixels.
#' @param bpx Bad pixel map: a matrix or data frame containing coordinates of bad pixels and column "type" defining bad pixel type ("hot", "bright", "dim", "dead").
#' @return Data frame containing bad pixel coordinates and type, along with ID and size of clusters and superclusters to which each is assigned, and a shape classification based on these.
#' @export
#' @examples
#' zz <- cluster.px(bp)
cluster.px <- function(bpx) {
    
    #--------------------------------------------------------------------------------------------
    # IDENTIFY CLUSTERS OF PIXELS
    
    # get unique coords from bad pixel map
    px <- bpx[!duplicated(bpx[,1:2]),]
    
    # remove any pixels identified as OK (can occur when working with subset of images/SD classifier)
    px <- px[px$type != "-",]
    
    # create matrix to populate raster
    px.vals <- array(dim = c(1996, 1996))
    px.vals[as.matrix(px[,1:2])] <- px$type
    
    # convert to raster
    r <- raster(t(px.vals[,1996:1]),  xmn = 0.5, xmx = 1996.5, ymn = 0.5, ymx = 1996.5)
    
    # cluster bad pixels by type
    xy <- data.frame()
    for (val in unique(getValues(r)[!is.na(getValues(r))])) {
        
        # mask all but one type of bad pixel
        tmp <- r; 
        tmp[tmp != val] <- NA
        
        # clump remaining bad pixels together
        cc <- clump(tmp, dir = 4)
        
        xy <- rbind(xy,
                    merge(data.frame(xyFromCell(cc, which(!is.na(getValues(cc)))),
                                     type = val,
                                     id = getValues(cc)[!is.na(getValues(cc))]),
                          setNames(data.frame(val, freq(cc)), c("type", "id", "count"))))
    }
    xy$type <- factor(levels(px$type)[xy$type], levels = levels(px$type)[levels(px$type) != "-"])
    
    # rationalize cluster IDs (otherwise IDs are reused for each type)
    xy <- transform(xy, id = match(apply(xy[, c("type", "id")], 1, paste, collapse = "-"),
                                   unique(apply(xy[, c("type", "id")], 1, paste, collapse = "-"))))

    #--------------------------------------------------------------------------------------------
    # IDENTIFY SUPERCLUSTERS OF SUBSET OF PIXELS
    
    # create raster containing only hot/bright clusters and singletons
    sc <- r
    sc[cellFromXY(sc, xy = xy[(xy$count > 9) | (xy$type %in% c("dead", "dim")), c("x", "y")])] <- NA
    
    # clump these smaller clusters of pixels
    sc.cc <- clump(sc, dir = 4)
    sc.xy <- merge(data.frame(xyFromCell(sc.cc, which(!is.na(getValues(sc.cc)))),
                              sc.id = getValues(sc.cc)[!is.na(getValues(sc.cc))]),
                   setNames(data.frame(freq(sc.cc)), c("sc.id", "sc.count")))
    
    # remove superclusters of size 1
    sc.xy <- sc.xy[sc.xy$sc.count > 1,]

    xy <- merge(xy, sc.xy, all = T)
    
    # classify shape according to cluster size/type
    xy$shape <- factor("Other", levels = c("Singleton", "Cluster", "Supercluster", "Other"))
        xy$shape[xy$count <= 9] <- "Cluster"
        xy$shape[xy$count == 1] <- "Singleton"
        xy$shape[xy$sc.count > xy$count] <- "Supercluster"
        xy$shape[xy$sc.count > 1 & xy$type == "hot"] <- "Supercluster"
        
    # rationalise supercluster IDs
    xy$sc.id[xy$shape != "Supercluster"] <- NA
    xy <- transform(xy, sc.id = match(sc.id, unique(sc.id)))
    
    # return data frame of pixel coords, pixel type, cluster size & id and supercluster id
    return(xy[,c("x", "y", "type", "shape", "id", "count", "sc.id", "sc.count")])
}