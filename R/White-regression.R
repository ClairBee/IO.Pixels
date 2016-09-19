
#' Linear regression of observed pixel values
#' 
#' Perform linear regression of specified variables
#' @param im Acquisition array on which regression is to be performed. Must contain black, white and grey images.
#' @param terms Text string providing formula to be fitted. Default is "g ~ b * w", predicting grey value based on black and white.
#' @param res.only Boolean: return only a matrix of residuals (T) or (F) a list containing the fitted data, residuals and original data, along with the RMSE and adjusted RMSE.
#' @return If res.only == T, returns a matrix of residuals. If F, returns a list containing a data.frame and the RMSE of the fitted model.
#' @export
#' 
fit.w.lm <- function(im, terms = "g ~ b * w", midline = 1024.5, res.only = T) {
    
    df <- setNames(data.frame(melt(im[, , "black"]), 
                              melt(im[, , "grey"]), 
                              melt(im[, , "white"]))[, c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    if (!is.na(midline)) {
        df$upper <- df$y > midline
        terms <- paste0(terms, " + upper")
    }
    
    # fit linear model to central part of image only (excludes edge effects)
    w.lm <- lm(as.formula(terms), 
               data = df[findInterval(df$x, c(40.5, 2008.5)) == 1 & 
                             findInterval(df$y, c(40.5, 2008.5)) == 1, ])
    
    df$fv <- predict(w.lm, df)
    
    target <- gsub(" ~.*$", "", terms)
    
    df$res <- eval(parse(text = paste0("df$", target))) - df$fv
    
    if (res.only) {
        return(array(df$res, dim = dim(im[,,"black"])))
    } else {
        return(list(df = df, r2 = round(summary(w.lm)$adj.r.squared, 3), rmse = round(summary(w.lm)$sigma, 2)))
    }
}


#' Scatter plot of linear regression
#' 
#' @export
#' 
wlm.plot <- function(dt) {
    
    dt <- toString(dt)
    
    wlm <- eval(parse(text = paste0("lm.", dt)))
    
    pdf(paste0(fpath, "white-fit-", dt, ".pdf"))
    par(mar = c(4, 4, 1, 1))
    smoothScatter(wlm$df$w, wlm$df$fv, nrpoints = 0, xlim = c(0,65535),
                  colramp = colorRampPalette(c(NA, "gold", "red", "blue"), space = "Lab"),
                  xlab = "Observed", ylab = "Fitted value")
    abline(line(wlm$df$w, wlm$df$fv), col = adjustcolor("darkred", alpha = 0.4), lty = 2)
    dev.off()
    
    write(paste0("Adj $r^2$ ", round(wlm$r2, 3), "; ",
                 "RMSE ", round(wlm$rmse, 2)),
          paste0(fpath, "fitted-wv-all-", dt, ".txt"))
}
