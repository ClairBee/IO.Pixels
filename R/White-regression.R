
#' Linear regression of white value against grey & black
#' 
#' @export
#' 
fit.w.lm <- function(im) {
    
    # fit linear model excluding edge pixels
    df <- setNames(data.frame(melt(im[,,"black"]), 
                              melt(im[,,"grey"]),
                              melt(im[,,"white"]))[,c("X1", "X2", "value", "value.1", "value.2")],
                   nm = c("x", "y", "b", "g", "w"))
    w.lm <- lm(w ~ b * g,
               data = df[findInterval(df$x, c(40.5, 2008.5)) == 1 &
                             findInterval(df$y, c(40.5, 2008.5)) == 1,])
    
    df$fv <- predict(w.lm, df[,c("b", "g")])
    df$res <- df$w - df$fv
    
    list(df = df, r2 = round(summary(w.lm)$adj.r.squared, 3), rmse = round(summary(w.lm)$sigma, 2))
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
