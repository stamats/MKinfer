h0plot <- function(x, ...){ 
  UseMethod("h0plot")
}
h0plot.default <- function(x, statistic, pval, alternative = "two.sided", 
                           sig.level = 0.05, hist.alpha = 0.2, 
                           stat.col = "darkgreen", dens.col = "black",
                           dens.alpha = 0.5, cval.col = "darkred", lwd = 1, 
                           breaks = NULL, digits = 3, ...){
  stopifnot(is.numeric(statistic))
  stopifnot(is.numeric(pval))
  stopifnot(is.numeric(sig.level))
  
  if(is.numeric(x)){
    dens <- density(x)
    if(alternative == "two.sided"){
      cval <- quantile(x, probs = c(sig.level/2, 1-sig.level/2))
      DFcval.lo <- data.frame(x = dens$x[dens$x < cval[1]], 
                              y = dens$y[dens$x < cval[1]])
      DFcval.up <- data.frame(x = dens$x[dens$x > cval[2]],
                              y = dens$y[dens$x > cval[2]])
    }else if(alternative == "less"){
      cval <- quantile(x, probs = sig.level)
      DFcval <- data.frame(x = dens$x[dens$x < cval],
                           y = dens$y[dens$x < cval])
    }else{
      cval <- quantile(x, probs = 1-sig.level)
      DFcval <- data.frame(x = dens$x[dens$x > cval],
                           y = dens$y[dens$x > cval])
    }
    
    DF <- data.frame(stat = x)
    gg <- ggplot(DF, aes(x = .data$stat)) +
      geom_histogram(aes(y = after_stat(density)), breaks = breaks, alpha = hist.alpha) +
      geom_density(color = dens.col, lwd = lwd) + ylab("density") +
      annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5, 
               label = paste0("p value = ", signif(pval, digits = digits)),
               fontface = 2, color = stat.col) +
      xlab("test statistic")
    if(alternative == "two.sided"){
      gg <- gg + geom_area(aes(x = x, y = .data$y), data = DFcval.lo, fill = cval.col, alpha = dens.alpha) +
        geom_area(aes(x = x, y = .data$y), data = DFcval.up, fill = cval.col, alpha = dens.alpha)
    }else{
      gg <- gg + geom_area(aes(x = x, y = .data$y), data = DFcval, fill = cval.col, alpha = dens.alpha)
    }
  }else{
    stop("not yet implemented!")
  }
  gg <- gg + geom_vline(xintercept = statistic, color = stat.col, lwd = lwd)
  gg
}
h0plot.perm.htest <- function(x, sig.level = 0.05, hist.alpha = 0.2, 
                              stat.col = "darkgreen", dens.col = "black",
                              dens.alpha = 0.5, cval.col = "darkred", lwd = 1, 
                              breaks = NULL, digits = 3, ...){
  gg <- h0plot(x = x$perm.statistic, statistic = x$statistic, 
               pval = x$perm.p.value, alternative = x$alternative, 
               sig.level = sig.level, hist.alpha = hist.alpha, 
               stat.col = stat.col, dens.col = dens.col,
               dens.alpha = dens.alpha, cval.col = cval.col, lwd = lwd, 
               breaks = breaks, digits = digits, ...)
  gg <- gg + ggtitle("Permutation Distribution")
  gg
}
h0plot.boot.htest <- function(x, sig.level = 0.05, hist.alpha = 0.2, 
                              stat.col = "darkgreen", dens.col = "black",
                              dens.alpha = 0.5, cval.col = "darkred", lwd = 1, 
                              breaks = NULL, digits = 3, ...){
  gg <- h0plot(x = x$boot.statistic, statistic = x$statistic, 
               pval = x$boot.p.value, alternative = x$alternative, 
               sig.level = sig.level, hist.alpha = hist.alpha, 
               stat.col = stat.col, dens.col = dens.col,
               dens.alpha = dens.alpha, cval.col = cval.col, lwd = lwd, 
               breaks = breaks, digits = digits, ...)
  gg <- gg + ggtitle("Bootstrap Distribution")
  gg
}
