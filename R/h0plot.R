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
      annotate(geom = "text", x = statistic, y = Inf, hjust = -0.15, vjust = 1.5, 
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
  gg <- gg + ggtitle("Permutation Distribution of Test Statistic under H0")
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
  gg <- gg + ggtitle("Bootstrap Distribution of Test Statistic under H0")
  gg
}
h0plot.htest <- function(x, sig.level = 0.05, hist.alpha = 0.2, 
                         stat.col = "darkgreen", dens.col = "black",
                         dens.alpha = 0.5, cval.col = "darkred", lwd = 1, 
                         breaks = NULL, digits = 3, qtail = 1e-3, n = 501, ...){
  if(names(x$statistic) == "t"){
    dfun <- function(x){ }
    body(dfun) <- substitute({ dt(x, df = para) },
                             list(para = x$parameter))
    if(x$alternative == "two.sided"){
      cval <- qt(p = c(sig.level/2, 1-sig.level/2), df = x$parameter)
    }else if(x$alternative == "less"){
      cval <- qt(p = sig.level, df = x$parameter)
    }else{
      cval <- qt(p = 1-sig.level, df = x$parameter)
    }
    MIN <- min(-abs(x$statistic), qt(qtail, df = x$parameter))
    MAX <- max(abs(x$statistic), qt(1-qtail, df = x$parameter))
    xlim <- c(MIN, MAX)
  }
  if(names(x$statistic) == "gt"){
    dfun <- function(x){ }
    body(dfun) <- substitute({ dgt(x, n1, n2, sd1, sd2) },
                             list(n1 = x$n[1], n2 = x$n[2], 
                                  sd1 = x$estimate[3], 
                                  sd2 = x$estimate[4]))
    if(x$alternative == "two.sided"){
      cval <- c(qgt(p = sig.level/2, n1 = x$n[1], n2 = x$n[2], 
                    sd1 = x$estimate[3], sd2 = x$estimate[4]), 
                qgt(p = 1-sig.level/2, n1 = x$n[1], n2 = x$n[2], 
                    sd1 = x$estimate[3], sd2 = x$estimate[4]))
    }else if(x$alternative == "less"){
      cval <- qgt(p = sig.level, n1 = x$n[1], n2 = x$n[2], 
                  sd1 = x$estimate[3], sd2 = x$estimate[4])
    }else{
      cval <- qgt(p = 1-sig.level, n1 = x$n[1], n2 = x$n[2], 
                  sd1 = x$estimate[3], sd2 = x$estimate[4])
    }
    MIN <- min(-abs(x$statistic), qgt(qtail, n1 = x$n[1], n2 = x$n[2], 
                                      sd1 = x$estimate[3], sd2 = x$estimate[4]))
    MAX <- max(abs(x$statistic), qgt(1-qtail, n1 = x$n[1], n2 = x$n[2], 
                                     sd1 = x$estimate[3], sd2 = x$estimate[4]))
    xlim <- c(MIN, MAX)
  }
  if(names(x$statistic) == "F"){
    dfun <- function(x){ }
    body(dfun) <- substitute({ df(x, df1 = para1, df2 = para2) },
                             list(para1 = x$parameter["num df"],
                                  para2 = x$parameter["denom df"]))
    if("alternative" %in% names(x)){
      if(x$alternative == "two.sided"){
        cval <- qf(p = c(sig.level/2, 1-sig.level/2), df1 = x$parameter["num df"],
                   df2 = x$parameter["denom df"])
      }else if(x$alternative == "less"){
        cval <- qf(p = sig.level, df1 = x$parameter["num df"],
                   df2 = x$parameter["denom df"])
      }else{
        cval <- qf(p = 1-sig.level, df1 = x$parameter["num df"],
                   df2 = x$parameter["denom df"])
      }
    }else{
      cval <- qf(p = 1-sig.level, df1 = x$parameter["num df"],
                 df2 = x$parameter["denom df"])
    }
    MIN <- 0
    MAX <- max(abs(x$statistic), qf(1-qtail, df1 = x$parameter["num df"],
                                    df2 = x$parameter["denom df"]))
    xlim <- c(MIN, MAX)
  }
  if(names(x$statistic) == "X-squared"){
    dfun <- function(x){ }
    body(dfun) <- substitute({ dchisq(x, df = para) },
                             list(para = x$parameter))
    if("alternative" %in% names(x)){
      if(x$alternative == "two.sided"){
        cval <- qchisq(p = c(sig.level/2, 1-sig.level/2), df = x$parameter)
      }else if(x$alternative == "less"){
        cval <- qchisq(p = sig.level, df = x$parameter)
      }else{
        cval <- qchisq(p = 1-sig.level, df = x$parameter)
      }
    }else{
      cval <- qchisq(p = 1-sig.level, df = x$parameter)
    }
    MIN <- ifelse(x$parameter == 1, qchisq(qtail, df = x$parameter), 0)
    MAX <- max(abs(x$statistic), qchisq(1-qtail, df = x$parameter))
    xlim <- c(MIN, MAX)
  }
  if(!names(x$statistic) %in% c("t", "gt", "F", "X-squared")){
    stop("Not yet implemented!")
  }
  gg <- ggplot(data = data.frame(x = xlim), aes(x)) +
    stat_function(fun = dfun, n = n) + ylab("density") + xlab("test statistic") +
    annotate(geom = "text", x = x$statistic, y = Inf, hjust = -0.15, vjust = 1.5, 
             label = paste0("p value = ", signif(x$p.value, digits = digits)),
             fontface = 2, color = stat.col) +
    ggtitle(paste0("Distribution of Test Statistic under H0\n",
                   x$method))
  
  if("alternative" %in% names(x)){
    if(x$alternative == "two.sided"){
      gg <- gg + stat_function(fun = dfun, xlim = c(MIN, cval[1]), geom = "area",
                               fill = cval.col, alpha = dens.alpha) + 
        stat_function(fun = dfun, xlim = c(cval[2], MAX), geom = "area",
                      fill = cval.col, alpha = dens.alpha)
    }else if(x$alternative == "less"){
      gg <- gg + stat_function(fun = dfun, xlim = c(MIN, cval), geom = "area",
                               fill = cval.col, alpha = dens.alpha)
    }else{
      gg <- gg + stat_function(fun = dfun, xlim = c(cval, MAX), geom = "area",
                               fill = cval.col, alpha = dens.alpha)
    }
  }else{
    gg <- gg + stat_function(fun = dfun, xlim = c(cval, MAX), geom = "area",
                             fill = cval.col, alpha = dens.alpha)
  }
  
  gg <- gg + geom_vline(xintercept = x$statistic, color = stat.col, lwd = lwd)
  gg
}
