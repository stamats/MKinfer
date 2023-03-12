baplot <- function(x, y, loa.level = 0.95, conf.level = 0.95, 
                   xlab = "Difference", ylab = "Mean", 
                   title = "Bland-Altman Plot", xlim = NULL, ylim = NULL,
                   type = c("parametric", "nonparametric"), 
                   loa.type = c("unbiased", "biased"), ci.diff = TRUE, ci.loa = TRUE, 
                   ci.type = c("exact", "approximate", "boot"), 
                   bootci.type = NULL, R = 9999, print.res = TRUE,
                   color.low = "#4575B4", color.upp = "#D73027", 
                   alpha = 1, shape = 19, na.rm = TRUE, ...){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(length(x) == length(y))
  if(length(loa.level) != 1)
    stop("'loa.level' has to be of length 1 (limit of agreement level)")
  if(loa.level < 0.5 | loa.level > 1)
    stop("'loa.level' has to be in [0.5, 1]")
  if(length(conf.level) != 1)
    stop("'conf.level' has to be of length 1 (confidence level)")
  if(conf.level < 0.5 | conf.level > 1)
    stop("'conf.level' has to be in [0.5, 1]")
  loa.alpha <- 1 - loa.level
  conf.alpha <- 1 - conf.level
  
  type <- match.arg(type)
  ci.type <- match.arg(ci.type)
  loa.type <- match.arg(loa.type)
  if(na.rm){
    ind.na <- is.na(x) | is.na(y)
    x <- x[!ind.na]
    y <- y[!ind.na]
  }
  if(type == "parametric" && is.null(bootci.type)){
    bootci.type <- "stud"
  }
  if(type == "nonparametric" && is.null(bootci.type)){
    bootci.type <- "perc"
  }
  stopifnot(length(bootci.type) == 1)
  stopifnot(bootci.type %in% c("stud", "perc", "bca"))
  
  Mean <- 0.5*(x+y)
  Diff <- x-y
  DF <- data.frame(Mean, Diff)
  
  if(type == "parametric"){
    m.diff <- mean(DF$Diff)
    sd.diff <- sd(DF$Diff)
    k <- qnorm(1-loa.alpha/2)
    n <- nrow(DF)
    df <- n - 1
    if(loa.type == "unbiased"){
      logc <- log(sqrt(df/2)) + lgamma(df/2) - lgamma(n/2)
      sd.cor <- exp(logc)
      loa <- m.diff + c(-k, k)*sd.cor*sd.diff
    }
    if(loa.type == "biased"){
      loa <- m.diff + c(-k, k)*sd.diff
    }
  }
  if(type == "nonparametric"){
    m.diff <- median(DF$Diff)
    loa <- quantile(DF$Diff, probs = c(loa.alpha/2, 1-loa.alpha/2))
  }
  gg <- ggplot(DF, aes(x = Mean, y = Diff)) +
    geom_point(alpha = alpha, shape = shape) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = m.diff) +
    geom_hline(yintercept = loa, color = c(color.low, color.upp)) +
    xlab(xlab) + ylab(ylab) + ggtitle(title)
  
  if(ci.diff){
    if(ci.type == "exact"){
      if(type == "parametric"){
        m.diff.ci <- meanCI(DF$Diff, conf.level = conf.level)$conf.int[1,]
      }
      if(type == "nonparametric"){
        m.diff.ci <- medianCI(DF$Diff, conf.level = conf.level)$conf.int[1,]
      }
    }
    if(ci.type == "approximate"){
      if(type == "parametric"){
        m.diff.ci <- meanCI(DF$Diff, conf.level = conf.level)$conf.int[1,]
      }
      if(type == "nonparametric"){
        m.diff.ci <- medianCI(DF$Diff, method = "asymptotic", 
                              conf.level = conf.level)$conf.int[1,]
      }
    }
    if(ci.type == "boot"){
      if(type == "parametric"){
        m.diff.ci <- meanCI(DF$Diff, conf.level = conf.level, boot = TRUE, 
                            bootci.type = bootci.type, R = R, ...)$conf.int[[4]][4:5]
      }
      if(type == "nonparametric"){
        m.diff.ci <- medianCI(DF$Diff, conf.level = conf.level, method = "boot", 
                              bootci.type = bootci.type, R = R, ...)$conf.int[[4]][4:5]
      }
    }
    gg <- gg + geom_hline(yintercept = m.diff.ci, linetype = "dotted")
  }
  if(ci.loa){
    if(ci.type == "exact"){
      if(type == "parametric"){
        zp <- qnorm(loa.alpha/2)
        tl <- qt(conf.alpha/2, df, zp*sqrt(n))
        tu <- qt(1-conf.alpha/2, df, zp*sqrt(n))
        loa.low.ci <- m.diff + c(tl, tu)*sd.diff/sqrt(n)
        zp <- qnorm(1-loa.alpha/2)
        tl <- qt(conf.alpha/2, df, zp*sqrt(n))
        tu <- qt(1-conf.alpha/2, df, zp*sqrt(n))
        loa.upp.ci <- m.diff + c(tl, tu)*sd.diff/sqrt(n)
      }
      if(type == "nonparametric"){
        loa.low.ci <- quantileCI(DF$Diff, prob = loa.alpha/2, 
                                 conf.level = conf.level)$conf.int[1,]
        loa.upp.ci <- quantileCI(DF$Diff, prob = 1-loa.alpha/2,
                                 conf.level = conf.level)$conf.int[1,]
      }
    }
    if(ci.type == "approximate"){
      if(type == "parametric"){
        n <- nrow(DF)
        sd.loa <- sqrt(1/n + k^2/(2*(n-1)))*sd.diff
        k.loa <- qt(1-conf.alpha/2, df = n-1)
        loa.low.ci <- loa[1] + c(-k.loa, k.loa)*sd.loa
        loa.upp.ci <- loa[2] + c(-k.loa, k.loa)*sd.loa
      }
      if(type == "nonparametric"){
        loa.low.ci <- quantileCI(DF$Diff, prob = loa.alpha/2, conf.level = conf.level, 
                                 method = "asymptotic")$conf.int[1,]
        loa.upp.ci <- quantileCI(DF$Diff, prob = 1-loa.alpha/2, conf.level = conf.level, 
                                 method = "asymptotic")$conf.int[1,]
      }
    }
    if(ci.type == "boot"){
      if(type == "parametric"){
        boot.loa.low <- function(x, i){
          n <- length(x[i])
          df <- n-1
          logc <- log(sqrt(df/2)) + lgamma(df/2) - lgamma(n/2)
          sd.cor <- exp(logc)
          k <- qnorm(loa.alpha/2)
          m.diff <- mean(x[i])
          sd.diff <- sd(x[i])
          loa <- m.diff + k*sd.cor*sd.diff
          var.loa <- (1 + n*k^2*(sd.cor^2-1))*sd.diff^2/n
          c(loa, var.loa)
        }
        boot.loa.upp <- function(x, i){
          n <- length(x[i])
          df <- n-1
          logc <- log(sqrt(df/2)) + lgamma(df/2) - lgamma(n/2)
          sd.cor <- exp(logc)
          k <- qnorm(1-loa.alpha/2)
          m.diff <- mean(x[i])
          sd.diff <- sd(x[i])
          loa <- m.diff + k*sd.cor*sd.diff
          var.loa <- (1 + n*k^2*(sd.cor^2-1))*sd.diff^2/n
          c(loa, var.loa)
        }
        boot.out.low <- boot(DF$Diff, statistic = boot.loa.low, R = R, ...)
        boot.out.upp <- boot(DF$Diff, statistic = boot.loa.upp, R = R, ...)
        loa.low.ci <- boot.ci(boot.out.low, type = bootci.type, conf = conf.level)[[4]][4:5]
        loa.upp.ci <- boot.ci(boot.out.upp, type = bootci.type, conf = conf.level)[[4]][4:5]
      }
      if(type == "nonparametric"){
        loa.low.ci <- quantileCI(DF$Diff, prob = loa.alpha/2, conf.level = conf.level, 
                                 method = "boot", bootci.type = bootci.type, R = R, ...)$conf.int[[4]][4:5]
        loa.upp.ci <- quantileCI(DF$Diff, prob = 1-loa.alpha/2, conf.level = conf.level, 
                                 method = "boot", bootci.type = bootci.type, R = R, ...)$conf.int[[4]][4:5]
      }
    }
    gg <- gg + geom_hline(yintercept = loa.low.ci, linetype = "dotted", color = color.low) +
      geom_hline(yintercept = loa.upp.ci, linetype = "dotted", color = color.upp)
  }
  if(!is.null(xlim)) gg <- gg + xlim(xlim)
  if(!is.null(ylim)) gg <- gg + ylim(ylim)
  
  if(print.res){
    if(ci.diff){
      names(m.diff.ci) <- c(paste(conf.alpha/2*100, "%"), paste((1-conf.alpha/2)*100, "%"))
      if(type == "parametric"){
        res <- list("mean of differences" = c("estimate" = m.diff, m.diff.ci))
      }
      if(type == "nonparametric"){
        res <- list("median of differences" = c("estimate" = m.diff, m.diff.ci))
      }
    }else{
      if(type == "parametric"){
        res <- c("mean of differences" = m.diff)
      }
      if(type == "nonparametric"){
        res <- c("median of differences" = m.diff)
      }
    }
    if(ci.loa){
      names(loa.low.ci) <- c(paste(conf.alpha/2*100, "%"), paste((1-conf.alpha/2)*100, "%"))
      names(loa.upp.ci) <- c(paste(conf.alpha/2*100, "%"), paste((1-conf.alpha/2)*100, "%"))
      res <- c(res, list(c("estimate" = loa[1], loa.low.ci)), 
               list(c("estimate" = loa[2], loa.upp.ci)))
    }else{
      res <- c(res, loa)
    }
    names(res)[2] <- paste("lower LoA (", loa.alpha/2*100, " %)", sep = "")
    names(res)[3] <- paste("upper LoA (", (1-loa.alpha/2)*100, " %)", sep = "")
    print(res)
  }
  
  gg
}
