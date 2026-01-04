## confidence interval for difference of means
normDiffCI <- function(x, y, conf.level = 0.95, paired = FALSE,
                       method = "welch",  boot = FALSE, R = 9999, 
                       bootci.type = "all", na.rm = TRUE,
                       alternative = c("two.sided", "less", "greater"), ...){
  if(!is.na(pmatch(method, "welch"))) method <- "welch"

  METHODS <- c("classical", "welch", "hsu", "xiao")
  method <- pmatch(method, METHODS)
  
  if (is.na(method))
    stop("invalid method")
  if (method == -1)
    stop("ambiguous method")

  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(conf.level) != 1)
    stop("'conf.level' has to be of length 1 (confidence level)")
  if(conf.level < 0.5 | conf.level > 1)
    stop("'conf.level' has to be in [0.5, 1]")

  alternative <- match.arg(alternative)
  
  Infos <- NULL
  alpha <- 1 - conf.level
  
  if(alternative != "two.sided") alpha <- 2*alpha
  
  if(paired){
    CIres <- normCI(x = x-y, conf.level = conf.level, boot = boot, R = R, 
                    bootci.type = bootci.type, na.rm = na.rm, 
                    alternative = alternative, ...)
    d <- CIres$estimate
    names(d) <- c("mean of differences", "sd of differences")
    if(boot){
      CI <- CIres$conf.int$mean
      Infos <- NULL
      method <- "Bootstrap confidence interval (paired)"
    }else{
      CI.lower <- CIres$conf.int[1,1]
      CI.upper <- CIres$conf.int[1,2]
      Infos <- CIres$Infos
      names(Infos) <- "SE of mean of differences"
      method <- "Confidence interval (paired)"
    }
  }else{
    mx <- mean(x, na.rm = na.rm)
    my <- mean(y, na.rm = na.rm)
    d <- mx - my
    if(na.rm) nx <- length(x[!is.na(x)]) else nx <- length(x)
    if(na.rm) ny <- length(y[!is.na(y)]) else ny <- length(y)
    vx <- var(x, na.rm = na.rm)
    vy <- var(y, na.rm = na.rm)
    if(method == 1){
      df <- nx + ny - 2
      s <- sqrt(((nx-1)*vx + (ny-1)*vy)/df)
      se <- s*sqrt(1/nx + 1/ny)
      if(boot){
        if(na.rm){
          DATA <- data.frame(vals = c(x[!is.na(x)], y[!is.na(y)]), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }else{
          DATA <- data.frame(vals = c(x, y), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }
        boot.diff.class <- function(x, i){
          xi <- x[i,]
          AMs <- tapply(xi$vals, xi$group, mean)
          d <- -diff(AMs)
          Vs <- tapply(xi$vals, xi$group, mean)
          Ns <- table(xi$group)
          df <- Ns[1] + Ns[2] - 2
          v <- ((Ns[1]-1)*Vs[1] + (Ns[2]-1)*Vs[2])/df
          VAR <- v*(1/Ns[1] + 1/Ns[2])
          c(d, VAR)
        } 
        boot.out <- boot(DATA, statistic = boot.diff.class, 
                         strata=DATA$group, R = R, ...)
        CI <- suppressWarnings(boot.ci(boot.out, type = bootci.type, conf = 1 - alpha))
        method <- "Bootstrap confidence interval (equal variances, unpaired)"
      }else{
        method <- "Classical confidence interval (unpaired)"
      }
    }
    if(method == 2){
      sex <- sqrt(vx/nx)
      sey <- sqrt(vy/ny)
      s <- sqrt(vx + vy)
      se <- sqrt(vx/nx + vy/ny)
      df <- se^4/(sex^4/(nx - 1) + sey^4/(ny - 1))
      if(boot){
        if(na.rm){
          DATA <- data.frame(vals = c(x[!is.na(x)], y[!is.na(y)]), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }else{
          DATA <- data.frame(vals = c(x, y), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }
        boot.diff.welch <- function(x, i){
          xi <- x[i,]
          AMs <- tapply(xi$vals, xi$group, mean)
          d <- -diff(AMs)
          Vs <- tapply(xi$vals, xi$group, mean)
          Ns <- table(xi$group)
          VAR <- Vs[1]/Ns[1] + Vs[2]/Ns[2]
          c(d, VAR)
        } 
        boot.out <- boot(DATA, statistic = boot.diff.welch, 
                         strata=DATA$group, R = R, ...)
        CI <- suppressWarnings(boot.ci(boot.out, type = bootci.type, conf = 1 - alpha))
        method <- "Bootstrap confidence interval (unequal variances, unpaired)"
      }else{
        method <- "Welch confidence interval (unpaired)"
      }
    }
    if(method == 3){
      if(nx < 6 || ny < 6)
        warning("For method of Hsu the sample size per group should be > 5.")
      s <- sqrt(vx + vy)
      se <- sqrt(vx/nx + vy/ny)
      df <- min(nx, ny) - 1
      if(boot){
        if(na.rm){
          DATA <- data.frame(vals = c(x[!is.na(x)], y[!is.na(y)]), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }else{
          DATA <- data.frame(vals = c(x, y), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }
        ## same function as in case of welch
        boot.diff.hsu <- function(x, i){
          xi <- x[i,]
          AMs <- tapply(xi$vals, xi$group, mean)
          d <- -diff(AMs)
          Vs <- tapply(xi$vals, xi$group, mean)
          Ns <- table(xi$group)
          VAR <- Vs[1]/Ns[1] + Vs[2]/Ns[2]
          c(d, VAR)
        } 
        boot.out <- boot(DATA, statistic = boot.diff.hsu, 
                         strata=DATA$group, R = R, ...)
        CI <- suppressWarnings(boot.ci(boot.out, type = bootci.type, conf = 1 - alpha))
        method <- "Bootstrap confidence interval (unequal variances, unpaired)"
      }else{
        method <- "Hsu confidence interval (unpaired)"
      }
    }
    if(method == 4){
      s <- sqrt(vx + vy)
      se <- sqrt(vx/nx + vy/ny)
      if(nx*(nx-1)/vx <= ny*(ny-1)/vy){
        df <- (ny-1)*(1 + ny/nx*vx/vy)
      }else{
        df <- (nx-1)*(1 + nx/ny*vy/vx) 
      }
      if(boot){
        if(na.rm){
          DATA <- data.frame(vals = c(x[!is.na(x)], y[!is.na(y)]), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }else{
          DATA <- data.frame(vals = c(x, y), 
                             group = factor(rep(c(1,2), c(nx, ny))))
        }
        ## same function as in case of welch
        boot.diff.xiao <- function(x, i){
          xi <- x[i,]
          AMs <- tapply(xi$vals, xi$group, mean)
          d <- -diff(AMs)
          Vs <- tapply(xi$vals, xi$group, mean)
          Ns <- table(xi$group)
          VAR <- Vs[1]/Ns[1] + Vs[2]/Ns[2]
          c(d, VAR)
        } 
        boot.out <- boot(DATA, statistic = boot.diff.xiao, 
                         strata=DATA$group, R = R, ...)
        CI <- suppressWarnings(boot.ci(boot.out, type = bootci.type, conf = 1 - alpha))
        method <- "Bootstrap confidence interval (unequal variances, unpaired)"
      }else{
        method <- "Xiao confidence interval (unpaired)"
      }
    }
    if(!boot){
      if(method == 4){
        t.alpha <- qgt(1-alpha/2, n1 = nx, n2 = ny, v1tov2 = vx/vy)
      }else{
        t.alpha <- qt(1-alpha/2, df = df)
      }
      ## confidence bounds
      CI.lower <- ifelse(alternative == "less", -Inf, d - t.alpha*se)
      CI.upper <- ifelse(alternative == "greater", Inf, d + t.alpha*se)
    }else{
      if(alternative == "less" && !paired){
        if("normal" %in% names(CI)){ 
          CI$normal[1,1] <- conf.level
          CI$normal[1,2] <- -Inf
        }
        if("basic" %in% names(CI)){ 
          CI$basic[1,1] <- conf.level
          CI$basic[1,4] <- -Inf
        }
        if("student" %in% names(CI)){ 
          CI$student[1,1] <- conf.level
          CI$student[1,4] <- -Inf
        }
        if("percent" %in% names(CI)){ 
          CI$percent[1,1] <- conf.level
          CI$percent[1,4] <- -Inf
        }
        if("bca" %in% names(CI)){ 
          CI$bca[1,1] <- conf.level
          CI$bca[1,4] <- -Inf
        }
      }
      if(alternative == "greater" && !paired){
        if("normal" %in% names(CI)){ 
          CI$normal[1,1] <- conf.level
          CI$normal[1,3] <- Inf
        }
        if("basic" %in% names(CI)){ 
          CI$basic[1,1] <- conf.level
          CI$basic[1,5] <- Inf
        }
        if("student" %in% names(CI)){ 
          CI$student[1,1] <- conf.level
          CI$student[1,5] <- Inf
        }
        if("percent" %in% names(CI)){ 
          CI$percent[1,1] <- conf.level
          CI$percent[1,5] <- Inf
        }
        if("bca" %in% names(CI)){ 
          CI$bca[1,1] <- conf.level
          CI$bca[1,5] <- Inf
        }
      }
    }
    names(d) <- "difference in means"
    Infos <- list(c(se, d/s), c(mx, sqrt(vx), my, sqrt(vy)))
    names(Infos[[1]]) <- c("SE of difference in means", "Cohen's d (SMD)")
    names(Infos[[2]]) <- c("mean of x", "SD of x", "mean of y", "SD of y")
  }

  if(!boot){
    CI <- matrix(c(CI.lower, CI.upper), nrow = 1)
    if(paired){
      rownames(CI) <- "mean of differences"
    }else{
      rownames(CI) <- "difference in means"
    }
    if(alternative == "two.sided")
      colnames(CI) <- c(paste(alpha/2*100, "%"), paste((1-alpha/2)*100, "%"))
    if(alternative == "less")
      colnames(CI) <- c("0 %", paste((1-alpha/2)*100, "%"))
    if(alternative == "greater")
      colnames(CI) <- c(paste(alpha/2*100, "%"), "100 %")
    attr(CI, "conf.level") <- conf.level
  }

  return(structure(list("estimate" = d, "conf.int" = CI, "Infos" = Infos,
                        method = method),
                   class = "confint"))
}
meanDiffCI <- function(x, y, conf.level = 0.95, paired = FALSE,
                       method = "welch",  boot = FALSE, R = 9999, 
                       bootci.type = "all", na.rm = TRUE,
                       alternative = c("two.sided", "less", "greater"), ...){
  normDiffCI(x = x, y = y, conf.level = conf.level, paired = paired, 
             method = method, boot = boot, R = R, 
             bootci.type = bootci.type, na.rm = na.rm, 
             alternative = alternative, ...)
}
