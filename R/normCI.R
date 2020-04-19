## Confidence Intervals for normal mean and standard deviation
normCI <- function(x, mean = NULL, sd = NULL, conf.level = 0.95, 
                   boot = FALSE, R = 9999, bootci.type = "all", na.rm = TRUE,
                   alternative = c("two.sided", "less", "greater")){
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.null(mean))
    if(length(mean) != 1) stop("'mean' has to be of length 1")
  if(!is.null(sd))
    if(length(sd) != 1) stop("'sd' has to be of length 1")
  if(!is.null(sd))
    if(sd <= 0) stop("'sd' has to be positive")
  if(length(conf.level) != 1)
    stop("'conf.level' has to be of length 1 (confidence level)")
  if(conf.level < 0.5 | conf.level > 1)
    stop("'conf.level' has to be in [0.5, 1]")

  alternative <- match.arg(alternative)
  
  alpha <- 1 - conf.level
  if(alternative != "two.sided") alpha <- 2*alpha
  
  Infos <- NULL
  if(is.null(mean)){
    m <- mean(x, na.rm = na.rm)
  }else{
    m <- mean
  }
  if(is.null(sd)){
    s <- sd(x, na.rm = na.rm)
  }else{
    s <- sd
  }
  if(is.null(mean) & is.null(sd)){
    est <- c(m, s)
    names(est) <- c("mean", "sd")
  }else{
    if(is.null(mean)){
      est <- m
      names(est) <- "mean"
    }else{
      est <- s
      names(est) <- "sd"
    }
  }
  if(boot){
    if(na.rm) x <- x[!is.na(x)]
    boot.mean <- function(x, i){
      AM <- mean(x[i]) 
      n <- length(i)
      VAR <- (n-1)*var(x[i])/n^2
      c(AM, VAR)
    } 
    boot.sd <- function(x, i){
      SD <- sd(x[i]) 
      n <- length(i)
      VAR <- SD^2*((n-1)/2*gamma((n-1)/2)/gamma(n/2) - 1)
      c(SD, VAR)
    } 
    if(is.null(mean)){
      boot.out <- boot(x, statistic = boot.mean, R = R)
      CI.AM <- try(boot.ci(boot.out, type = bootci.type, conf = 1-alpha),
                silent = TRUE)
      if(inherits(CI.AM, "try-error"))
        stop("Function 'boot.ci' returned an error. Please try a different 'bootci.type'.")
      if(alternative == "less"){
        if("normal" %in% names(CI.AM)){ 
          CI.AM$normal[1,1] <- conf.level
          CI.AM$normal[1,2] <- -Inf
        }
        if("basic" %in% names(CI.AM)){ 
          CI.AM$basic[1,1] <- conf.level
          CI.AM$basic[1,4] <- -Inf
        }
        if("student" %in% names(CI.AM)){ 
          CI.AM$student[1,1] <- conf.level
          CI.AM$student[1,4] <- -Inf
        }
        if("percent" %in% names(CI.AM)){ 
          CI.AM$percent[1,1] <- conf.level
          CI.AM$percent[1,4] <- -Inf
        }
        if("bca" %in% names(CI.AM)){ 
          CI.AM$bca[1,1] <- conf.level
          CI.AM$bca[1,4] <- -Inf
        }
      }
      if(alternative == "greater"){
        if("normal" %in% names(CI.AM)){ 
          CI.AM$normal[1,1] <- conf.level
          CI.AM$normal[1,3] <- Inf
        }
        if("basic" %in% names(CI.AM)){ 
          CI.AM$basic[1,1] <- conf.level
          CI.AM$basic[1,5] <- Inf
        }
        if("student" %in% names(CI.AM)){ 
          CI.AM$student[1,1] <- conf.level
          CI.AM$student[1,5] <- Inf
        }
        if("percent" %in% names(CI.AM)){ 
          CI.AM$percent[1,1] <- conf.level
          CI.AM$percent[1,5] <- Inf
        }
        if("bca" %in% names(CI.AM)){ 
          CI.AM$bca[1,1] <- conf.level
          CI.AM$bca[1,5] <- Inf
        }
      }
    }    
    if(is.null(sd)){
      boot.out <- boot(x, statistic = boot.sd, R = R)
      CI.SD <- try(boot.ci(boot.out, type = bootci.type, conf = 1-alpha),
                silent = TRUE)
      if(inherits(CI.SD, "try-error"))
        stop("Function 'boot.ci' returned an error. Please try a different 'bootci.type'.")
      if(alternative == "less"){
        if("normal" %in% names(CI.SD)){ 
          CI.SD$normal[1,1] <- conf.level
          CI.SD$normal[1,2] <- 0
        }
        if("basic" %in% names(CI.SD)){ 
          CI.SD$basic[1,1] <- conf.level
          CI.SD$basic[1,4] <- 0
        }
        if("student" %in% names(CI.SD)){ 
          CI.SD$student[1,1] <- conf.level
          CI.SD$student[1,4] <- 0
        }
        if("percent" %in% names(CI.SD)){ 
          CI.SD$percent[1,1] <- conf.level
          CI.SD$percent[1,4] <- 0
        }
        if("bca" %in% names(CI.SD)){ 
          CI.SD$bca[1,1] <- conf.level
          CI.SD$bca[1,4] <- 0
        }
      }
      if(alternative == "greater"){
        if("normal" %in% names(CI.SD)){ 
          CI.SD$normal[1,1] <- conf.level
          CI.SD$normal[1,3] <- Inf
        }
        if("basic" %in% names(CI.SD)){ 
          CI.SD$basic[1,1] <- conf.level
          CI.SD$basic[1,5] <- Inf
        }
        if("student" %in% names(CI.SD)){ 
          CI.SD$student[1,1] <- conf.level
          CI.SD$student[1,5] <- Inf
        }
        if("percent" %in% names(CI.SD)){ 
          CI.SD$percent[1,1] <- conf.level
          CI.SD$percent[1,5] <- Inf
        }
        if("bca" %in% names(CI.SD)){ 
          CI.SD$bca[1,1] <- conf.level
          CI.SD$bca[1,5] <- Inf
        }
      }
    }
    
    if(is.null(mean) && is.null(sd))
      CI <- list("mean" = CI.AM, "standard deviation" = CI.SD)
    if(is.null(mean) && !is.null(sd)) CI <- CI.AM
    if(!is.null(mean) && is.null(sd)) CI <- CI.SD
    
    METHOD <- "Bootstrap confidence interval(s)"
  }else{
    if(na.rm) n <- length(x[!is.na(x)]) else n <- length(x)
    if(is.null(sd)){
      k <- qt(1-alpha/2, df = n-1)
    }else{
      k <- qnorm(1-alpha/2)
    }
    if(is.null(mean)){
      sem <- s/sqrt(n)
      names(sem) <- "SE of mean"
      CI.lower.mean <- ifelse(alternative == "less", -Inf, m - k*sem)
      CI.upper.mean <- ifelse(alternative == "greater", Inf, m + k*sem)
      Infos <- sem
    }
    if(is.null(sd)){
      CI.lower.sd <- ifelse(alternative == "less", 0, sqrt(n-1)*s/sqrt(qchisq(1-alpha/2, df = n-1)))
      CI.upper.sd <- ifelse(alternative == "greater", Inf, sqrt(n-1)*s/sqrt(qchisq(alpha/2, df = n-1)))
    }
    if(is.null(mean) & is.null(sd)){
      CI <- rbind(c(CI.lower.mean, CI.upper.mean),
                  c(CI.lower.sd, CI.upper.sd))
      rownames(CI) <- c("mean", "sd")
    }else{
      if(is.null(mean)){
        CI <- matrix(c(CI.lower.mean, CI.upper.mean), nrow = 1)
        rownames(CI) <- "mean"
      }else{
        CI <- matrix(c(CI.lower.sd, CI.upper.sd), nrow = 1)
        rownames(CI) <- "sd"
      }
    }
    if(alternative == "two.sided")
      colnames(CI) <- c(paste(alpha/2*100, "%"), paste((1-alpha/2)*100, "%"))
    if(alternative == "less")
      colnames(CI) <- c("0 %", paste((1-alpha/2)*100, "%"))
    if(alternative == "greater")
      colnames(CI) <- c(paste(alpha/2*100, "%"), "100 %")
    attr(CI, "conf.level") <- conf.level
    METHOD <- "Exact confidence interval(s)"
  }

  return(structure(list("estimate" = est, "conf.int" = CI, "Infos" = Infos,
                        method = METHOD),
                   class = "confint"))
}
meanCI <- function(x, conf.level = 0.95, boot = FALSE, R = 9999, 
                   bootci.type = "all", na.rm = TRUE,
                   alternative = c("two.sided", "less", "greater")){
  res <- normCI(x = x, conf.level = conf.level, boot = boot, R = R,
                bootci.type = bootci.type, na.rm = na.rm, 
                alternative = alternative)
  if(boot){
    res$conf.int <- res$conf.int[[1]]
  }else{
    res$conf.int <- res$conf.int[1,,drop = FALSE]
  }
  attr(res$conf.int, "conf.level") <- conf.level
  res
}
sdCI <- function(x, conf.level = 0.95, boot = FALSE, R = 9999, 
                 bootci.type = "all", na.rm = TRUE,
                 alternative = c("two.sided", "less", "greater")){
  res <- normCI(x = x, conf.level = conf.level, boot = boot, R = R,
                bootci.type = bootci.type, na.rm = na.rm, 
                alternative = alternative)
  if(boot){
    res$conf.int <- res$conf.int[[2]]
  }else{
    res$conf.int <- res$conf.int[2,,drop = FALSE]
  }
  attr(res$conf.int, "conf.level") <- conf.level
  res
}

print.confint <- function(x, digits = getOption("digits"), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  out <- x$conf.int
  attr(out, "conf.level") <- NULL
  if(is.list(out) | inherits(out, "bootci")){
    print(out, digits = digits, ...)
  }else{
    if (nrow(x$conf.int) > 1) {
      cat(format(100 * attr(x$conf.int, "conf.level")),
          " percent confidence intervals:\n", sep = "")
    }else{
      cat(format(100 * attr(x$conf.int, "conf.level")),
          " percent confidence interval:\n", sep = "")
    }
    print(out, digits = digits, ...)
  }
  if (!is.null(x$estimate)) {
    cat("\n")
    if(length(x$estimate) == 1) cat("sample estimate:\n")
    if(length(x$estimate) > 1) cat("sample estimates:\n")
    print(x$estimate, digits = digits, ...)
  }
  if (!is.null(x$Infos)) {
    cat("\n")
    cat("additional information:\n")
    if(is.list(x$Infos) & length(x$Infos) > 1){
      for(i in seq_len(length(x$Infos))){
        if(i > 1) cat("\n")
        print(x$Infos[[i]], digits = digits, ...)
      }
    }else{
      print(x$Infos, digits = digits, ...)
    }
  }
  cat("\n")
  invisible(x)
}
