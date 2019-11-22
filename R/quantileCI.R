## Confidence Intervals for quantiles
quantileCI <- function(x, prob = 0.5, conf.level = 0.95, method = "exact",
                       R = 9999, bootci.type = c("norm", "basic", "perc", "bca"), 
                       minLength = FALSE, na.rm = FALSE,
                       alternative = c("two.sided", "less", "greater")){
    if(!is.na(pmatch(method, "exact")))
        method <- "exact"

    METHODS <- c("exact", "asymptotic", "boot")
    method <- pmatch(method, METHODS)
    alternative <- match.arg(alternative)
    
    if(is.na(method))
        stop("invalid method")

    if(method == -1)
        stop("ambiguous method")

    stopifnot(is.numeric(x), is.numeric(prob), is.numeric(conf.level))
    if(length(x) <= 1)
        stop("'x' has to be of at least length 2")
    if(length(prob) != 1)
      stop("'prob' has to be of length 1 (quantile)")
    if(prob <= 0 | prob >= 1)
      stop("'prob' has to be in (0, 1)")
    if(length(conf.level) != 1)
        stop("'conf.level' has to be of length 1 (confidence level)")
    if(conf.level < 0.5 | conf.level > 1)
        stop("'conf.level' has to be in [0.5, 1]")

    alpha <- 1 - conf.level
    
    if(alternative != "two.sided") alpha <- 2*alpha
    
    z <- qnorm(1-alpha/2)
    if(na.rm) x <- x[!is.na(x)]
    n <- length(x)
    est <- quantile(x, prob = prob)
    xs <- sort(x)

    if(method == 1){ # exact
        CI.mat <- matrix(NA, ncol = 2, nrow = n-1)
        pcov.vec <- numeric(n-1)
        for(i in 1:(n-1)){
          for(j in (i+1):n){
            pcov <- pbinom(j-1, size = n, prob = prob)-pbinom(i-1, size = n, prob = prob)
            if(pcov > conf.level){
              pcov.vec[i] <- pcov
              CI.mat[i,] <- c(xs[i], xs[j])
              break
            }
          }
        }
        if(all(pcov.vec == 0)){
          CI.lower <- ifelse(alternative == "less", -Inf, xs[1])
          CI.upper <- ifelse(alternative == "greater", Inf, xs[n])
          CI <- matrix(c(CI.lower, CI.upper), nrow = 1)
          attr(CI, "conf.level") <- 1
          rownames(CI) <- rep(paste(100*prob, "% quantile"), nrow(CI))
          colnames(CI) <- c("lower", "upper")
        }else{
          CI.mat <- CI.mat[pcov.vec > 0,,drop = FALSE]
          pcov.vec <- pcov.vec[pcov.vec > 0]
          pcov.min <- min(pcov.vec)
          CI <- CI.mat[pcov.vec == pcov.min,,drop = FALSE]
          if(minLength){
            CI <- CI[which.min(diff(t(CI))),,drop = FALSE]
          }
          if(alternative == "less") CI.mat[,1] <- -Inf
          if(alternative == "greater") CI.mat[,2] <- Inf
          attr(CI, "conf.level") <- pcov.min
          rownames(CI) <- rep(paste(100*prob, "% quantile"), nrow(CI))
          colnames(CI) <- c("lower", "upper")
        }
        if(minLength){
          meth <- paste("minimum length", METHODS[method], "confidence interval")
        }else{
          meth <- paste(METHODS[method], "confidence interval")
        }
    }
    if(method == 2){ # approx
        prob.sd <- sqrt(n*prob*(1-prob))
        k.lo <- max(1, floor(n*prob - z*prob.sd))
        k.up <- min(n, ceiling(n*prob + z*prob.sd))
        CI.lower <- ifelse(alternative == "less", -Inf, xs[k.lo])
        CI.upper <- ifelse(alternative == "greater", Inf, xs[k.up])
        CI <- matrix(c(CI.lower, CI.upper), nrow = 1)
        attr(CI, "conf.level") <- conf.level
        rownames(CI) <- rep(paste(100*prob, "% quantile"), nrow(CI))
        if(alternative == "two.sided")
          colnames(CI) <- c(paste(alpha/2*100, "%"), paste((1-alpha/2)*100, "%"))
        if(alternative == "less")
          colnames(CI) <- c("0 %", paste((1-alpha/2)*100, "%"))
        if(alternative == "greater")
          colnames(CI) <- c(paste(alpha/2*100, "%"), "100 %")
        meth <- "asymptotic confidence interval"
    }
    if(method == 3){ # boot
      boot.quant <- function(x, i){ 
        quantile(x[i], probs = prob) 
      } 
      boot.out <- boot(x, statistic = boot.quant, R = R)
      CI <- try(boot.ci(boot.out, type = bootci.type, conf = 1 - alpha),
                silent = TRUE)
      if(inherits(CI, "try-error"))
        stop("Function 'boot.ci' returned an error. Please try a different 'bootci.type'.")
      if(alternative == "less"){
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
      if(alternative == "greater"){
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
      meth <- "bootstrap confidence interval"
    }

    names(est) <- paste(100*prob, "% quantile")
    return(structure(list("estimate" = est, "conf.int" = CI,
                          "method" = meth),
                     class = "confint"))
}

medianCI <- function(x, conf.level = 0.95, method = "exact", 
                     R = 9999, bootci.type = c("norm", "basic", "perc", "bca"), 
                     minLength = FALSE, na.rm = FALSE,
                     alternative = c("two.sided", "less", "greater")){
    res <- quantileCI(x, prob = 0.5, conf.level = conf.level, method = method,
                      R = R, bootci.type = bootci.type, minLength = minLength,
                      na.rm = na.rm, alternative = alternative)
    if(method != "boot"){
      rownames(res$conf.int) <- rep("median", nrow(res$conf.int))
    }
    names(res$estimate) <- "median"
    res
}

madCI <- function(x, conf.level = 0.95, method = "exact", minLength = FALSE,
                  R = 9999, bootci.type = c("norm", "basic", "perc", "bca"),
                  na.rm = FALSE, constant = 1.4826,
                  alternative = c("two.sided", "less", "greater")){
  M <- median(x, na.rm = na.rm)
  res <- medianCI(constant*abs(x-M), conf.level = conf.level,
                  method = method, R = R, bootci.type = bootci.type, 
                  minLength = minLength, na.rm = na.rm, 
                  alternative = alternative)
  if(method != "boot"){
    rownames(res$conf.int) <- rep("MAD", nrow(res$conf.int))
  }
  names(res$estimate) <- "MAD"
  res
}
