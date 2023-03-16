## Confidence Intervals for quantiles
quantileCI <- function(x, prob = 0.5, conf.level = 0.95, method = "exact",
                       R = 9999, bootci.type = c("norm", "basic", "perc", "bca"), 
                       na.rm = FALSE, alternative = c("two.sided", "less", "greater"), ...){
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
        if(prob^n + (1-prob)^n > alpha){
          fun.n <- function(n, prob, alpha){
            prob^n + (1-prob)^n - alpha
          }
          n.min <- ceiling(uniroot(f = fun.n, lower = 1, upper = 5e4, prob = prob, 
                                   alpha = alpha, extendInt = "yes")$root)
          warning("There is no exact confidence interval that achieves the requested confidence level of at least ", conf.level, 
                  "\n", "A sample size of at least ", n.min, " would be required!")
        }
        fun.cover <- function(l, u, n, prob){
          pbinom(u-1, size = n, prob = prob) - pbinom(l-1, size = n, prob = prob) 
        }
        res <- outer(1:n, 1:n, fun.cover, n = n, prob = prob)
        res[lower.tri(res, diag = TRUE)] <- NA
        res <- as.vector(res)
        names.res <- outer(1:n, 1:n, function(x, y) paste(x, y, sep = " | "))
        names(res) <- names.res
        
        if(any(res > conf.level, na.rm = TRUE)){
          res <- res[res > 1-alpha]
          ind.min <- which.min(res - (1-alpha))
        }else{
          ind.min <- which.min(abs(res-(1-alpha)))
        }
        conf.level.exact <- res[ind.min]
        names(conf.level.exact) <- NULL
        res.sel <- res[which(res == conf.level.exact)]
        if(length(res.sel) > 1){
          bounds.all <- do.call("rbind", strsplit(names(res.sel), " \\| "))
          CI.length <- numeric(nrow(bounds.all))
          for(i in 1:nrow(bounds.all)){
            ind <- as.numeric(bounds.all[i,])
            CI.length[i] <- diff(xs[ind])
          }
          ind.min <- which.min(CI.length)
          bounds <- as.numeric(bounds.all[ind.min,])
        }else{
          bounds <- as.numeric(strsplit(names(res.sel), " \\| ")[[1]])
        }
        lower <- bounds[1]
        upper <- bounds[2]

        # plev <- 0
        # upper <- ceiling((n+1)*prob)
        # lower <- floor((n+1)*prob)
        # while(plev < conf.level){
        #   upper <- upper + 1
        #   lower <- lower - 1
        #   plev <- pbinom(upper-1, size = n, prob = prob) - pbinom(lower-1, size = n, prob = prob) 
        # }
        # lower <- max(lower, 1)
        # upper <- min(upper, n)
        
        CI.lower <- ifelse(alternative == "less", -Inf, xs[lower])
        CI.upper <- ifelse(alternative == "greater", Inf, xs[upper])
        CI <- matrix(c(CI.lower, CI.upper), nrow = 1)
        if(alternative != "two.sided") conf.level.exact <- 1-(1-conf.level.exact)/2
        alpha.exact <- 1 - round(conf.level.exact, 3)
        attr(CI, "conf.level") <- round(conf.level.exact, 3)
        rownames(CI) <- paste(100*prob, "% quantile")
        if(alternative == "two.sided")
          colnames(CI) <- c(paste(alpha.exact/2*100, "%"), paste((1-alpha.exact/2)*100, "%"))
        if(alternative == "less")
          colnames(CI) <- c("0 %", paste((1-alpha.exact)*100, "%"))
        if(alternative == "greater")
          colnames(CI) <- c(paste(alpha.exact*100, "%"), "100 %")
        meth <- "exact confidence interval"
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
      boot.out <- boot(x, statistic = boot.quant, R = R, ...)
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

medianCI <- function(x, conf.level = 0.95, method = "exact", R = 9999, 
                     bootci.type = c("norm", "basic", "perc", "bca"), na.rm = FALSE, 
                     alternative = c("two.sided", "less", "greater"), ...){
    res <- quantileCI(x, prob = 0.5, conf.level = conf.level, method = method,
                      R = R, bootci.type = bootci.type, 
                      na.rm = na.rm, alternative = alternative, ...)
    if(method != "boot"){
      rownames(res$conf.int) <- rep("median", nrow(res$conf.int))
    }
    names(res$estimate) <- "median"
    res
}

madCI <- function(x, conf.level = 0.95, method = "exact", 
                  R = 9999, bootci.type = c("norm", "basic", "perc", "bca"),
                  na.rm = FALSE, constant = 1.4826,
                  alternative = c("two.sided", "less", "greater"), ...){
  alternative <- match.arg(alternative)
  if(method %in% c("exact", "asymptotic")){
    M <- median(x, na.rm = na.rm)
    res <- medianCI(constant*abs(x-M), conf.level = conf.level,
                    method = method, R = R, bootci.type = bootci.type, 
                    na.rm = na.rm, alternative = alternative, ...)
    if(method != "boot"){
      rownames(res$conf.int) <- rep("MAD", nrow(res$conf.int))
    }
    names(res$estimate) <- "MAD"
  }else{
    if(method == "boot"){ # boot
      if(na.rm) x <- x[!is.na(x)]
      alpha <- 1 - conf.level
      if(alternative != "two.sided") alpha <- 2*alpha
      
      est <- mad(x, constant = constant)
      boot.mad <- function(x, i, constant){ 
        mad(x[i], constant = constant) 
      } 
      boot.out <- boot(x, statistic = boot.mad, R = R, constant = constant, ...)
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
    names(est) <- "MAD"
    res <- structure(list("estimate" = est, "conf.int" = CI,
                          "method" = meth),
                     class = "confint")
  }
  res
}
