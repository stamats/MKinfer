pairwise.ext.t.test <- function(x, g, method = "t.test", p.adjust.method = "holm", 
                                paired = FALSE, ...){
  DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(g)))
  if (!is.na(pmatch(method, "t.test")))
    method <- "t.test"
  METHODS <- c("t.test", "boot.t.test", "perm.t.test", "hsu.t.test")
  method <- pmatch(method, METHODS)
  
  if (is.na(method))
    stop("invalid method")
  if (method == -1)
    stop("ambiguous method")
  
  if(method == 1){ # t.test
    tfun <- function(x, y, paired, ...){ 
      t.test(x, y, paired = paired, ...)
    }
  }
  if(method == 2){ # boot.t.test
    tfun <- function(x, y, paired, ...){ 
      boot.t.test(x, y, paired = paired, ...)
    }
  }
  if(method == 3){ # perm.t.test
    tfun <- function(x, y, paired, ...){ 
      perm.t.test(x, y, paired = paired, ...)
    }
  }
  if(method == 4){ # hsu.t.test
    tfun <- function(x, y, paired, ...){ 
      hsu.t.test(x, y, paired = paired, ...)
    }
  }
  ttests <- pairwise.fun(x = x, g = g, fun = tfun, paired = paired, ...)
  
  METH <- unlist(sapply(ttests, "[", "method"))[1]
  names(METH) <- NULL
  null.value <- unlist(sapply(ttests, "[", "null.value"))[1]
  if(!paired)
    names(null.value) <- "difference in means"
  else
    names(null.value) <- "mean of differences"
  alternative <- unlist(sapply(ttests, "[", "alternative"))[1]
  names(alternative) <- NULL
  statistic <- unlist(sapply(ttests, "[", "statistic"))
  names(statistic) <- NULL
  parameter <- unlist(sapply(ttests, "[", "parameter"))
  names(parameter) <- NULL
  if(method == 1 || method == 4){
    if(!paired){
      estimate <- (sapply(sapply(ttests, "[", "estimate"), "[", 1) - 
                     sapply(sapply(ttests, "[", "estimate"), "[", 2))
    } else {
      estimate <- unlist(sapply(ttests, "[", "estimate"))
    } 
  }
  if(method == 2)
    estimate <- unlist(sapply(ttests, "[", "boot.estimate"))
  if(method == 3)
    estimate <- unlist(sapply(ttests, "[", "perm.estimate"))
  names(estimate) <- NULL
  if(method == 1 || method == 4)
    stderr <- unlist(sapply(ttests, "[", "stderr"))
  if(method == 2)
    stderr <- unlist(sapply(ttests, "[", "boot.stderr"))
  if(method == 3)
    stderr <- unlist(sapply(ttests, "[", "perm.stderr"))
  names(stderr) <- NULL
  if(method == 1 || method == 4){
    CI.low <- sapply(sapply(ttests, "[", "conf.int"), "[", 1)
    names(CI.low) <- NULL
    CI.upp <- sapply(sapply(ttests, "[", "conf.int"), "[", 2)
    names(CI.upp) <- NULL
  }
  if(method == 2){
    CI.low <- sapply(sapply(ttests, "[", "boot.conf.int"), "[", 1)
    names(CI.low) <- NULL
    CI.upp <- sapply(sapply(ttests, "[", "boot.conf.int"), "[", 2)
    names(CI.upp) <- NULL
  }
  if(method == 3){
    CI.low <- sapply(sapply(ttests, "[", "perm.conf.int"), "[", 1)
    names(CI.low) <- NULL
    CI.upp <- sapply(sapply(ttests, "[", "perm.conf.int"), "[", 2)
    names(CI.upp) <- NULL
  }
  conf.level <- attr(ttests[[1]]$conf.int, which="conf.level")
  names(conf.level) <- NULL
  if(method == 1 || method == 4){
    p.value <- unlist(sapply(ttests, "[", "p.value"))
  }
  if(method == 2){
    p.value <- unlist(sapply(ttests, "[", "boot.p.value"))
  }
  if(method == 3){
    p.value <- unlist(sapply(ttests, "[", "perm.p.value"))
  }
  names(p.value) <- NULL
  
  res <- list(data.name = DNAME,
              method = paste0("Pairwise ", METH, "s"),
              null.value = null.value,
              alternative = alternative,
              conf.level = conf.level
              p.adjust.method = p.adjust.method)
  if(method == 1 || method == 4)
    res$results <- data.frame(groups = names(ttests),
                              t = statistic,
                              df = parameter,
                              estimate = estimate,
                              stderr = stderr,
                              CI.low = CI.low,
                              CI.upp = CI.upp,
                              p.value = p.value)
  else
    res$results <- data.frame(groups = names(ttests),
                              estimate = estimate,
                              stderr = stderr,
                              CI.low = CI.low,
                              CI.upp = CI.upp,
                              p.value = p.value)
  
  if(method == 1 || method == 4){
    if(!paired){
      names(res$results)[4] <- "diff. in means"
    } else {
      names(res$results)[4] <- "mean of diffs"
    }
  }else{
    if(!paired){
      names(res$results)[2] <- "diff. in means"
    } else {
      names(res$results)[2] <- "mean of diffs"
    }
  }
  res$results$adj.p.value <- p.adjust(res$results$p.value, 
                                      method = p.adjust.method)
  if(method == 1 || method == 4){
    names(res$results)[5] <- "SE"
    names(res$results)[6] <- paste0(100*(1-conf.level)/2, "%")
    names(res$results)[7] <- paste0(100*(1-(1-conf.level)/2), "%")
    names(res$results)[8] <- "p-value"
    names(res$results)[9] <- "adj. p-value"
  }else{
    names(res$results)[3] <- "SE"
    names(res$results)[4] <- paste0(100*(1-conf.level)/2, "%")
    names(res$results)[5] <- paste0(100*(1-(1-conf.level)/2), "%")
    names(res$results)[6] <- "p-value"
    names(res$results)[7] <- "adj. p-value"
  }
  class(res) <- "pw.htest"
  res
  
}