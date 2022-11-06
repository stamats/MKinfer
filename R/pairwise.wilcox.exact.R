pairwise.wilcox.exact <- function(x, g, p.adjust.method = "holm", paired = FALSE, ...){
  DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(g)))
  wfun <- function(x, y, paired, ...){ 
    wilcox.exact(x,y, conf.int = TRUE, paired = paired, ...)
  }
  wtests <- pairwise.fun(x = x, g = g, fun = wfun, paired = paired, ...)

  method <- unlist(sapply(wtests, "[", 6))[1]
  names(method) <- NULL
  null.value <- unlist(sapply(wtests, "[", 4))[1]
  names(null.value) <- "location shift"
  alternative <- unlist(sapply(wtests, "[", 5))[1]
  names(alternative) <- NULL
  statistic <- unlist(sapply(wtests, "[", 1))
  names(statistic) <- NULL
  estimate <- unlist(sapply(wtests, "[", 9))
  estimate.name <- names(estimate)[1]
  names(estimate) <- NULL
  CI.low <- sapply(sapply(wtests, "[", 8), "[", 1)
  names(CI.low) <- NULL
  CI.upp <- sapply(sapply(wtests, "[", 8), "[", 2)
  names(CI.upp) <- NULL
  conf.level <- attr(wtests[[1]]$conf.int, which="conf.level")
  names(conf.level) <- NULL
  p.value <- unlist(sapply(wtests, "[", 3))
  names(p.value) <- NULL
  
  res <- list(data.name = DNAME,
              method = paste0("Pairwise ", method, "s"),
              null.value = null.value,
              alternative = alternative,
              conf.level = conf.level)
  
  res$results <- data.frame(groups = names(wtests),
                            W = statistic,
                            estimate = estimate,
                            CI.low = CI.low,
                            CI.upp = CI.upp,
                            p.value = p.value)
  
  if(!paired)
    names(res$results)[3] <- "diff. in location"
  else
    names(res$resulte)[3] <- "(pseudo)median"
  
  res$results$adj.p.value <- p.adjust(res$results$p.value, 
                                      method = p.adjust.method)
  names(res$results)[4] <- paste0(100*(1-conf.level)/2, "%")
  names(res$results)[5] <- paste0(100*(1-(1-conf.level)/2), "%")
  names(res$results)[6] <- "p-value"
  names(res$results)[7] <- "adj. p-value"
  class(res) <- "pw.htest"
  res
}
print.pw.htest <- function(x, digits = getOption("digits"), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative, two.sided = "not equal to", 
                           less = "less than", greater = "greater than")
        cat("true ", names(x$null.value), " is ", alt.char, 
            " ", x$null.value, "\n", sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, digits = digits, ...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  cat("\n")
  nc <- ncol(x$results)
  x$results[,(nc-4):nc] <- signif(x$results[,(nc-4):nc], digits = digits)
  print(x$results)
  cat("\n")
  invisible(x)
}
