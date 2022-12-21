rm.oneway.test <- function (x, g, id, method = "aov") {
  if (!is.na(pmatch(method, "aov")))
    method <- "aov"
  METHODS <- c("aov", "lme", "friedman", "quade")
  method <- pmatch(method, METHODS)
  
  stopifnot(is.numeric(x))
  if(length(x) != length(g) || length(x) != length(id)) 
    stop("'x', 'g' and 'id' must have the same length")
  DNAME <- paste(deparse1(substitute(x)), ",", deparse1(substitute(g)), "and", 
                 deparse1(substitute(id)))
  OK <- complete.cases(x, g, id)
  x <- x[OK]
  g <- g[OK]
  id <- id[OK]
  g <- factor(g)
  k <- nlevels(g)
  id <- factor(id)
  if (k < 2L) 
    stop("all observations are in the same group")
  n <- length(x)
  if (n < 2L) 
    stop("not enough observations")
  
  if(method == 1){ ## aov
    res <- summary(aov(x ~ g + Error(id/g)))
    STATISTIC <- res$`Error: id:g`[[1]][1,"F value"]
    names(STATISTIC) <- "F"
    PARAMETER <- res$`Error: id:g`[[1]][,"Df"]
    names(PARAMETER) <- c("num df", "denom df")
    PVAL <- res$`Error: id:g`[[1]][1,"Pr(>F)"]
    METHOD <- "Repeated measures 1-way ANOVA"
    RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME)
    class(RVAL) <- "htest"
  }
  if(method == 2){ ## lme
    res <- anova(lme(x ~ g, random = ~1 | id))
    STATISTIC <- res$`F-value`[2]
    names(STATISTIC) <- "F"
    PARAMETER <- c(res$numDF[2], res$denDF[2])
    names(PARAMETER) <- c("num df", "denom df")
    PVAL <- res$`p-value`[2]
    METHOD <- "Mixed-effects 1-way ANOVA"
    RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME)
    class(RVAL) <- "htest"
  }
  if(method == 3){ ## friedman
    RVAL <- friedman.test(y = x, groups = g, blocks = id)
  }
  if(method == 4){ ## quade
    RVAL <- quade.test(y = x, groups = g, blocks = id)
  }
  RVAL
}
