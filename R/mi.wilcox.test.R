mi.wilcox.test <- function(miData, ...){
  UseMethod("mi.wilcox.test")
}
mi.wilcox.test.default <- function(miData, x, y = NULL, 
                                   alternative = c("two.sided", "less", "greater"),
                                   mu = 0, paired = FALSE, exact = NULL, conf.int = TRUE,
                                   conf.level = 0.95, subset = NULL, ...){
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if(!is.list(miData))
    stop("'miData' must be a list of imputed datasets")
  if(is.null(y) & paired)
    stop("'y' is missing for paired test")

  nrImp <- length(miData)
  resTests <- vector(mode = "list", length = nrImp)
  for(i in 1:nrImp){
    if(is.null(subset))
      xi <- miData[[i]][,x]
    else
      xi <- miData[[i]][subset, x]
    if(is.null(y))
      yi <- NULL
    else{
      if(is.null(subset)){
        yi <- miData[[i]][,y]
      }else{
        yi <- miData[[i]][subset,y]
        if(is.factor(yi)) yi <- factor(yi)
      }
    }
    if(paired){
      xi <- xi-yi
      yi <- NULL
    }
    if(is.null(yi)){
      resTests[[i]] <- wilcox.exact(xi, alternative = alternative, mu = mu, 
                                    paired = FALSE, exact = exact, conf.int = conf.int, 
                                    conf.level = conf.level, ...)
    }else{
      resTests[[i]] <- wilcox.exact(xi ~ yi, alternative = alternative, mu = mu, 
                                    paired = FALSE, exact = exact, conf.int = conf.int, 
                                    conf.level = conf.level, ...)
    }
  }
  pvals <- sapply(resTests, function(x) x$p.value)
  ind.med <- which.min(abs(pvals - median(pvals)))
  res <- resTests[[ind.med]]
  res$method <- paste("Multiple Imputation", res$method)
  if(is.null(yi)){
    if(paired){
      dname <- paste("Variables ", x, " and ", y, "\n", 
                     "number of imputations: ", nrImp, sep = "")
    }else{
      dname <- paste("Variable ", x, "\n", 
                     "number of imputations: ", nrImp, sep = "")
    }
  }else{
    dname <- paste("Variable ", x, ": ", paste(paste("group", levels(yi)),
                                               collapse = " vs "), 
                   "\n", "number of imputations: ", nrImp, sep = "")
  }
  res$data.name <- dname
  res
}

mi.wilcox.test.amelia <- function(miData, x, y = NULL, 
                                  alternative = c("two.sided", "less", "greater"),
                                  mu = 0, paired = FALSE, exact = NULL, conf.int = TRUE,
                                  conf.level = 0.95, subset = NULL, ...){
  mi.wilcox.test(miData$imputations, x = x, y = y, alternative = alternative,
                 mu = mu, paired = paired, exact = exact, conf.int = conf.int, 
                 conf.level = conf.level, subset = subset, ...)
}

mi.wilcox.test.mids <- function(miData, x, y = NULL, 
                                alternative = c("two.sided", "less", "greater"),
                                mu = 0, paired = FALSE, exact = NULL, conf.int = TRUE,
                                conf.level = 0.95, subset = NULL, ...){
  mi.wilcox.test(mids2datlist(miData), x = x, y = y, alternative = alternative,
                 mu = mu, paired = paired, exact = exact, conf.int = conf.int, 
                 conf.level = conf.level, subset = subset, ...)
}
