mpe.t.test <- function(X, Y, conf.level = 0.975){
  if(missing(X))
    stop("'X' must be given")
  if(missing(Y))
    stop("'Y' must be given")
  if(ncol(X) !=  ncol(Y))
    stop("'X' and 'Y' must have the same number of columns")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  ncols <- ncol(X)
  pvals <- numeric(ncols)
  statistic <- numeric(ncols)
  parameter <- numeric(ncols)
  conf <- matrix(NA, ncol = 2, nrow = ncols)
  estimate <- matrix(NA, ncol = 2, nrow = ncols)
  for(i in seq_len(ncols)){
    temp <- t.test(X[,i], Y[,i], var.equal = TRUE,
                   conf.level = conf.level)
    pvals[i] <- temp$p.value
    statistic[i] <- temp$statistic
    parameter[i] <- temp$parameter
    conf[i,] <- c(temp$conf.int[1], Inf)
    estimate[i,] <- temp$estimate
  }
  colnames(estimate) <- c("X", "Y")
  rownames(estimate) <- paste("EDP", 1:ncols, sep = ".")
  colnames(conf) <- c(paste(1-conf.level), 1)
  rownames(conf) <- paste("EP", 1:ncols, sep = ".")
  alternative <- "true difference in means is larger than 0 for all endpoints"
  attr(conf, "conf.level") <- conf.level
  method <- "Intersection-union t-test"
  pval <- max(pvals)
  stat <- min(statistic)
  names(stat) <- "t"
  para <- parameter[1]
  names(para) <- "df"
  rval <- list(method = method, statistic = stat, parameter = para, p.value = pval,
               conf.int = conf, estimate = estimate, alternative= alternative, data.name = dname)
  class(rval) <- "mpe.test"
  rval
}




