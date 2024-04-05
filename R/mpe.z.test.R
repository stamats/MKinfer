mpe.z.test <- function(X, Y, Sigma, conf.level = 0.975){
  if(missing(X))
    stop("'X' must be given")
  if(missing(Y))
    stop("'Y' must be given")
  if(ncol(X) !=  ncol(Y))
    stop("'X' and 'Y' must have the same number of columns")
  if(missing(Sigma))
    stop("'Sigma' must be given")
  if(nrow(Sigma) != ncol(Sigma))
    stop("covariance matrix 'Sigma' must be quadratic")
  if(nrow(Sigma) != ncol(X))
    stop("covariance matrix has wrong dimension")
  if(max(abs(Sigma - t(Sigma))) > 1e-10)
    stop("matrix 'Sigma' must be symmetric")
  if(!all(eigen(Sigma)$values > 0))
    stop("matrix 'Sigma' must be positive definite")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  diff <- colMeans(X) - colMeans(Y)
  names(diff) <- paste("EP", 1:ncol(X), sep = ".")
  n <- nrow(X)
  m <- nrow(Y)
  SD <- sqrt(diag(Sigma))
  Z <- diff/(SD*sqrt(1/n+1/m))
  statistic <- min(Z)
  names(statistic) <- "z"
  pval <- max(pnorm(Z, lower.tail = FALSE))
  alpha <- 1 - conf.level
  cint <- rep(qnorm(1 - alpha), ncol(X))
  cint <- cbind(Z - cint, Inf)*SD*sqrt(1/n+1/m)

  colnames(cint) <- c(paste(1-conf.level), 1)
  rownames(cint) <- paste("EP", 1:ncol(X), sep = ".")

  alternative <- "true difference in means is larger than 0 for all endpoints"
  attr(cint, "conf.level") <- conf.level
  method <- "Intersection-union z-test"
  rval <- list(method = method, statistic = statistic, p.value = pval,
               conf.int = cint, estimate = diff,
               alternative = alternative, data.name = dname)
  class(rval) <- "mpe.test"
  rval
}

print.mpe.test <- function (x, digits = getOption("digits"), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=", format(signif(x$statistic,
                                                               max(1L, digits - 2L)))))
  if (!is.null(x$parameter))
    out <- c(out, paste(names(x$parameter), "=", format(signif(x$parameter,
                                                               max(1L, digits - 2L)))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits -
                                                3L))
    out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) ==
                                       "<") fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ", x$alternative, "\n")
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence intervals:\n", sep = "")
    attr(x$conf.int, "conf.level") <- NULL
    print(x$conf.int)
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}

