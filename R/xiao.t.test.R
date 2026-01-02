## generalized central t distribution as defined in Section 5.2 in Xiao (2018)
dgt <- function(x, n1, n2, sd1, sd2){
  T1 <- (n2-1)*(1 + n2/n1*sd1^2/sd2^2)
  T2 <- (n1-1)*n1*sd2^2/((n2-1)*n2*sd1^2)
  T3 <- Re(hypergeo((n1+n2-1)/2, (n1-1)/2, (n1+n2-2)/2, 1 - (T2+x^2/T1)/(1+x^2/T1)))
  1/sqrt(pi*T1)*gamma((n1+n2-1)/2)/gamma((n1+n2-2)/2)*
    (1+x^2/T1)^(-(n1+n2-1)/2)*T2^((n1-1)/2)*T3
}
pgt <- function(q, n1, n2, sd1, sd2, lower.tail = TRUE){
  ifelse(lower.tail, 
         integrate(f = dgt, lower = -Inf, upper = q, n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2)$value,
         integrate(f = dgt, lower = q, upper = Inf, n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2)$value)
}
qgt <- function(p, n1, n2, sd1, sd2, tol = .Machine$double.eps^0.5){
  uniroot(f = function(q, p, n1, n2, sd1, sd2){ pgt(q, n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2) - p },
          lower = -5, upper = 5, extendInt = "yes", 
          n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2, p = p, tol = tol)$root
}
xiao.t.test <- function(x, ...){
    UseMethod("xiao.t.test")
}
xiao.t.test.default <- function (x, y, alternative = c("two.sided", "less", "greater"),
                                mu = 0, conf.level = 0.95, ...){
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")

  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  yok <- !is.na(y)
  xok <- !is.na(x)
  y <- y[yok]
  x <- x[xok]
  nx <- length(x)
  if(nx < 2) stop("not enough 'x' observations")
  mx <- mean(x)
  vx <- var(x)
  ny <- length(y)
  if(ny < 2) stop("not enough 'y' observations")
  my <- mean(y)
  vy <- var(y)
  estimate <- c(mx, my, sqrt(vx), sqrt(vy))
  names(estimate) <- c("mean of x", "mean of y", "SD of x", "SD of y")
  s <- sqrt(vx + vy)
  se <- sqrt(vx/nx + vy/ny)
  if(nx*(nx-1)/vx <= ny*(ny-1)/vy){
    df <- (ny-1)*(1 + ny/nx*vx/vy)
  }else{
    df <- (nx-1)*(1 + nx/ny*vy/vx) 
  }
  if(se < 10 * .Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
  tstat <- (mx - my - mu)/se

  if(alternative == "less") {
    pval <- pgt(tstat, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy))
    cint <- c(-Inf, tstat + qgt(conf.level, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy)))
  }
  else if(alternative == "greater") {
    pval <- pgt(tstat, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy), lower.tail = FALSE)
    cint <- c(tstat - qgt(conf.level, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy)), Inf)
  }
  else{
    pval <- 2 * pgt(-abs(tstat), n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy))
    alpha <- 1 - conf.level
    cint <- qgt(1 - alpha/2, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy))
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * se
  names(tstat) <- "gt"
  names(df) <- "df"
  names(mu) <- "difference in means"
  n <- c(nx, ny)
  names(n) <- c("n of x", "n of y")
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               stderr = se, alternative = alternative, n = n,
               method = "Xiao Two Sample t-test",
               data.name = dname)
  class(rval) <- "htest"
  rval
}

xiao.t.test.formula <- function (formula, data, subset, na.action, ...){
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("xiao.t.test", c(DATA, list(...)))
  y$data.name <- DNAME
  if (length(y$estimate) == 2L)
    names(y$estimate) <- paste("mean in group", levels(g))
  y
}
