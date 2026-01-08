rmx.t.test <- function(x, ...){
    UseMethod("rmx.t.test")
}
rmx.t.test.default <- function (x, y = NULL, alternative = c("two.sided", "less", "greater"),
                                mu = 0, paired = FALSE, conf.level = 0.95, 
                                method = c("welch", "xiao", "hsu", "student"), 
                                eps.lower = 0, eps.upper = 0.05, k = 3L, ...){
  alternative <- match.arg(alternative)
  
  method <- match.arg(method)
  
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  
  
  if( !is.null(y) ) {
    dname <- paste(deparse1(substitute(x)),"and", deparse1(substitute(y)))
    if(paired)
      xok <- yok <- complete.cases(x,y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    dname <- deparse1(substitute(x))
    if (paired) stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if (paired) {
    x <- x-y
    y <- NULL
  }
  
  nx <- length(x)
  estx <- robrmx:::rmx.norm(x, eps.lower = eps.lower, eps.upper = eps.upper, k = k)
  mx <- estx$rmxEst[1]
  vx <- estx$rmxEst[2]^2
  
  if(is.null(y)) {
    if(nx < 2) stop("not enough 'x' observations")
    df <- nx-1
    se <- sqrt(vx/nx)
    if(!is.na(se) && se < 10 *.Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    tstat <- (mx-mu)/se
    METHOD <- if(paired) "RMX Paired t-test" else "RMX One Sample t-test"
    estimate <- c(mx, sqrt(vx))
    if(paired){
      names(estimate) <- c("mean difference", "SD of difference")
    }else{
      names(estimate) <- c("mean of x", "SD of x")
    }
  }else{
    ny <- length(y)
    if(nx < 1 || (method != "student" && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (method != "student" && ny < 2))
      stop("not enough 'y' observations")
    if(method == "student" && nx+ny < 3) stop("not enough observations")
    esty <- robrmx:::rmx.norm(y, eps.lower = eps.lower, eps.upper = eps.upper, k = k)
    my <- esty$rmxEst[1]
    vy <- esty$rmxEst[2]^2
    
    estimate <- c(mx, my, sqrt(vx), sqrt(vy))
    names(estimate) <- c("mean of x", "mean of y", "SD of x", "SD of y")
    
    if(method == "student"){
      df <- nx + ny - 2
      v <- 0
      if (nx > 1) 
        v <- v + (nx - 1) * vx
      if (ny > 1) 
        v <- v + (ny - 1) * vy
      v <- v/df
      se <- sqrt(v * (1/nx + 1/ny))
      METHOD <- "RMX Student Two Sample t-test"
    }else{
      se <- sqrt(vx/nx + vy/ny)
      if(method == "welch"){
        df <- se^4/(se^4/(nx - 1) + se^4/(ny - 1))
        METHOD <- "RMX Welch Two Sample t-test"
      }
      if(method == "hsu"){
        df <- min(nx, ny) - 1
        METHOD <- "RMX Hsu Two Sample t-test"
      }
      if(method == "xiao"){
        if(nx*(nx-1)/vx <= ny*(ny-1)/vy){
          df <- (ny-1)*(1 + ny/nx*vx/vy)
        }else{
          df <- (nx-1)*(1 + nx/ny*vy/vx) 
        }
        METHOD <- "RMX Xiao Two Sample t-test"
      }
    }
    if(se < 10 * .Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/se
  }
  
  if(alternative == "less") {
    if(method == "xiao"){
      pval <- pgt(tstat, n1 = nx, n2 = ny, v1tov2 = vx/vy)
      cint <- c(-Inf, tstat + qgt(conf.level, n1 = nx, n2 = ny, v1tov2 = vx/vy))
    }else{
      pval <- pt(tstat, df)
      cint <- c(-Inf, tstat + qt(conf.level, df))
    }
  }
  else if(alternative == "greater") {
    if(method == "xiao"){
      pval <- pgt(tstat, n1 = nx, n2 = ny, v1tov2 = vx/vy, lower.tail = FALSE)
      cint <- c(tstat - qgt(conf.level, n1 = nx, n2 = ny, v1tov2 = vx/vy), Inf)
    }else{
      pval <- pt(tstat, df, lower.tail = FALSE)
      cint <- c(tstat - qt(conf.level, df), Inf)
    }
  }
  else{
    alpha <- 1 - conf.level
    if(method == "xial"){
      pval <- 2 * pgt(-abs(tstat), n1 = nx, n2 = ny, v1tov2 = vx/vy)
      cint <- qgt(1 - alpha/2, n1 = nx, n2 = ny, v1tov2 = vx/vy)
    }else{
      pval <- 2 * pt(-abs(tstat), df)
      cint <- qt(1 - alpha/2, df)
    }
    cint <- tstat + c(-cint, cint)
  }
  
  cint <- mu + cint * se
  if(method == "xiao"){
    names(tstat) <- "gt"
  }else{
    names(tstat) <- "t"
  }
  names(df) <- "df"
  if(paired){ 
    names(mu) <- "mean difference" 
  } else if(!is.null(y)){
    names(mu) <- "difference in means"
  } else {
    names(mu) <- "mean"
  }
  if(is.null(y)){
    n <- nx
    names(n) <- "n of x"
  }else{
    n <- c(nx, ny)
    names(n) <- c("n of x", "n of y")
  }
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               stderr = se, alternative = alternative, n = n,
               method = METHOD, data.name = dname)
  class(rval) <- "htest"
  rval
}

rmx.t.test.formula <- function (formula, data, subset, na.action = na.pass, ...) {
  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  if ("paired" %in% ...names())
    stop("cannot use 'paired' in formula method")    
  oneSampleOrPaired <- FALSE
  if (length(attr(terms(formula[-2L]), "term.labels")) != 1L) 
    if (formula[[3L]] == 1L)
      oneSampleOrPaired <- TRUE
  else
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ") # works in all cases
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  if (! oneSampleOrPaired) {
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L) 
      stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    ## Call the default method.
    y <- rmx.t.test(x = DATA[[1L]], y = DATA[[2L]], ...)
    if (length(y$estimate) == 4L) {
      names(y$estimate)[1:2] <- paste("mean in group", levels(g))
      names(y$estimate)[3:4] <- paste("SD of group", levels(g))
      names(y$null.value) <-
        paste("difference in means between",
              paste("group", levels(g), collapse = " and "))
    }
  }
  else { # 1-sample and paired tests
    respVar <- mf[[response]]
    if (inherits(respVar, "Pair")) {
      ## Call the default method.
      y <- rmx.t.test(x = respVar[, 1L], y = respVar[, 2L],
                      paired = TRUE, ...)
    }
    else {
      ## Call the default method.            
      y <- rmx.t.test(x = respVar, ...)
    }
  }
  y$data.name <- DNAME
  y
}
