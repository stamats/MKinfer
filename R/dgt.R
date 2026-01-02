## generalized central t distribution as defined in Section 5.2 in Xiao (2018)
dgt <- function(x, n1, n2, sd1, sd2, log = FALSE, MAX = 1e4){
  T1 <- (n2-1)*(1 + n2/n1*sd1^2/sd2^2)
  T2 <- (n1-1)*n1*sd2^2/((n2-1)*n2*sd1^2)
  
  maxiter <- 2000
  suppressWarnings({
  T3 <- try(Re(hypergeo((n1+n2-1)/2, (n1-1)/2, (n1+n2-2)/2, 1 - (T2+x^2/T1)/(1+x^2/T1), 
                        maxiter = maxiter)), 
            silent = TRUE)
  while(inherits(T3, "try-error") || any(is.infinite(T3))){
    maxiter <- 2*maxiter
    T3 <- try(Re(hypergeo((n1+n2-1)/2, (n1-1)/2, (n1+n2-2)/2, 1 - (T2+x^2/T1)/(1+x^2/T1), 
                          maxiter = maxiter)), 
              silent = TRUE)
    if(maxiter > MAX) break
  }})
  if(inherits(T3, "try-error") || any(is.infinite(T3)))
    stop("Function hypergeo did not converge! You can try increasing 'MAX'.")
  
  res <- 1/sqrt(pi*T1)*gamma((n1+n2-1)/2)/gamma((n1+n2-2)/2)*(1+x^2/T1)^(-(n1+n2-1)/2)*T2^((n1-1)/2)*T3
  if(log) res <- log(res)
  res
}

pgt <- function(q, n1, n2, sd1, sd2, lower.tail = TRUE, log.p = FALSE, 
                rel.tol = .Machine$double.eps^0.25){
  if(max(c(length(q), length(n1), length(n2), length(sd1), length(sd2))) == 1){
    if(lower.tail){
      res <- integrate(f = dgt, lower = -Inf, upper = q, n1 = n1, n2 = n2, 
                       sd1 = sd1, sd2 = sd2, rel.tol = rel.tol)$value
    }else{
      res <- integrate(f = dgt, lower = q, upper = Inf, n1 = n1, n2 = n2, 
                       sd1 = sd1, sd2 = sd2, rel.tol = rel.tol)$value
    }
  }else{
    ARGS <- cbind(q, n1, n2, sd1, sd2)
    dfun <- function(args, lower.tail, rel.tol){
      if(lower.tail){
        res <- integrate(f = dgt, lower = -Inf, upper = args[1], n1 = args[2], 
                         n2 = args[3], sd1 = args[4], sd2 = args[5], rel.tol = rel.tol)$value
      }else{
        res <- integrate(f = dgt, lower = args[1], upper = Inf, n1 = args[2], 
                         n2 = args[3], sd1 = args[4], sd2 = args[5], rel.tol = rel.tol)$value
      }
    }
    res <- apply(ARGS, 1, dfun, lower.tail = lower.tail, rel.tol = rel.tol)
  }
  if(log.p) res <- log(res)
  res
}

qgt <- function(p, n1, n2, sd1, sd2, lower.tail = TRUE, log.p = FALSE, tol = .Machine$double.eps^0.5){
  MIN <- -15
  MAX <- 15
  if(log.p) p <- exp(p)
  fun <- function(q, p, n1, n2, sd1, sd2, lower.tail){ 
    pgt(q, n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2, lower.tail = lower.tail) - p 
  }
  if(max(c(length(p), length(n1), length(n2), length(sd1), length(sd2))) == 1){
    res <- uniroot(f = fun, lower = MIN, upper = MAX, extendInt = "yes", 
                   n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2, p = p, 
                   lower.tail = lower.tail, tol = tol)$root
  }else{
    ARGS <- cbind(p, n1, n2, sd1, sd2)
    qfun <- function(args, lower.tail, tol, MIN, MAX){
      uniroot(f = fun, lower = MIN, upper = MAX, extendInt = "yes", 
              n1 = args[2], n2 = args[3], sd1 = args[4], sd2 = args[5], 
              p = args[1], lower.tail = lower.tail, tol = tol)$root
    }
    res <- apply(ARGS, 1, qfun, lower.tail = lower.tail, tol = tol, MIN = MIN, MAX = MAX)
  }
  res
}
