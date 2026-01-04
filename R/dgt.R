## Gauss-Legendre integration if integrate fails
.GLInt <- function(f, lower, upper, ...){
  AW <- getFromNamespace(".AW.500", ns = "MKinfer")
  
  # transformation to [lower, upper]
  xl <- (upper - lower)/2
  W <- xl*AW[,2]
  A <- xl*AW[,1] + (lower + upper)/2
  
  sum(W*c(f(A, ...)))
}
## generalized central t distribution as defined in Section 5.2 in Xiao (2018)
dgt <- function(x, n1, n2, v1tov2, log = FALSE, MAX = 1e4){
  T1 <- (n2-1)*(1 + n2/n1*v1tov2)
  T2 <- (n1-1)/(n2-1)*n1/n2*(1/v1tov2)
  
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

pgt <- function(q, n1, n2, v1tov2, lower.tail = TRUE, log.p = FALSE, 
                rel.tol = .Machine$double.eps^0.25, parallel = FALSE, 
                cl = NULL){
  if(max(c(length(q), length(n1), length(n2), length(v1tov2))) == 1){
    if(lower.tail){
      res <- try(integrate(f = dgt, lower = -Inf, upper = q, n1 = n1, n2 = n2, 
                           v1tov2 = v1tov2, rel.tol = rel.tol)$value,
                 silent = TRUE)
      if(inherits(res, "try-error")){
        res <- .GLInt(f = dgt, lower = q-20, upper = q, n1 = n1, n2 = n2, 
                      v1tov2 = v1tov2)
      }
    }else{
      res <- try(integrate(f = dgt, lower = q, upper = Inf, n1 = n1, n2 = n2, 
                           v1tov2, rel.tol = rel.tol)$value,
                 silent = TRUE)
      if(inherits(res, "try-error")){
        res <- .GLInt(f = dgt, lower = q, upper = q+20, n1 = n1, n2 = n2, 
                      v1tov2 = v1tov2)
      }
    }
  }else{
    ARGS <- cbind(q, n1, n2, v1tov2)
    dfun <- function(args, lower.tail, rel.tol){
      if(lower.tail){
        res <- try(integrate(f = dgt, lower = -Inf, upper = args[1], n1 = args[2], 
                             n2 = args[3], v1tov2 = args[4], rel.tol = rel.tol)$value,
                   silent = TRUE)
        if(inherits(res, "try-error")){
          res <- .GLInt(f = dgt, lower = args[1]-20, upper = args[1], n1 = args[2], 
                        n2 = args[3], v1tov2 = args[4])
        }
      }else{
        res <- try(integrate(f = dgt, lower = args[1], upper = Inf, n1 = args[2], 
                             n2 = args[3], v1tov2 = args[4], rel.tol = rel.tol)$value,
                   silent = TRUE)
        if(inherits(res, "try-error")){
          res <- .GLInt(f = dgt, lower = args[1], upper = args[1]+20, n1 = args[2], 
                        n2 = args[3], v1tov2 = args[4])
        }
      }
      res
    }
    if(parallel){
      if(is.null(cl)){
        ncores <- detectCores()-1
        cl <- makePSOCKcluster(rep("localhost", ncores))
        clusterExport(cl, list("dgt"), 
                      envir = as.environment("package:MKinfer"))
        res <- parRapply(cl = cl, x = ARGS, FUN = dfun, lower.tail = lower.tail, 
                         rel.tol = rel.tol)
        stopCluster(cl)
      }else{
        clusterExport(cl, list("dgt"), 
                      envir = as.environment("package:MKinfer"))
        res <- parRapply(cl = cl, x = ARGS, FUN = dfun, lower.tail = lower.tail, 
                         rel.tol = rel.tol)
      }
    }else{
      res <- apply(ARGS, 1, dfun, lower.tail = lower.tail, rel.tol = rel.tol)
    }
  }
  if(log.p) res <- log(res)
  res
}

qgt <- function(p, n1, n2, v1tov2, lower.tail = TRUE, log.p = FALSE, 
                tol = .Machine$double.eps^0.5, parallel = FALSE, cl = NULL){
  if(any(p <= 0)) stop("p must be positive.")
  if(any(p >= 1)) stop("p must be smaller than 1.")
  MIN <- -15
  MAX <- 15
  if(log.p) p <- exp(p)
  fun <- function(q, p, n1, n2, v1tov2, lower.tail){ 
    pgt(q, n1 = n1, n2 = n2, v1tov2, lower.tail = lower.tail) - p 
  }
  if(max(c(length(p), length(n1), length(n2), length(v1tov2))) == 1){
    res <- try(uniroot(f = fun, lower = MIN, upper = MAX, extendInt = "yes", 
                       n1 = n1, n2 = n2, v1tov2 = v1tov2, p = p, 
                       lower.tail = lower.tail, tol = tol)$root,
               silent = TRUE)
    if(inherits(q, "try-error")){
      res <- NA
    }
  }else{
    ARGS <- cbind(p, n1, n2, v1tov2)
    qfun <- function(args, lower.tail, tol, MIN, MAX){
      q <- try(uniroot(f = fun, lower = MIN, upper = MAX, extendInt = "yes", 
                       n1 = args[2], n2 = args[3], v1tov2 = args[4], 
                       p = args[1], lower.tail = lower.tail, tol = tol)$root,
               silent = TRUE)
      if(inherits(q, "try-error")){
        q <- NA
      }
      q
    }
    if(parallel){
      if(is.null(cl)){
        ncores <- detectCores()-1
        cl <- makePSOCKcluster(rep("localhost", ncores))
        clusterExport(cl, list("pgt"), 
                      envir = as.environment("package:MKinfer"))
        res <- parRapply(cl = cl, x = ARGS, FUN = qfun, lower.tail = lower.tail, 
                         tol = tol, MIN = MIN, MAX = MAX)
        stopCluster(cl)
      }else{
        clusterExport(cl, list("pgt"), 
                      envir = as.environment("package:MKinfer"))
        res <- parRapply(cl = cl, x = ARGS, FUN = qfun, lower.tail = lower.tail, 
                         tol = tol, MIN = MIN, MAX = MAX)
      }
    }else{
      res <- apply(ARGS, 1, qfun, lower.tail = lower.tail, tol = tol, MIN = MIN, MAX = MAX)
    }
  }
  res
}
