md2sens <- function(delta, sd1 = NULL, sd2 = NULL){
  stopifnot(!missing(delta))
  stopifnot(is.numeric(delta))
  stopifnot(is.numeric(sd1))
  stopifnot(is.numeric(sd2))
  if(length(delta) != 1){
    stop("'delta' must be single number (vector of length 1)!")
  }
  if(delta < 0){
    stop("'delta' must be a non-negative real number!")
  } 
  
  var.equal <- ifelse(abs(sd1-sd2) < 1e-8, TRUE, FALSE)
  
  if(var.equal){
    sd <- sd1
    if(length(sd) != 1){
      stop("'sd' must be single number (vector of length 1)!")
    }
    if(is.null(sd)){
      stop("'sd' must be a positive real number")
    }
    spec <- sens <- pnorm(delta/2, mean = 0, sd = sd)
  }else{
    if(length(sd1) != 1 || length(sd2) != 1){
      stop("'sd1' and 'sd2' must single numbers (vectors of length 1)!")
    }
    if(is.null(sd1) || is.null(sd2)){
      stop("'sd1' and 'sd2' must be positive real numbers!")
    }
    Min <- min(-3*sd1, delta - 3*sd2)
    Max <- max(3*sd1, delta + 3*sd2)
    diff <- function(x, delta, sd1, sd2){
      dnorm(x, mean = 0, sd = sd1) - dnorm(x, mean = delta, sd = sd2)
    }
    if(sd1 < sd2){
      cutoff <- uniroot(diff, interval = c(0, Max), 
                        delta = delta, sd1 = sd1, sd2 = sd2)$root
    }else{
      cutoff <- uniroot(diff, interval = c(Min, delta), 
                        delta = delta, sd1 = sd1, sd2 = sd2)$root
    }
    sens <- pnorm(cutoff, mean = 0, sd = sd1)
    spec <- pnorm(cutoff, mean = delta, sd = sd2, lower.tail = FALSE)
  }
  c("sensitivity" = sens, "specificity" = spec)
}
