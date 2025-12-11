md2zfactor <- function(delta, sd1 = 1, sd2 = 1){
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
  if(length(sd1) != 1 || length(sd2) != 1){
    stop("'sd1' and 'sd2' must single numbers (vectors of length 1)!")
  }
  if(is.null(sd1) || is.null(sd2)){
    stop("'sd1' and 'sd2' must be positive real numbers!")
  }
  if(sd1 <= 0 || sd2 <= 0){
    stop("The values of 'sd1' and 'sd2' must be positive real numbers!")
  }
  
  1 - 3*(sd1 + sd2)/delta
}
