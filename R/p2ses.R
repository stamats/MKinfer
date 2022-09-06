p2ses <- function(p, alternative = c("two.sided", "less", "greater")){
  alternative <- match.arg(alternative)
  stopifnot(all(p <= 1))
  stopifnot(all(p >= 0))
  
  if(alternative == "less"){
    ses <- qnorm(p)  
  }else if (alternative == "greater"){
    ses <- qnorm(p, lower.tail = FALSE)
  }else{
    ses <- qnorm(p/2, lower.tail = FALSE) 
  }
  ses
}
