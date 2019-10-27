## unpaired case (Newcombe (1998), Statist. Med. 17, 873-890)
## a: number of successes in group 1
## b: number of successes in group 2
## c: number of failures in group 1
## d: number of failures in group 2
## paired case (Newcombe (1998), Statist. Med. 17, 2635-2650)
## a: number of successes in group 1 + 2
## b: number of sccesses in group 1 with failure in group 2
## c: number of failures in group 1 with successes in group 2
## d: number of failures in group 1 + 2
## conf.level: confidence level
## paired: logical
## method
binomDiffCI <- function(a, b, c, d, conf.level = 0.95, 
                        paired = FALSE, 
                        method = ifelse(paired, "wilson-cc", "wilson")){
  stopifnot(is.numeric(a), is.numeric(b), is.numeric(c), is.numeric(d),
            is.numeric(conf.level))
  stopifnot(length(a) == 1, length(b) == 1, length(c) == 1, 
            length(d) == 1, length(conf.level) == 1)
  stopifnot(a >= 0, b >= 0, c >= 0, d >= 0)
  if(!(0 < conf.level & conf.level < 1)) 
    stop("'conf.level' must be in (0,1)")
  
  a <- trunc(a)
  b <- trunc(b)
  c <- trunc(c)
  d <- trunc(d)
  alpha <- 1 - conf.level
  
  if(paired){
    METHODS <- c("wald", "wald-cc", "wilson", "wilson-cc")
    method <- pmatch(method, METHODS)
  }else{
    stopifnot(a+c > 0, b+d > 0)
    METHODS <- c("wald", "wald-cc", "wilson")
    method <- pmatch(method, METHODS)
  }
  if (is.na(method))
    stop("invalid method")
  if (method == -1)
    stop("ambiguous method")
  
  Infos <- NULL
  if(paired){
    n <- a + b + c + d
    p1 <- (a+b)/n
    p2 <- (a+c)/n
    D <- (b-c)/n
    names(D) <- "difference of proportions"
    if(method == 1){ # wald
      if(n < 30)
        warning("The sample size might be too low for the asymptotic method.")
      if(p1 < 0.1 | p2 < 0.1 | p1 > 0.9 | p2 > 0.9)
        warning("At least one of the observed proportions might be too extreme for the asymptotic method.")
      z <- qnorm(1-alpha/2)
      SED <- sqrt((b+c) - (b-c)^2/n)/n
      CI.lower <- max(-1, D - z*SED)
      CI.upper <- min(D + z*SED, 1)
      Infos <- list(SED, c(p1, p2))
      names(Infos[[1]]) <- "SE of difference of proportions"
      names(Infos[[2]]) <- c("proportion of group 1", "proportion of group 2")
    }
    if(method == 2){ # wald with continuity correction
      if(n < 30)
        warning("The sample size might be too low for the asymptotic method.")
      if(p1 < 0.1 | p2 < 0.1 | p1 > 0.9 | p2 > 0.9)
        warning("At least one of the observed proportions might be too extreme for the asymptotic method.")
      z <- qnorm(1-alpha/2)
      SED <- sqrt((b+c) - (b-c)^2/n)/n
      ## error in Newcombe (1998)
      CC <- 0.5/n
      CI.lower <- max(-1, D - z*SED + CC)
      CI.upper <- min(D + z*SED - CC, 1)
      Infos <- list(SED, c(p1, p2))
      names(Infos[[1]]) <- "SE of difference of proportions"
      names(Infos[[2]]) <- c("proportion of group 1", "proportion of group 2")
    }
    if(method == 3){ # wilson
      CI1 <- binomCI(a+b, n, method = "wilson")$conf.int
      l1 <- CI1[1]
      u1 <- CI1[2]
      CI2 <- binomCI(a+c, n, method = "wilson")$conf.int
      l2 <- CI2[1]
      u2 <- CI2[2]
      A <- (a+b)*(c+d)*(a+c)*(b+d)
      if(A == 0){
        phi <- 0
      }else{
        phi <- (a*d-b*c)/sqrt(A)
      }
      CI.lower <- D - sqrt((p1-l1)^2 - 2*phi*(p1-l1)*(u2-p2) + (u2-p2)^2)
      CI.upper <- D + sqrt((p2-l2)^2 - 2*phi*(p2-l2)*(u1-p1) + (u1-p1)^2)
      Infos <- c(p1, p2)
      names(Infos) <- c("proportion of group 1", "proportion of group 2")
    }
    if(method == 4){ # wilson with continuity correction
      CI1 <- binomCI(a+b, n, method = "wilson")$conf.int
      l1 <- CI1[1]
      u1 <- CI1[2]
      CI2 <- binomCI(a+c, n, method = "wilson")$conf.int
      l2 <- CI2[1]
      u2 <- CI2[2]
      A <- (a+b)*(c+d)*(a+c)*(b+d)
      if(A == 0){
        phi <- 0
      }else{
        B <- a*d - b*c
        if(B <= 0) C <- B
        if(B > 0) C <- max(B-n/2, 0)
        phi <- C/sqrt(A)
      }
      CI.lower <- D - sqrt((p1-l1)^2 - 2*phi*(p1-l1)*(u2-p2) + (u2-p2)^2)
      CI.upper <- D + sqrt((p2-l2)^2 - 2*phi*(p2-l2)*(u1-p1) + (u1-p1)^2)
      Infos <- c(p1, p2)
      names(Infos) <- c("proportion of group 1", "proportion of group 2")
    }
    nameCI <- "difference of proportions (paired data)"
    nameMethod <- paste(METHODS[method], 
                        "confidence interval (paired data)")
  }else{
    m <- a+c
    n <- b+d
    p1 <- a/m
    p2 <- b/n
    D <- p1 - p2
    names(D) <- "difference of proportions"
    if(method == 1){ ## Wald, simple asymptotic method
      if(m < 30 | n < 30)
        warning("The sample sizes might be too low for the asymptotic method.")
      if(p1 < 0.1 | p2 < 0.1 | p1 > 0.9 | p2 > 0.9)
        warning("At least one of the observed proportions might be too extreme for the asymptotic method.")
      z <- qnorm(1-alpha/2)
      SED <- sqrt(a*c/m^3 + b*d/n^3)
      CI.lower <- max(-1, D - z*SED)
      CI.upper <- min(D + z*SED, 1)
      Infos <- list(SED, c(p1, p2))
      names(Infos[[1]]) <- "SE of difference of proportions"
      names(Infos[[2]]) <- c("proportion of group 1", "proportion of group 2")
    }
    if(method == 2){ ## Wald, simple asymptotic method with continuity correction
      if(m < 30 | n < 30)
        warning("The sample sizes might be too low for the asymptotic method.")
      if(p1 < 0.1 | p2 < 0.1 | p1 > 0.9 | p2 > 0.9)
        warning("At least one of the observed proportions might be too extreme for the asymptotic method.")
      z <- qnorm(1-alpha/2)
      ## sign error in Newcombe (1998)
      SED <- sqrt(a*c/m^3 + b*d/n^3)
      CC <- 0.5*(1/m + 1/n)
      CI.lower <- max(-1, D - z*SED + CC)
      CI.upper <- min(D + z*SED - CC, 1)
      Infos <- list(SED, c(p1, p2))
      names(Infos[[1]]) <- "SE of difference of proportions"
      names(Infos[[2]]) <- c("proportion of group 1", "proportion of group 2")
    }
    if(method == 3){ ## Wilson
      CI1 <- binomCI(a, m, method = "wilson")$conf.int
      l1 <- CI1[1]
      u1 <- CI1[2]
      CI2 <- binomCI(b, n, method = "wilson")$conf.int
      l2 <- CI2[1]
      u2 <- CI2[2]
      CI.lower <- D - sqrt((p1-l1)^2 + (u2-p2)^2)
      CI.upper <- D + sqrt((p2-l2)^2 + (u1-p1)^2)
      Infos <- c(p1, p2)
      names(Infos) <- c("proportion of group 1", "proportion of group 2")
    }
    nameCI <- "difference of independent proportions"
    nameMethod <- paste(METHODS[method], 
                        "confidence interval (independent proportions)")
  }

  CI <- matrix(c(CI.lower, CI.upper), nrow = 1)
  rownames(CI) <- nameCI
  colnames(CI) <- c(paste(alpha/2*100, "%"),
                    paste((1-alpha/2)*100, "%"))
  attr(CI, "conf.level") <- conf.level
  
  return(structure(list("estimate" = D, "conf.int" = CI, "Infos" = Infos,
                        method = nameMethod),
                   class = "confint"))
}
