\name{binomDiffCI}
\alias{binomDiffCI}
\title{ Confidence Intervals for Difference of Two Binomial Proportions }
\description{
  This function can be used to compute confidence intervals for the difference 
  of two binomial proportions. It includes methods for the independent and the 
  paired case.
}
\usage{
binomDiffCI(a, b, c, d, conf.level = 0.95, paired = FALSE, 
            method = ifelse(paired, "wilson-cc", "wilson"), 
            R = 9999, type = "all")
}
\arguments{
  \item{a}{ independent: number of successes of group 1; 
            paired: number of cases with success in group 1 and 2. }
  \item{b}{ independent: number of successes of group 2; 
            paired: number of cases with success in group 1 and failure in group 2. }
  \item{c}{ independent: number of failures of group 1; 
            paired: number of cases with failure in group 1 and success in group 2. }
  \item{d}{ independent: number of failures of group 2; 
            paired: number of cases with failure in group 1 and 2. }
  \item{conf.level}{ confidence level }
  \item{paired}{ a logical value indicating whether the two groups are paired. }
  \item{method}{ character string specifing which method to use; see details. }
  \item{R}{ number of bootstrap replicates. }
  \item{type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
}
\details{
  The Wald intervals (independent and paired) are obtained by applying the normal 
  approximation. There are also Wald intervals with continuity correction.

  The Wilson intervals are recommended by Newcombe and Altman (2000); see
  Chapter 6 of Altman et al. (2000). In the paired case, the continuity corrected
  version of the interval is recommended. The intervals are proposed in 
  Newcombe (1998a) and Newcombe (1998b).
  
  The bootstrap interval is calculated by using function \code{\link[boot]{boot.ci}}.
}
\value{
  A list with class \code{"confint"} containing the following components:
  \item{estimate}{ the estimated probability of success. }
  \item{conf.int}{ a confidence interval for the probability of success. }
}
\references{
  D. Altman, D. Machin, T. Bryant, M. Gardner (eds). \emph{Statistics with Confidence:
  Confidence Intervals and Statistical Guidelines}, 2nd edition.
  John Wiley and Sons 2000.
  
  R.G. Newcombe (1998a). Interval estimation for the difference between independent 
  proportions: comparison of eleven methods. \emph{Stat Med}, \bold{17}(8), 873-890.
  
  R.G. Newcombe (1998b). Improved confidence intervals for the difference between 
  binomial proportions based on paired data. \emph{Stat Med}, \bold{17}(22), 2635-2650.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[stats]{prop.test}}, \code{\link[boot]{boot.ci}}}
\examples{
## Example 1: Altman et al. (2000, p. 49)
## the confidence interval computed by prop.test
prop.test(c(63, 38), c(93, 92))$conf.int
## wald / simple asymptotic interval
binomDiffCI(a = 63, b = 38, c = 30, d = 54, method = "wald")
## wald / simple asymptotic interval with continuity correction
binomDiffCI(a = 63, b = 38, c = 30, d = 54, method = "wald-cc")
## wilson
binomDiffCI(a = 63, b = 38, c = 30, d = 54)
## bootstrap intervals (R = 999 to reduce computation time for R checks)
binomDiffCI(a = 63, b = 38, c = 30, d = 54, method = "boot", R = 999)

## Example 2: Altman et al. (2000, p. 50)
## the confidence interval computed by prop.test
prop.test(c(5, 0), c(56, 29))$conf.int
## wald / simple asymptotic interval
binomDiffCI(a = 5, b = 0, c = 51, d = 29, method = "wald")
## wald / simple asymptotic interval with continuity correction
binomDiffCI(a = 5, b = 0, c = 51, d = 29, method = "wald-cc")
## wilson
binomDiffCI(a = 5, b = 0, c = 51, d = 29)
## bootstrap intervals (R = 999 to reduce computation time for R checks)
binomDiffCI(a = 5, b = 0, c = 51, d = 29, method = "boot", R = 999)


## Example 3: Altman et al. (2000, p. 51)
## wald / simple asymptotic interval
binomDiffCI(a = 14, b = 5, c = 0, d = 22, paired = TRUE, method = "wald")
## wald / simple asymptotic interval with continuity correction
binomDiffCI(a = 14, b = 5, c = 0, d = 22, paired = TRUE, method = "wald-cc")
## wilson
binomDiffCI(a = 14, b = 5, c = 0, d = 22, paired = TRUE, method = "wilson")
## wilson with continuity correction
binomDiffCI(a = 14, b = 5, c = 0, d = 22, paired = TRUE)
## bootstrap intervals (R = 999 to reduce computation time for R checks)
binomDiffCI(a = 14, b = 5, c = 0, d = 22, paired = TRUE, method = "boot", R = 999)

## Example 4: Altman et al. (2000, p. 51)
## wald / simple asymptotic interval
binomDiffCI(a = 212, b = 144, c = 256, d = 707, paired = TRUE, method = "wald")
## wald / simple asymptotic interval with continuity correction
binomDiffCI(a = 212, b = 144, c = 256, d = 707, paired = TRUE, method = "wald-cc")
## wilson
binomDiffCI(a = 212, b = 144, c = 256, d = 707, paired = TRUE, method = "wilson")
## wilson with continuity correction
binomDiffCI(a = 212, b = 144, c = 256, d = 707, paired = TRUE)
## bootstrap intervals (R = 999 to reduce computation time for R checks)
binomDiffCI(a = 212, b = 144, c = 256, d = 707, paired = TRUE, method = "boot",
            type = c("norm", "basic", "stud", "perc"), R = 999) ## type = "bca" gives error
}
\keyword{univar}
