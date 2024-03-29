\name{quantileCI}
\alias{quantileCI}
\alias{medianCI}
\alias{madCI}
\title{ Confidence Intervals for Quantiles }
\description{
  These functions can be used to compute confidence intervals for quantiles
  (including median and MAD).
}
\usage{
quantileCI(x, prob = 0.5, conf.level = 0.95, method = "exact",
           R = 9999, bootci.type = c("norm", "basic", "perc", "bca"), 
           na.rm = FALSE, alternative = c("two.sided", "less", "greater"), ...)
medianCI(x, conf.level = 0.95, method = "exact",
         R = 9999, bootci.type = c("norm", "basic", "perc", "bca"), 
         na.rm = FALSE, alternative = c("two.sided", "less", "greater"), ...)
madCI(x, conf.level = 0.95, method = "exact", 
      R = 9999, bootci.type = c("norm", "basic", "perc", "bca"), 
      na.rm = FALSE, constant = 1.4826, 
      alternative = c("two.sided", "less", "greater"), ...)
}
\arguments{
  \item{x}{ numeric data vector }
  \item{prob}{ quantile }
  \item{conf.level}{ confidence level }
  \item{method}{ character string specifing which method to use; see details. }
  \item{R}{ number of bootstrap replicates. }
  \item{bootci.type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
  \item{na.rm}{ logical, remove \code{NA} values. }
  \item{constant}{ scale factor (see \code{\link[stats]{mad}}). }
  \item{alternative}{ a character string specifying one- or two-sided confidence intervals. 
  Must be one of "two.sided" (default), "greater" or "less" (one-sided intervals). 
  You can specify just the initial letter. }
  \item{\dots}{ further arguments passed to function \code{\link[boot]{boot}},
  e.g. for parallel computing. }
}
\details{
  The exact confidence interval (\code{method = "exact"}) is computed using binomial
  probabilities; see Section 6.8.1 in Sachs and Hedderich (2009). If the result is not
  unique, i.e. there is more than one interval with coverage proability closest to
  \code{conf.level}, the shortest confidence interval is returned.

  The asymptotic confidence interval (\code{method = "asymptotic"}) is based on the
  normal approximation of the binomial distribution; see Section 6.8.1 in Sachs and Hedderich (2009).
  
  In case of discrete data, there are alternative bootstrap approaches that might 
  give better results; see Jentsch and Leucht (2016).
  
  Since \code{madCI} is computed as the median confidence interval of the 
  absolut deviations from the sample median and ignores the variablity of the 
  sample median, the exact and asymptotic confidence intervals might be too short.
  We recommend to use bootstrap confidence intervals.
}
\value{
  A list with components
  \item{estimate}{ the sample quantile. }
  \item{CI}{ a confidence interval for the sample quantile. }
}
\references{
  L. Sachs and J. Hedderich (2009). \emph{Angewandte Statistik}. Springer.
  
  C. Jentsch and A. Leucht (2016). Bootstrapping sample quantiles of discrete data.
  \emph{Ann Inst Stat Math}, \bold{68}: 491-539.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[stats]{quantile}}}
\examples{
## To get a non-trivial exact confidence interval for the median
## one needs at least 6 observations
x <- rnorm(5)
medianCI(x)

## asymptotic confidence interval
medianCI(x, method = "asymptotic")
madCI(x, method = "asymptotic")

## bootstrap confidence interval
x <- rnorm(50)
medianCI(x)
medianCI(x, method = "asymptotic")
## (R = 999 to reduce computation time for R checks)
medianCI(x, method = "boot", R = 999)

madCI(x)
madCI(x, method = "asymptotic")
## (R = 999 to reduce computation time for R checks)
madCI(x, method = "boot", R = 999)

## confidence interval for quantiles
quantileCI(x, prob = 0.25)
quantileCI(x, prob = 0.25, method = "asymptotic")

quantileCI(x, prob = 0.75)
## (R = 999 to reduce computation time for R checks)
quantileCI(x, prob = 0.75, method = "boot", R = 999)

## one-sided
quantileCI(x, prob = 0.75, alternative = "greater")
medianCI(x, alternative = "less", method = "asymptotic")
madCI(x, alternative = "greater", method = "boot", R = 999)

\donttest{
## parallel computing for bootstrap
medianCI(x, method = "boot", R = 9999, parallel = "multicore", 
         ncpus = 2)
}
}
\keyword{univar}
