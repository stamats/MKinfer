\name{normCI}
\alias{normCI}
\alias{meanCI}
\alias{sdCI}
\title{ Confidence Intervals for Mean and Standard Deviation }
\description{
  This function can be used to compute confidence intervals for mean and
  standard deviation of a normal distribution.
}
\usage{
normCI(x, mean = NULL, sd = NULL, conf.level = 0.95, 
       boot = FALSE, R = 9999, bootci.type = "all", na.rm = TRUE, 
       alternative = c("two.sided", "less", "greater"), ...)
meanCI(x, conf.level = 0.95, boot = FALSE, R = 9999, bootci.type = "all", 
       na.rm = TRUE, alternative = c("two.sided", "less", "greater"), ...)
sdCI(x, conf.level = 0.95, boot = FALSE, R = 9999, bootci.type = "all", 
     na.rm = TRUE, alternative = c("two.sided", "less", "greater"), ...)
}
\arguments{
  \item{x}{ vector of observations. }
  \item{mean}{ mean if known otherwise \code{NULL}. }
  \item{sd}{ standard deviation if known otherwise \code{NULL}. }
  \item{conf.level}{ confidence level. }
  \item{boot}{ a logical value indicating whether bootstrapped confidence intervals shall be computed. }
  \item{R}{ number of bootstrap replicates. }
  \item{bootci.type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
  \item{na.rm}{ a logical value indicating whether NA values should be stripped before the computation proceeds. }
  \item{alternative}{ a character string specifying one- or two-sided confidence intervals. 
  Must be one of "two.sided" (default), "greater" or "less" (one-sided intervals). 
  You can specify just the initial letter. }
  \item{\dots}{ further arguments passed to function \code{\link[boot]{boot}},
  e.g. for parallel computing. }
}
\details{
  The standard confidence intervals for mean and standard deviation are
  computed that can be found in many textbooks, e.g. Chapter 4 in Altman et al. (2000).
  
  In addition, bootstrap confidence intervals for mean and/or SD can be computed,
  where function \code{\link[boot]{boot.ci}} is applied.
}
\value{
  A list with class \code{"confint"} containing the following components:
  \item{estimate}{ the estimated mean and sd. }
  \item{conf.int}{ confidence interval(s) for mean and/or sd. }
  \item{Infos}{ additional information. }
}
\references{
  D. Altman, D. Machin, T. Bryant, M. Gardner (eds). Statistics with Confidence:
  Confidence Intervals and Statistical Guidelines, 2nd edition 2000.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
\examples{
x <- rnorm(50)
## mean and sd unknown
normCI(x)
meanCI(x)
sdCI(x)

## one-sided
normCI(x, alternative = "less")
meanCI(x, alternative = "greater")
sdCI(x, alternative = "greater")

## bootstrap intervals (R = 999 to reduce computation time for R checks)
normCI(x, boot = TRUE, R = 999)
meanCI(x, boot = TRUE, R = 999)
sdCI(x, boot = TRUE, R = 999)

normCI(x, boot = TRUE, R = 999, alternative = "less")
meanCI(x, boot = TRUE, R = 999, alternative = "less")
sdCI(x, boot = TRUE, R = 999, alternative = "greater")

## sd known
normCI(x, sd = 1)
## bootstrap intervals only for mean (sd is ignored)
## (R = 999 to reduce computation time for R checks)
normCI(x, sd = 1, boot = TRUE, R = 999)

## mean known
normCI(x, mean = 0)
## bootstrap intervals only for sd (mean is ignored)
## (R = 999 to reduce computation time for R checks)
normCI(x, mean = 0, boot = TRUE, R = 999)

\donttest{
## parallel computing for bootstrap
normCI(x, boot = TRUE, R = 9999, parallel = "multicore", ncpus = 2)
}
}
\keyword{univar}
