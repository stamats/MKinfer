\name{normDiffCI}
\alias{normDiffCI}
\alias{meanDiffCI}
\title{ Confidence Intervals for Difference of Means }
\description{
  This function can be used to compute confidence intervals for difference
  of means assuming normal distributions.
}
\usage{
normDiffCI(x, y, conf.level = 0.95, paired = FALSE, method = "welch", 
           boot = FALSE, R = 9999, bootci.type = "all", na.rm = TRUE, 
           alternative = c("two.sided", "less", "greater"), ...)
meanDiffCI(x, y, conf.level = 0.95, paired = FALSE, method = "welch", 
           boot = FALSE, R = 9999, bootci.type = "all", na.rm = TRUE, 
           alternative = c("two.sided", "less", "greater"), ...)
}
\arguments{
  \item{x}{ numeric vector of data values of group 1. }
  \item{y}{ numeric vector of data values of group 2. }
  \item{conf.level}{ confidence level. }
  \item{paired}{ a logical value indicating whether the two groups are paired. }
  \item{method}{ a character string specifing which method to use in the unpaired case; see details. }
  \item{boot}{ a logical value indicating whether bootstrapped confidence intervals shall be computed. }
  \item{R}{ number of bootstrap replicates. }
  \item{bootci.type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be stripped before the computation proceeds. }
  \item{alternative}{ a character string specifying one- or two-sided confidence intervals. 
  Must be one of "two.sided" (default), "greater" or "less" (one-sided intervals). 
  You can specify just the initial letter. }
  \item{\dots}{ further arguments passed to function \code{\link[boot]{boot}},
  e.g. for parallel computing. }
}
\details{
  The standard confidence intervals for the difference of means are
  computed that can be found in many textbooks, e.g. Chapter 4 in Altman et al. (2000).

  The method \code{"classical"} assumes equal variances whereas methods
  \code{"welch"} and \code{"hsu"} allow for unequal variances. The latter two
  methods use different formulas for computing the degrees of freedom of the
  respective t-distribution providing the quantiles in the confidence interval.
  Instead of the Welch-Satterhwaite equation the method of Hsu uses the minimum
  of the group sample sizes minus 1; see Section 6.8.3 of Hedderich and Sachs (2018).
}
\value{
  A list with class \code{"confint"} containing the following components:
  \item{estimate}{ point estimate (mean of differences or difference in means). }
  \item{conf.int}{ confidence interval. }
  \item{Infos}{ additional information. }
}
\references{
  D. Altman, D. Machin, T. Bryant, M. Gardner (eds). \emph{Statistics with Confidence:
  Confidence Intervals and Statistical Guidelines}, 2nd edition.
  John Wiley and Sons 2000.

  J. Hedderich, L. Sachs. \emph{Angewandte Statistik: Methodensammlung mit R}.
  Springer 2018.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
\examples{
x <- rnorm(100)
y <- rnorm(100, sd = 2)
## paired
normDiffCI(x, y, paired = TRUE)
## (R = 999 to reduce computation time for R checks)
normDiffCI(x, y, paired = TRUE, boot = TRUE, R = 999)
## compare
normCI(x-y)
## (R = 999 to reduce computation time for R checks)
normCI(x-y, boot = TRUE, R = 999)

## unpaired
y <- rnorm(90, mean = 1, sd = 2)
## classical
normDiffCI(x, y, method = "classical")
## (R = 999 to reduce computation time for R checks)
normDiffCI(x, y, method = "classical", boot = TRUE, R = 999)
## Welch (default as in case of function t.test)
normDiffCI(x, y, method = "welch")
## (R = 999 to reduce computation time for R checks)
normDiffCI(x, y, method = "welch", boot = TRUE, R = 999)
## Hsu
normDiffCI(x, y, method = "hsu")
## In case of bootstrap there is no difference between welch and hsu
## (R = 999 to reduce computation time for R checks)
normDiffCI(x, y, method = "hsu", boot = TRUE, R = 999)

## one-sided
normDiffCI(x, y, alternative = "less")
normDiffCI(x, y, boot = TRUE, R = 999, alternative = "greater")

\donttest{
## parallel computing for bootstrap
normDiffCI(x, y, method = "welch", boot = TRUE, R = 9999, 
           parallel = "multicore", ncpus = 2)
}

\donttest{
## Monte-Carlo simulation: coverage probability
M <- 100 # increase for more stable/realistic results!
CIhsu <- CIwelch <- CIclass <- matrix(NA, nrow = M, ncol = 2)
for(i in 1:M){
  x <- rnorm(10)
  y <- rnorm(30, sd = 0.1)
  CIclass[i,] <- normDiffCI(x, y, method = "classical")$conf.int
  CIwelch[i,] <- normDiffCI(x, y, method = "welch")$conf.int
  CIhsu[i,] <- normDiffCI(x, y, method = "hsu")$conf.int
}
## coverage probabilies
## classical
sum(CIclass[,1] < 0 & 0 < CIclass[,2])/M
## Welch
sum(CIwelch[,1] < 0 & 0 < CIwelch[,2])/M
## Hsu
sum(CIhsu[,1] < 0 & 0 < CIhsu[,2])/M
}
}
\keyword{univar}
