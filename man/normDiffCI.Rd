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
           boot = FALSE, R = 1000, type = "all", na.rm = TRUE)
meanDiffCI(x, y, conf.level = 0.95, paired = FALSE, method = "welch", 
           boot = FALSE, R = 1000, type = "all", na.rm = TRUE)
}
\arguments{
  \item{x}{ numeric vector of data values of group 1. }
  \item{y}{ numeric vector of data values of group 2. }
  \item{conf.level}{ confidence level. }
  \item{paired}{ a logical value indicating whether the two groups are paired. }
  \item{method}{ a character string specifing which method to use in the unpaired case; see details. }
  \item{boot}{ a logical value indicating whether bootstrapped confidence intervals shall be computed. }
  \item{R}{ number of bootstrap replicates; see details. }
  \item{type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be stripped before the computation proceeds. }
}
\details{
  The standard confidence intervals for the difference of means are
  computed that can be found in many textbooks, e.g. Chapter 4 in Altman et al. (2000).

  The method \code{"classical"} assumes equal variances whereas methods
  \code{"welch"} and \code{"hsu"} allow for unequal variances. The latter two
  methods use different formulas for computing the degrees of freedom of the
  respective t-distribution providing the quantiles in the confidence interval.
  Instead of the Welch-Satterhwaite equation the method of Hsu uses the minimum
  of the group sample sizes minus 1; see Section 6.8.3 of Hedderich and Sachs (2016).
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
  Springer 2016.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
\examples{
x <- rnorm(100)
y <- rnorm(100, sd = 2)
## paired
normDiffCI(x, y, paired = TRUE)
normDiffCI(x, y, paired = TRUE, boot = TRUE)
## compare
normCI(x-y)
normCI(x-y, boot = TRUE)

## unpaired
y <- rnorm(90, mean = 1, sd = 2)
## classical
normDiffCI(x, y, method = "classical")
normDiffCI(x, y, method = "classical", boot = TRUE)
## Welch (default is in case of function t.test)
normDiffCI(x, y, method = "welch")
normDiffCI(x, y, method = "welch", boot = TRUE)
## Hsu
normDiffCI(x, y, method = "hsu")
## In case of bootstrap there is no difference between welch and hsu
normDiffCI(x, y, method = "hsu", boot = TRUE)

\donttest{
## Monte-Carlo simulation: coverage probability
M <- 100
CIbootbca <- CIhsu <- CIwelch <- CIclass <- matrix(NA, nrow = M, ncol = 2)
for(i in 1:M){
  x <- rnorm(10)
  y <- rnorm(30, sd = 0.1)
  CIclass[i,] <- normDiffCI(x, y, method = "classical")$conf.int
  CIwelch[i,] <- normDiffCI(x, y, method = "welch")$conf.int
  CIhsu[i,] <- normDiffCI(x, y, method = "hsu")$conf.int
  CIbootbca[i,] <- normDiffCI(x, y, method = "welch", boot = TRUE, R = 1000, 
                              type = "bca")$conf.int$bca[c(4,5)]
}
## coverage probabilies
## classical
sum(CIclass[,1] < 0 & 0 < CIclass[,2])/M
## Welch
sum(CIwelch[,1] < 0 & 0 < CIwelch[,2])/M
## Hsu
sum(CIhsu[,1] < 0 & 0 < CIhsu[,2])/M
## Bootstrap bca
sum(CIbootbca[,1] < 0 & 0 < CIbootbca[,2])/M
}
}
\keyword{univar}
