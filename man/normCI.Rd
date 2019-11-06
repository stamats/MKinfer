\name{normCI}
\alias{normCI}
\title{ Confidence Intervals for Mean and Standard Deviation }
\description{
  This function can be used to compute confidence intervals for mean and
  standard deviation of a normal distribution.
}
\usage{
normCI(x, mean = NULL, sd = NULL, conf.level = 0.95, 
       boot = FALSE, R = 1000, type = "all", na.rm = TRUE)
}
\arguments{
  \item{x}{ vector of observations. }
  \item{mean}{ mean if known otherwise \code{NULL}. }
  \item{sd}{ standard deviation if known otherwise \code{NULL}. }
  \item{conf.level}{ confidence level. }
  \item{boot}{ a logical value indicating whether bootstrapped confidence intervals shall be computed. }
  \item{R}{ number of bootstrap replicates; see details. }
  \item{type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
  \item{na.rm}{ a logical value indicating whether NA values should be stripped before the computation proceeds. }
}
\details{
  The standard confidence intervals for mean and standard deviation are
  computed that can be found in many textbooks, e.g. Chapter 4 in Altman et al. (2000).
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
## bootstrap interval
normCI(x, boot = TRUE)

## sd known
normCI(x, sd = 1)
## bootstrap interval only for mean (sd is ignored)
normCI(x, sd = 1, boot = TRUE)

## mean known
normCI(x, mean = 0)
## bootstrap interval only for sd (mean is ignored)
normCI(x, mean = 0, boot = TRUE)
}
\keyword{univar}
