\name{rm.oneway.test}
\alias{rm.oneway.test}
\title{Test for Equal Means in a Repeated Measures One-Way Layout}
\description{
  Test whether two or more samples have the same locations in a repeated
  measures setting. 
}
\usage{
rm.oneway.test(x, g, id, method = "aov")
}
\arguments{
  \item{x}{ numeric, response (outcome, dependent) variable.}
  \item{g}{ factor, grouping (independent) variable.}
  \item{id}{ factor, subject id (blocking variable).}
  \item{method}{ name of method, possible methods are \code{"aov"}, \code{"lme"}, 
  \code{"friedman"}, \code{"quade"}}
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the degrees of freedom of the exact or approximate F
    distribution of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating the test performed.}
  \item{data.name}{a character string giving the names of the data.}
}
\details{
  The function wraps the functions \code{\link{aov}}, \code{\link[nlme]{lme}}, 
  \code{\link{friedman.test}} and \code{\link{quade.test}} into one function 
  for a repeated measures one-way layout.
}
\references{
  Chambers, J. M., Freeny, A and Heiberger, R. M. (1992), 
  \emph{Analysis of variance; designed experiments. Chapter 5 of Statistical Models in S},
  eds J. M. Chambers and T. J. Hastie, Wadsworth and Brooks/Cole.
  
  Pinheiro, J.C., and Bates, D.M. (2000), \emph{Mixed-Effects Models in S and S-PLUS}, 
  Springer.
  
  Myles Hollander and Douglas A. Wolfe (1973), \emph{Nonparametric Statistical Methods}. 
  New York: John Wiley and Sons. Pages 139-146.

  D. Quade (1979), Using weighted rankings in the analysis of complete blocks 
  with additive block effects. \emph{Journal of the American Statistical Association} 
  \bold{74}, 680-683.

  William J. Conover (1999), \emph{Practical nonparametric statistics}. 
  New York: John Wiley and Sons. Pages 373-380.
}
\seealso{
  \code{\link{aov}}, \code{\link[nlme]{lme}}, \code{\link{friedman.test}}, 
  \code{\link{quade.test}}
}
\examples{
set.seed(123)
outcome <- c(rnorm(10), rnorm(10, mean = 1.5), rnorm(10, mean = 1))
timepoints <- factor(rep(1:3, each = 10))
patients <- factor(rep(1:10, times = 3))
rm.oneway.test(outcome, timepoints, patients)
rm.oneway.test(outcome, timepoints, patients, method = "lme")
rm.oneway.test(outcome, timepoints, patients, method = "friedman")
rm.oneway.test(outcome, timepoints, patients, method = "quade")
}
\keyword{htest}
