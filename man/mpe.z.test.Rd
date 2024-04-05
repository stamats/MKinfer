\name{mpe.z.test}
\alias{mpe.z.test}
\title{Intersection-Union z-Test for Testing Multiple Co-Primary Endpoints}
\description{
The function computes the intersection-union z-test.
}
\usage{
mpe.z.test(X, Y, Sigma, conf.level = 0.975)
}
\arguments{
  \item{X}{matrix with observations of group 1 in rows}
  \item{Y}{matrix with obersvations of group 2 in rows}
  \item{Sigma}{known covariance matrix.}
  \item{conf.level}{confidence level of the interval.}
}
\details{
  The function computes the intersection-union z-test.
  The implementation is based on the formulas given in the references below.

  The null hypothesis reads \eqn{\mu_{Tk}-\mu_{Ck}\le 0}{mu_Tk-mu_Ck <= 0} for
  at least one \eqn{k\in\{1,\ldots,K\}}{k in {1,...,K}} where Tk is treatment k,
  Ck is control k and K is the number of co-primary endpoints (i.e. number of
  columns of \code{X} and \code{Y}).
}
\value{
Object of class \code{"mpe.test"}.
}
\references{
  Sugimoto, T. and Sozu, T. and Hamasaki, T. (2012). A convenient formula for sample
  size calculations in clinical trials with multiple co-primary continuous endpoints.
  \emph{Pharmaceut. Statist.}, \bold{11}: 118-128. doi:10.1002/pst.505

  Sozu, T. and Sugimoto, T. and Hamasaki, T. and Evans, S.R. (2015). \emph{Sample
  Size Determination in Clinical Trials with Multiple Endpoints}. Springer Briefs in
  Statistics, ISBN 978-3-319-22005-5.
}
\author{Srinath Kolampally, Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{The function first appeared in package \pkg{mpe}, which is now archived on CRAN.}
\seealso{\code{\link{mpe.t.test}}}
\examples{
library(mvtnorm)
delta <- c(0.25, 0.5)
Sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
n <- 50
X <- rmvnorm(n=n, mean = delta, sigma = Sigma)
Y <- rmvnorm(n=n, mean = rep(0, length(delta)), sigma = Sigma)
mpe.z.test(X = X, Y = Y, Sigma = Sigma)
}
\keyword{htest}
\keyword{multivariate}
