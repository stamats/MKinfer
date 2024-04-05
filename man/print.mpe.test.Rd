\name{print.mpe.test}
\alias{print.mpe.test}
\title{Print Method for MPE Hypothesis Test}
\description{
  Printing objects of class \code{"mpe.test"} by simple \code{\link{print}} methods.
}
\usage{
  \method{print}{mpe.test}(x, digits = getOption("digits"), prefix = "\t", ...)
}
\arguments{
  \item{x}{object of class \code{"mpe.test"}.}
  \item{digits}{number of significant digits to be used.}
  \item{prefix}{string, passed to \code{\link{strwrap}} for displaying
    the \code{method} component of the \code{mpe.test} object.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  The \code{print} method is based on the respective method \code{print.htest} 
  of package \pkg{stats}.
}
\value{
  the argument \code{x}, invisibly, as for all \code{\link{print}} methods.
}
\author{Srinath Kolampally, Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{The function first appeared in package \pkg{mpe}, which is now archived on CRAN.}
\seealso{\code{\link[stats]{print.power.htest}}, \code{\link{mpe.z.test}},
\code{\link{mpe.t.test}}.
}
\examples{
library(mvtnorm)
delta <- c(0.25, 0.5)
Sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
n <- 50
X <- rmvnorm(n=n, mean = delta, sigma = Sigma)
Y <- rmvnorm(n=n, mean = rep(0, length(delta)), sigma = Sigma)
mpe.t.test(X = X, Y = Y)
}
\keyword{htest}
\keyword{power.htest}
