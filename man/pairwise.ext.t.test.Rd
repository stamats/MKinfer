\name{pairwise.ext.t.test}
\alias{pairwise.ext.t.test}
\title{ Compute Pairwise t Tests }
\description{
  The function computes pairwise t tests using functions \code{t.test}, 
  \code{boot.t.test} or \code{perm.t.test}.
}
\usage{
pairwise.ext.t.test(x, g, method = "t.test", p.adjust.method = "holm", 
                    paired = FALSE, ...)
}
\arguments{
  \item{x}{ numeric vector. }
  \item{g}{ grouping vector or factor }
  \item{method}{ character giving the name of the function to be applied; that is, 
                 \code{"t.test"}, \code{"boot.t.test"} or \code{"perm.t.test"}. }
  \item{p.adjust.method}{ method for adjusting p values (see \code{\link[stats]{p.adjust}}). 
                          Can be abbreviated. }
  \item{paired}{ a logical indicating whether you want a paired test.}
  \item{\dots}{ additional arguments to fun. }
}
\details{
  The function computes pairwise t tests using function \code{\link[stats]{t.test}},
  \code{\link{boot.t.test}} or \code{\link{perm.t.test}}.

  The implementation is in certain aspects analogously to 
  \code{\link[stats]{pairwise.t.test}}. However, a more detailed 
  output is generated.
}
\value{
  Object of class \code{"pw.htest"} containing the following components:
  \item{data.name}{ a character string giving the names of the data. }
  \item{method}{ the type of test applied. }
  \item{null.value}{ the location parameter mu. }
  \item{alternative}{ a character string describing the alternative hypothesis. }
  \item{conf.level}{ confidence level of the confidence interval. }
  \item{results}{ a data.frame containing the results of function 
                  \code{\link[stats]{t.test}}, \code{\link{boot.t.test}} 
                  or \code{\link{perm.t.test}}. }
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{ \code{\link[stats]{t.test}}, \code{\link[stats]{pairwise.t.test}} }
\examples{
set.seed(13)
x <- rnorm(100)
g <- factor(sample(1:4, 100, replace = TRUE))
levels(g) <- c("a", "b", "c", "d")
pairwise.ext.t.test(x, g)
## in contrast to
pairwise.t.test(x, g)
## moreover
pairwise.ext.t.test(x, g, method = "boot.t.test")
pairwise.ext.t.test(x, g, method = "perm.t.test")
}
\keyword{univar}
