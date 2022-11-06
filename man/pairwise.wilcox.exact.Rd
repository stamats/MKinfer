\name{pairwise.wilcox.exact}
\alias{pairwise.wilcox.exact}
\title{ Compute Pairwise Wilcoxon Tests }
\description{
  The function computes pairwise Wilcoxon rank sum and signed rank tests
  using function \code{wilcox.exact} of package exactRankTests.
}
\usage{
pairwise.wilcox.exact(x, g, p.adjust.method = "holm", paired = FALSE, ...)
}
\arguments{
  \item{x}{ numeric vector. }
  \item{g}{ grouping vector or factor }
  \item{p.adjust.method}{ method for adjusting p values (see \code{\link[stats]{p.adjust}}). 
                          Can be abbreviated. }
  \item{paired}{ a logical indicating whether you want a paired test.}
  \item{\dots}{ additional arguments to fun. }
}
\details{
  The function computes pairwise Wilcoxon rank sum and signed rank tests.

  The implementation is in certain aspects analogously to 
  \code{\link[stats]{pairwise.wilcox.test}}. However, a more detailed 
  output is generated.
}
\value{
  Object of class \code{"pw.htest"} containing the following components:
  \item{data.name}{ a character string giving the names of the data.}
  \item{method}{ the type of test applied.}
  \item{null.value}{ the location parameter mu.}
  \item{alternative}{ a character string describing the alternative hypothesis.}
  \item{conf.level}{ confidence level of the confidence interval.}
  \item{results}{ a data.frame containing the results of function 
                  \code{\link[exactRankTests]{wilcox.exact}}}
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{ \code{\link[exactRankTests]{wilcox.exact}}, 
          \code{\link[stats]{pairwise.wilcox.test}} }
\examples{
set.seed(13)
x <- rnorm(100)
g <- factor(sample(1:4, 100, replace = TRUE))
levels(g) <- c("a", "b", "c", "d")
pairwise.wilcox.exact(x, g)
## in contrast to
pairwise.wilcox.test(x, g)
}
\keyword{univar}
