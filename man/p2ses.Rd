\name{p2ses}
\alias{p2ses}
\title{ Compute SES from p value. }
\description{
  The function computes the SES (standardized effect size) from the p value 
  for permutation/randomisation tests as proposed by Botta-Dukat (2018).
}
\usage{
p2ses(p, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{p}{ numeric vector of p values. }
  \item{alternative}{ a character string specifying the alternative hypothesis, 
  must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. 
  You can specify just the initial letter.}
}
\details{
  The function uses the probit transformation (\code{qnorm}) to compute 
  an alternative SES based on p values from a permutation/randomization test 
  as proposed by Botta-Dukat (2018) for skewed distributions. 
}
\value{SES}
\references{
  Botta-Dukat, Z (2018). Cautionary note on calculating standardized effect size (SES) 
  in randomization test. Community Ecology 19(1), 77-83.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\examples{
## symmetric case
x <- rnorm(100)
res <- perm.t.test(x)
p2ses(res$p.value)
abs(res$statistic)

## skewed case
x <- rgamma(100, shape = 5)
res <- perm.t.test(x, mu = 5)
p2ses(res$p.value)
abs(res$statistic)
}
\keyword{univar}
