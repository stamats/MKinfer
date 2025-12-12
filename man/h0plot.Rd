\name{h0plot}
\alias{h0plot}
\alias{h0plot.default}
\alias{h0plot.perm.htest}
\alias{h0plot.boot.htest}
\title{ Visualize MD. }
\description{
  The function visualizes the distribution of the test statistic.
}
\usage{
h0plot(x, \dots)

\S3method{h0plot}{default}(x, statistic, pval, alternative = "two.sided", sig.level = 0.05, 
                           hist.alpha = 0.2, stat.col = "darkgreen", dens.col = "black",
                           dens.alpha = 0.5, cval.col = "darkred", lwd = 1, 
                           breaks = NULL, digits = 3, \dots)

\S3method{h0plot}{perm.htest}(x, sig.level = 0.05, hist.alpha = 0.2, stat.col = "darkgreen", 
                              dens.col = "black", dens.alpha = 0.5, 
                              cval.col = "darkred", lwd = 1, breaks = NULL, 
                              digits = 3, \dots)

\S3method{h0plot}{boot.htest}(x, sig.level = 0.05, hist.alpha = 0.2, stat.col = "darkgreen", 
                              dens.col = "black", dens.alpha = 0.5, cval.col = "darkred", 
                              lwd = 1, breaks = NULL, digits = 3, \dots)
}
\arguments{
  \item{x}{ bootstrap/permuation results of the test statistic or result of
  class \code{perm.htest} or \code{boot.htest}. }
  \item{statistic}{ test statistic of the test on the original data. }
  \item{pval}{ p value of the test. }
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{sig.level}{ significance level of the test. }
  \item{hist.alpha}{ opacity used for \code{geom_histogram}. }
  \item{stat.col}{ color used for line representing the (final) test statistic. }
  \item{dens.col}{ color used for the density line. }
  \item{dens.alpha}{ opacity used for \code{geom_area} to color the rejection
    region. }
  \item{cval.col}{ color used for rejection region (area under the density beyond 
    the critical value). }
  \item{lwd}{ line width used for line representing the (final) test statistic 
    as well ase for \code{geom_density}. }
  \item{breaks}{ argument \code{breaks} of \code{geom_histogram}. }
  \item{digits}{ number of significant digits used for p values. }
  \item{\dots}{ further arguments that may be passed through.}
}
\details{
  The function can be used to plot the result of a permutation or bootstrap 
  test.
}
\value{\code{ggplot} object}
\references{
  B. Efron, R.J. Tibshirani. \emph{An Introduction to the Bootstrap}.
  Chapman and Hall/CRC 1993.
  
  A. Janssen (1997). Studentized permutation tests for non-i.i.d, hypotheses 
  and the generalized Behrens-Fisher problem. 
  \emph{Statistics and Probability Letters}, \bold{36}, 9-21.
  
  E. Chung, J.P. Romano (2013). Exact and asymptotically robust permutation
  tests. \emph{The Annals of Statistics}, \bold{41}(2), 484-507.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\examples{
## Bootstrap based Welch t-Test
h0plot(boot.t.test(extra ~ group, data = sleep, bootStat = TRUE))
h0plot(boot.t.test(extra ~ group, data = sleep, bootStat = TRUE, 
                   alternative = "less"))

## Permutation based Welch t-Test
h0plot(perm.t.test(extra ~ group, data = sleep, permStat = TRUE))
h0plot(perm.t.test(extra ~ group, data = sleep, permStat = TRUE, 
                   alternative = "less"))
}
\keyword{hplot}
