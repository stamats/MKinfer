\name{hsu.t.test}
\alias{hsu.t.test}
\alias{hsu.t.test.default}
\alias{hsu.t.test.formula}
\title{Hsu Two-Sample t-Test}
\description{
  Performs Hsu two sample t-tests on vectors of data.
}
\usage{
hsu.t.test(x, \dots)

\method{hsu.t.test}{default}(x, y,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, conf.level = 0.95, \dots)

\method{hsu.t.test}{formula}(formula, data, subset, na.action, \dots)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{y}{a (non-empty) numeric vector of data values.}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{mu}{a number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).}
  \item{conf.level}{confidence level of the interval.}
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs}
    is a numeric variable giving the data values and \code{rhs} a factor
    with two levels giving the corresponding groups.}
  \item{data}{an optional matrix or data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  The function and its documentation was adapted from \code{\link{t.test}}.

  \code{alternative = "greater"} is the alternative that \code{x} has a
  larger mean than \code{y}.

  If the input data are effectively constant (compared to the larger of the
  two means) an error is generated.

  Following the recommendation of Section 6.8.3 and 7.4.4.2 of Hedderich and 
  Sachs (2018), one should at least have six observations per group to apply 
  the test.
  
  The test is conservative, but is a valid alternative when the sample sizes
  are very small in one group and larger in the second group as considered 
  in Adusah and Brooks (2011).
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the t-statistic.}
  \item{parameter}{the degrees of freedom for the t-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the mean appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{the estimated means and standard deviations.}
  \item{null.value}{the specified hypothesized value of the mean or mean
    difference depending on whether it was a one-sample test or a
    two-sample test.}
  \item{stderr}{the standard error of the difference in means, used as
    denominator in the t-statistic formula.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string indicating what type of t-test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
  A.K. Adusah, G.P. Brooks (2011). Type I Error Inflation of the Separate-Variances 
  Welch t test with Very Small Sample Sizes when Assumptions Are Met. 
  \emph{Journal of Modern Applied Statistical Methods} \bold{10}, 
  362-372.

  J. Hedderich, L. Sachs. \emph{Angewandte Statistik: Methodensammlung mit R}.
  Springer 2018.
  
  P. Hsu (1938). Contribution to the theory of ``student's'' t-test as 
  applied to the problem of two samples. \emph{Statistical Research Memoirs} \bold{2}, 
  1-24.
}
\seealso{\code{\link{t.test}}}
\examples{
## Examples taken and adapted from function t.test
t.test(1:10, y = c(7:20))      # P = .00001855
t.test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
hsu.t.test(1:10, y = c(7:20))
hsu.t.test(1:10, y = c(7:20, 200))

## Traditional interface
with(mtcars, t.test(mpg[am == 0], mpg[am == 1]))
with(mtcars, hsu.t.test(mpg[am == 0], mpg[am == 1]))
## Formula interface
t.test(mpg ~ am, data = mtcars)
hsu.t.test(mpg ~ am, data = mtcars)
}
\keyword{htest}
