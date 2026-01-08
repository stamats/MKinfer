\name{xiao.t.test}
\alias{xiao.t.test}
\alias{xiao.t.test.default}
\alias{xiao.t.test.formula}
\title{Xiao Two-Sample t-Test}
\description{
  Performs Xiao two sample t-tests on vectors of data.
}
\usage{
xiao.t.test(x, \dots)

\method{xiao.t.test}{default}(x, y,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, conf.level = 0.95, \dots)

\method{xiao.t.test}{formula}(formula, data, subset, na.action, \dots)
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
  The function and its documentation was adapted from \code{\link{hsu.t.test}}.

  \code{alternative = "greater"} is the alternative that \code{x} has a
  larger mean than \code{y}.

  If the input data are effectively constant (compared to the larger of the
  two means) an error is generated.

  The test is based on the solution of a generalized Behrens-Fisher problem 
  derived by Xiao (2018) and corresponds to the two-sample test derived in 
  Section 5.3 of Xiao (2018). The test is based on a generalized t-distribution.
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
  Y. Xiao (2018). On the Solution of a Generalized Behrens-Fisher Problem. 
  \emph{Far East Journal of Theoretical Statistics}, \bold{54} (1), 21-140.
}
\seealso{\code{\link{t.test}}, \code{\link{hsu.t.test}}}
\examples{
## Examples taken and adapted from function t.test
t.test(1:10, y = c(7:20))      # P = .00001855
t.test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
xiao.t.test(1:10, y = c(7:20))
xiao.t.test(1:10, y = c(7:20, 200))

## Traditional interface
with(mtcars, t.test(mpg[am == 0], mpg[am == 1]))
with(mtcars, xiao.t.test(mpg[am == 0], mpg[am == 1]))
## Formula interface
t.test(mpg ~ am, data = mtcars)
xiao.t.test(mpg ~ am, data = mtcars)

## Example 6.1 in Xiao (2018)
x <- c(134, 146, 104, 119, 124, 161, 107, 83, 113, 129, 97, 123)
y <- c(70, 118, 101, 85, 107, 132, 94)
xiao.t.test(x, y, alternative = "greater")
t.test(x, y, var.equal = TRUE, alternative = "greater")
t.test(x, y, alternative = "greater")
hsu.t.test(x, y, alternative = "greater")
}
\keyword{htest}
