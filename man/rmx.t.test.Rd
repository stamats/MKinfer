\name{rmx.t.test}
\alias{rmx.t.test}
\alias{rmx.t.test.default}
\alias{rmx.t.test.formula}
\title{RMX based t-Test}
\description{
  Performs one and two sample t-tests on vectors of data using optimally robust
  RMX estimators.
}
\usage{
rmx.t.test(x, \dots)

\method{rmx.t.test}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, conf.level = 0.95, 
                             method = c("welch", "xiao", "hsu", "student"), 
                             eps.lower = 0, eps.upper = 0.05, k = 3L, \dots)

\method{rmx.t.test}{formula}(formula, data, subset, na.action = na.pass, \dots)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{y}{an optional (non-empty) numeric vector of data values.}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{mu}{a number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).}
  \item{paired}{a logical indicating whether you want a paired
    t-test.}
  \item{conf.level}{confidence level of the interval.}
  \item{method}{ method that should be used to estimate the variance and the 
    degrees of freedom in the two-sample case; see also \code{\link[stats]{t.test}},
    \code{\link{hsu.t.test}} and \code{\link{xiao.t.test}}.}
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
    lower bound for the amount of gross errors; see \code{\link[robrmx]{rmx}}. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
    upper bound for the amount of gross errors; see \code{\link[robrmx]{rmx}}. }
  \item{k}{ positive integer: k-step is used to compute the optimally robust estimator;
    see \code{\link[robrmx]{rmx}}. }
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs}
    is a numeric variable giving the data values and \code{rhs} either 
    \code{1} for a one-sample or paired test or a factor
    with two levels giving the corresponding groups. If \code{lhs} is of 
    class \code{"\link{Pair}"} and \code{rhs} is \code{1}, a paired test
    is done, see Examples.}
  \item{data}{an optional matrix or data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{\link{NA}}s.}
  \item{\dots}{further arguments to be passed to or from methods.
    For the \code{formula} method, this includes arguments of the
    default method, but not \code{paired}.}
}
\details{
  \code{alternative = "greater"} is the alternative that \code{x} has a
  larger mean than \code{y}. For the one-sample case: that the mean is positive.

  If \code{paired} is \code{TRUE} then both \code{x} and \code{y} must
  be specified and they must be the same length.  Missing values are
  silently removed (in pairs if \code{paired} is \code{TRUE}).  

  If the input data are effectively constant (compared to the larger of the
  two means) an error is generated.
  
  The function uses optimally robust RMX estimators to estimate mean and
  variance; 
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the t-statistic.}
  \item{parameter}{the degrees of freedom for the t-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the mean appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{the estimated mean or difference in means depending on
    whether it was a one-sample test or a two-sample test.}
  \item{null.value}{the specified hypothesized value of the mean or mean
    difference depending on whether it was a one-sample test or a
    two-sample test.}
  \item{stderr}{the standard error of the mean (difference), used as
    denominator in the t-statistic formula.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string indicating what type of t-test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\seealso{ \code{\link[stats]{t.test}} \code{\link{hsu.t.test}}, \code{\link{xiao.t.test}} }
\examples{
## Two-sample t-test
t.test(1:10, y = c(7:20))      # P = .00001855
t.test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
## Welch df
rmx.t.test(1:10, y = c(7:20))      # P = .003322
rmx.t.test(1:10, y = c(7:20, 200)) # P = .003838  -- STILL significant
## Hsu df
rmx.t.test(1:10, y = c(7:20), method = "hsu")      # P = .00006961
rmx.t.test(1:10, y = c(7:20, 200), method = "hsu") # P = .0009753  -- STILL significant
## Xiao generalized t
rmx.t.test(1:10, y = c(7:20), method = "xiao")      # P = .00004425
rmx.t.test(1:10, y = c(7:20, 200), method = "xiao") # P = .00007251  -- STILL significant


## Traditional interface
with(mtcars, t.test(mpg[am == 0], mpg[am == 1]))
with(mtcars, rmx.t.test(mpg[am == 0], mpg[am == 1]))

## Formula interface
t.test(mpg ~ am, data = mtcars)
rmx.t.test(mpg ~ am, data = mtcars)

## One-sample t-test
## Traditional interface
t.test(sleep$extra)
rmx.t.test(sleep$extra)

## Formula interface
t.test(extra ~ 1, data = sleep)
rmx.t.test(extra ~ 1, data = sleep)

## Paired t-test
## The sleep data is actually paired, so could have been in wide format:
sleep2 <- reshape(sleep, direction = "wide",
                  idvar = "ID", timevar = "group")

## Traditional interface
t.test(sleep2$extra.1, sleep2$extra.2, paired = TRUE)
rmx.t.test(sleep2$extra.1, sleep2$extra.2, paired = TRUE)

## Formula interface
t.test(Pair(extra.1, extra.2) ~ 1, data = sleep2)
rmx.t.test(Pair(extra.1, extra.2) ~ 1, data = sleep2)
}
\keyword{htest}
