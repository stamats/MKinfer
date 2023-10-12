\name{mi.wilcox.test}
\alias{mi.wilcox.test}
\alias{mi.wilcox.test.default}
\alias{mi.wilcox.test.amelia}
\alias{mi.wilcox.test.mids}
\title{Multiple Imputation Wilcoxon Rank Sum and Signed Rank Tests}
\description{
  Performs one and two sample Wilcoxon tests on multiple imputed datasets.
}
\usage{
mi.wilcox.test(miData, \dots)

\method{mi.wilcox.test}{default}(miData, x, y = NULL,
        alternative = c("two.sided", "less", "greater"), mu = 0,
        paired = FALSE, exact = NULL, conf.int = TRUE,
        conf.level = 0.95, subset = NULL, \dots)

\method{mi.wilcox.test}{amelia}(miData, x, y = NULL,
        alternative = c("two.sided", "less", "greater"), mu = 0,
        paired = FALSE, exact = NULL, conf.int = TRUE,
        conf.level = 0.95, subset = NULL, \dots)

\method{mi.wilcox.test}{mids}(miData, x, y = NULL,
        alternative = c("two.sided", "less", "greater"), mu = 0,
        paired = FALSE, exact = NULL, conf.int = TRUE,
        conf.level = 0.95, subset = NULL, \dots)
}
\arguments{
  \item{miData}{list of multiple imputed datasets.}
  \item{x}{name of a variable that shall be tested.}
  \item{y}{an optional name of a variable that shall be tested (paired test)
  or a variable that shall be used to split into groups (unpaired test).}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{mu}{a number indicating the true value of the mean (or
    difference in means if you are performing a two sample test).}
  \item{paired}{a logical indicating whether you want a paired
    t-test.}
  \item{exact}{a logical indicating whether an exact p-value should be computed.}
  \item{conf.int}{a logical indicating whether a confidence interval should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  For details about the tests see \code{\link[exactRankTests]{wilcox.exact}}

  We use the median p rule (MPR) for the computation of the p value of the test; see 
  Section 5.3.2 of van Buuren (2018) or Section 13.3 in Heymans and Eekhout (2019).
  The approach seems to work well in many situations such as logistic regression
  (Eekhout et al. (2017)) or GAM (Bolt et al. (2022)). However, we are not aware 
  of any work that has investigated the MPR approach for Wilcoxon tests. Hence,
  this function should be regarded as experimental.
  
  We recommend to use an odd number of imputations.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{he value of the test statistic with a name describing it.}
  \item{p.value}{the p-value for the test.}
  \item{pointprob}{this gives the probability of observing the test statistic itself 
    (called point-prob).}
  \item{null.value}{the location parameter mu.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{conf.int}{a confidence interval for the location parameter. 
    (Only present if argument conf.int = TRUE.)}
  \item{estimate}{Hodges-Lehmann estimate of the location parameter. 
    (Only present if argument conf.int = TRUE.)}
}
\references{
  van Buuren, S. (2018). \emph{Flexible Imputation of Missing Data}.
  Chapman & Hall/CRC. \url{https://stefvanbuuren.name/fimd/}.

  Heymans, M.W. and Eekhout, I. (2019). \emph{Applied Missing Data Analysis With SPSS 
  and (R)Studio}. Self-publishing. \url{https://bookdown.org/mwheymans/bookmi/}.
  
  Eekhout, I, van de Wiel, MA, Heymans, MW (2017). Methods for significance testing 
  of categorical covariates in logistic regression models after multiple imputation: 
  power and applicability analysis. \emph{BMC Med Res Methodol}, \bold{17}, 1:129.
  \doi{10.1186/s12874-017-0404-7}
  
  Bolt, MA, MaWhinney, S, Pattee, JW, Erlandson, KM, Badesch, DB, Peterson, RA (2022). 
  Inference following multiple imputation for generalized additive models: an 
  investigation of the median p-value rule with applications to the Pulmonary 
  Hypertension Association Registry and Colorado COVID-19 hospitalization data. 
  \emph{BMC Med Res Methodol}, \bold{22}, 1:148. \doi{10.1186/s12874-022-01613-w}.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{
  \code{\link[exactRankTests]{wilcox.exact}}
}
\examples{
## Generate some data
set.seed(123)
x <- rnorm(25, mean = 1)
x[sample(1:25, 5)] <- NA
y <- rnorm(20, mean = -1)
y[sample(1:20, 4)] <- NA
pair <- c(rnorm(25, mean = 1), rnorm(20, mean = -1))
g <- factor(c(rep("yes", 25), rep("no", 20)))
D <- data.frame(ID = 1:45, response = c(x, y), pair = pair, group = g)

## Use Amelia to impute missing values
library(Amelia)
res <- amelia(D, m = 9, p2s = 0, idvars = "ID", noms = "group")

## Per protocol analysis (Exact Wilcoxon rank sum test)
library(exactRankTests)
wilcox.exact(response ~ group, data = D, conf.int = TRUE)
## Intention to treat analysis (Multiple Imputation Exact Wilcoxon rank sum test)
mi.wilcox.test(res, x = "response", y = "group")

## Specifying alternatives
mi.wilcox.test(res, x = "response", y = "group", alternative = "less")
mi.wilcox.test(res, x = "response", y = "group", alternative = "greater")

## One sample test
wilcox.exact(D$response[D$group == "yes"], conf.int = TRUE)
mi.wilcox.test(res, x = "response", subset = D$group == "yes")
mi.wilcox.test(res, x = "response", mu = -1, subset = D$group == "yes",
               alternative = "less")
mi.wilcox.test(res, x = "response", mu = -1, subset = D$group == "yes",
               alternative = "greater")

## paired test
wilcox.exact(D$response, D$pair, paired = TRUE, conf.int = TRUE)
mi.wilcox.test(res, x = "response", y = "pair", paired = TRUE)

## Use mice to impute missing values
library(mice)
res.mice <- mice(D, m = 9, print = FALSE)
mi.wilcox.test(res.mice, x = "response", y = "group")
}
\keyword{htest}
