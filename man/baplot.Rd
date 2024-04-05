\name{baplot}
\alias{baplot}
\title{ Bland-Altman Plots }
\description{
  Produce classical parametric as well as nonparametric Bland-Altman plots 
  including confidence intervals for the mean of the differences and 
  the lower and upper limit of agreement.
}
\usage{
baplot(x, y, loa.level = 0.95, conf.level = 0.95, 
       xlab = "Mean", ylab = "Difference", 
       title = "Bland-Altman Plot", xlim = NULL, ylim = NULL,
       type = c("parametric", "nonparametric"), 
       loa.type = c("unbiased", "biased"), ci.diff = TRUE, ci.loa = TRUE, 
       ci.type = c("exact", "approximate", "boot"), 
       bootci.type = NULL, R = 9999, print.res = TRUE,
       color.low = "#4575B4", color.upp = "#D73027", 
       alpha = 1, shape = 19, na.rm = TRUE, ...)
}
\arguments{
  \item{x}{ numeric, measurements of method 1.}
  \item{y}{ numeric, measurements of method 2.}
  \item{loa.level}{ single numeric, level of the level of agreement.}
  \item{conf.level}{ single numeric, significance level.}
  \item{xlab}{ label of x-axis.}
  \item{ylab}{ label of y-axis.}
  \item{title}{ title of plot.}
  \item{xlim}{ optional numeric vector of length 2, limits of x-axis.}
  \item{ylim}{ optional numeric vector of length 2, limits of y-axis.}
  \item{type}{ single character, either \code{"parametric"} (default) or \code{"nonparametric"}; see also Details below.}
  \item{loa.type}{ single character, either \code{"unbiased"} (default) or \code{"biased"}; see also Details below.}
  \item{ci.diff}{ single logical, plot confidence interval for the mean of the differences.}
  \item{ci.loa}{ single logical, plot confidence intervals for the lower and upper limit of agreement.}
  \item{ci.type}{ single character, either \code{"exact"} (default) or \code{"approximate"} or \code{"boot"}; see also Details below.}
  \item{bootci.type}{ single character, type of bootstrap interval; see \code{\link[boot]{boot.ci}} and also Details below. }
  \item{R}{ single numeric, number of bootstrap replicates. }
  \item{print.res}{ single logical, print results of computations in addition to plotting.}
  \item{color.low}{ single color (character), color for lower limit of agreement.}
  \item{color.upp}{ single color (character), color for upper limit of agreement.}
  \item{alpha}{ blending factor (default: no blending.}
  \item{shape}{ point shape used.}
  \item{na.rm}{ single logical, remove \code{NA} values before plotting.}
  \item{\dots}{ further arguments passed to function \code{\link[boot]{boot}},
  e.g. for parallel computing. }
}
\details{
  The plot generates a \code{ggplot2} object that is shown.
  
  Setting \code{type = "parametric"} (default), a classical Bland-Altman plot 
  with mean of the differences and standard lower and upper limit of agreement 
  (mean +/- 1.96 SD) will be generated (\code{loa.level = 0.95}); 
  see Altman and Bland (1983,1986).
  
  Setting \code{type = "nonparametric"}, a nonparametric Bland-Altman plot with 
  median of the differences and empirical 2.5\% and 97.5\% quantiles 
  (\code{loa.level = 0.95}) as lower and upper limit of agreement will be generated.
  
  By changing \code{loa.level} the lower and upper limit of agreement will 
  correspond to the \code{(1-loa.level)/2} and the \code{1 - (1-low.level)/2}
  quantile.
  
  
  Setting \code{loa.type = "unbiased"}, the unbiased estimator of the standard
  deviation will be used; see Shieh (2018).
  
  Setting \code{loa.type = "biased"}, the standard deviation will be estimated
  by function \code{sd}.
  
  
  Setting \code{ci.type = "exact"}, the exact confidence intervals of Shieh (2018) 
  will be used for the lower and upper limit of agreement.
  
  Setting \code{ci.type = "approximate"}, the approximate confidence intervals
  of Bland and Altman (1986, 1999) will be used.
  
  Setting \code{ci.type = "boot"}, bootstrap confidence intervals will be used.
  In the case of parametric Bland-Altman plots, the studentized bootstrap method will
  be used by default. In the case of nonparametric Bland-Altman plots, the 
  bootstrap percentile method will be used by default.
  
  The argument \code{bootci.type} can be set to \code{"stud"} (studentized bootstrap),
  \code{"perc"} (bootstrap percentile) or \code{"bca"} (adjusted bootstrap percentile).
}
\value{
  Object of class \code{gg} and \code{ggplot}.
}
\references{
  Altman DG, Bland JM (1983). Measurement in Medicine: the Analysis of Method Comparison Studies. 
  \emph{The Statistician} 32:307-317.
  
  Bland JM, Altman DG (1986). Statistical methods for assessing agreement between two methods 
  of clinical measurement. \emph{Lancet} 327(8476):307-10

  Bland JM, Altman DG (1999). Measuring agreement in method comparison studies. 
  \emph{Stat Methods Med Res} 8:135-160.
  
  Shieh G (2018). The appropriateness of Bland-Altman's approximate confidence intervals 
  for limits of agreement. \emph{BMC Med Res Methodol} 18:45.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
\examples{
data(fingsys)
## classical parametric Bland-Altman plot with exact confidence intervals
baplot(x = fingsys$armsys, y = fingsys$fingsys, 
       xlab = "Mean of arm and finger [mm Hg]",
       ylab = "Difference (arm - finger) [mm Hg]")

## nonparametric Bland-Altman plot with exact confidence intervals
baplot(x = fingsys$armsys, y = fingsys$fingsys, type = "nonparametric",
       xlab = "Mean of arm and finger [mm Hg]",
       ylab = "Difference (arm - finger) [mm Hg]",
       title = "Nonparametric Bland-Altman Plot")
}
\keyword{hplot}
