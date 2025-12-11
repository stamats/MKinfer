\name{mdplot}
\alias{mdplot}
\title{ Visualize MD. }
\description{
  The function visualizes the mean difference (MD) of two groups, assuming two 
  normal distributions. In addition, the sensitivity and specificity based on the 
  optimal cutoff (Youden's J statistic) can be added.
}
\usage{
mdplot(delta, sd1 = 1, sd2 = 1, constant = 3, n = 501, alpha = 0.25, 
        lwd = 1, col1 = "darkblue", col2 = "darkred", xlim, 
        addSens = TRUE, digits = 3)
}
\arguments{
  \item{delta}{ (absolute) mean difference of the two groups. }
  \item{sd1}{ standard deviation of first group. }
  \item{sd2}{ standard deviation of second group. }
  \item{constant}{ used to define the plot range (\code{constant*sd}). }
  \item{n}{integer; the number of values at which the two normal densities are evaluated. }
  \item{alpha}{ opacity used for \code{geom_area}. }
  \item{lwd}{ linewidth of the plotted lines. }
  \item{col1}{ color of density of first group. }
  \item{col2}{ color of density of second group. }
  \item{xlim}{ limits of x-axis. If missing, calculated based on \code{constant} 
  and the specified SDs. }
  \item{addSens}{ logical; should sensitivity and specificity be added to the plot. }
  \item{digits}{ number of digits used for sensitivity and specificity. }
}
\details{
  The function can be used to demonstrate (standardized) effect sizes and how 
  these effect sizes transform into sensitivities and specificities.
  
  The two standard diviations are assumed to be equal if their difference is
  smaller than \code{1e-08}.
}
\value{\code{ggplot} object}
\references{
  B\"ohning D, B\"ohning W, Holling H (2008). Revisiting Youden's index as a useful 
  measure of the misclassification error in meta-analysis of diagnostic studies. 
  Stat Methods Med Res. 17(6), 543-54.
  
  J.H. Zhang, T.D.Y. Chung and K.R. Oldenburg (1999). A simple statistical 
  parameter for use in evaluation and validation of high throughput screening assays. 
  \emph{Journal of Biomolecular Screening} \bold{4} (2): 67-73.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\examples{
## small effect size
mdplot(delta = 0.2)
## medium effect size
mdplot(delta = 0.5)
## large effect size
mdplot(delta = 0.8)
## z-factor = 0  (z-factor = 1 - 3*2*sd/delta); see Zhang et al. (1999)
mdplot(delta = 6)
## z-factor = 0.5  (z-factor = 1 - 3*2*sd/delta); see Zhang et al. (1999)
mdplot(delta = 12)

## unequal variances
mdplot(delta = 0.8, sd1 = 1, sd2 = 2)
mdplot(delta = 0.8, sd1 = 2, sd2 = 1)
}
\keyword{univar}
