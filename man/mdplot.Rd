\name{mdplot}
\alias{mdplot}
\alias{mdplot2}
\alias{md2sens}
\title{ Visualize MD. }
\description{
  The functions visualize the mean difference (MD) of two groups, assuming two 
  normal distributions, either with identical standard deviation (SD) or two 
  (different) SDs. In addition, the sensitivity and specificity based on the 
  optimal cutoff (Youden's J statistic) are computed.
}
\usage{
mdplot(delta, sd = 1, constant = 3, n = 501, alpha = 0.25, lwd = 1, 
       col1 = "darkblue", col2 = "darkred", xlim, addSens = TRUE, 
       digits = 3)
mdplot2(delta, sd1 = 1, sd2 = 1, constant = 3, n = 501, alpha = 0.25, 
        lwd = 1, col1 = "darkblue", col2 = "darkred", xlim, 
        addSens = TRUE, digits = 3)
md2sens(delta, sd = NULL, sd1 = NULL, sd2 = NULL, var.equal = TRUE)
}
\arguments{
  \item{delta}{ (absolute) mean difference of the two groups. }
  \item{sd}{ common standard deviation of both groups. }
  \item{sd1}{ standard deviation of first group. }
  \item{sd2}{ standard deviation of second group. }
  \item{constant}{ used to define the plot range (\code{constant*sd}). }
  \item{n}{integer; the number of values at which the two normal densities are evaluated. }
  \item{alpha}{ opacity used for \code{geom_area}. }
  \item{lwd}{ linewidth of the plotted lines. }
  \item{col1}{ color of density of first group. }
  \item{col2}{ color of density of second group. }
  \item{xlim}{ limits of x-axis. If missing, calculated based on \code{constant} 
  and the specified SD(s). }
  \item{addSens}{ logical; should sensitivity and specificity be added to the plot. }
  \item{digits}{ number of digits used for sensitivity and specificity. }
  \item{var.equal}{ logical; assume equal variance. If \code{TRUE} the value 
  of \code{sd} will be used for both groups, otherwise \code{sd1} and \code{sd2}
  have to be specified. }
}
\details{
  The functions can be used to demonstrate (standardized) effect sizes and how 
  these effect sizes transform into sensitivities and specificities.
}
\value{Functions \code{mdplot} and \code{mdplot2} return an ggplot object. 
Function \code{md2sens} returns a numeric vector with sensitivity and 
specificity.
}
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
mdplot(delta = 0.2, sd = 1)
## medium effect size
mdplot(delta = 0.5, sd = 1)
## large effect size
mdplot(delta = 0.8, sd = 1)
## z-factor = 0  (z-factor = 1 - 3*2*sd/delta); see Zhang et al. (1999)
mdplot(delta = 6, sd = 1)
## z-factor = 0.5  (z-factor = 1 - 3*2*sd/delta); see Zhang et al. (1999)
mdplot(delta = 12, sd = 1)

## unequal variances
mdplot2(delta = 0.8, sd1 = 1, sd2 = 2)
mdplot2(delta = 0.8, sd1 = 2, sd2 = 1)

## (standardized) mean difference to sensitivity/specificity
## equal variances
library(ggplot2)
delta <- seq(from = 0.0, to = 6, by = 0.05)
res <- sapply(delta, md2sens, sd = 1)
DF <- data.frame(SMD = delta, sensitivity = res[1,], 
                 specificity = res[2,])
ggplot(DF, aes(x = SMD, y = sensitivity)) +
  geom_line() + ylim(0.5, 1.0) + xlab("(standardized) mean difference") +
  ylab("sensitivity = specificity")

## unequal variances
delta <- seq(from = 0.0, to = 6, by = 0.05)
res <- sapply(delta, md2sens, sd1 = 1, sd2 = 2, var.equal = FALSE)
DF <- data.frame(MD = delta, performance = c(res[1,], res[2,]),
                 measure = c(rep("sensitivity", length(delta)),
                             rep("specificity", length(delta))))
ggplot(DF, aes(x = MD, y = performance, color = measure)) +
  geom_line() + ylim(0, 1.0) + xlab("mean difference") +
  scale_color_manual(values = c("darkblue", "darkred")) +
  ggtitle("SD1 = 1, SD2 = 2")
}
\keyword{univar}
