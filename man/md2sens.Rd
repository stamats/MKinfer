\name{md2sens}
\alias{md2sens}
\title{ Compute Sensitivity and Specificty based on MD. }
\description{
  The function calculates the sensitivity and specificity for a given mean 
  difference (MD) and specified standard deviations of the two groups. 
  This is achieved by using the optimal cut-off poin based on Youden's J 
  statistic.
}
\usage{
md2sens(delta, sd1 = NULL, sd2 = NULL)
}
\arguments{
  \item{delta}{ (absolute) mean difference of the two groups. }
  \item{sd1}{ standard deviation of first group. }
  \item{sd2}{ standard deviation of second group. }
}
\details{
  This function can be used to demonstrate how (standardized) effect sizes 
  translate into sensitivity and specifity, based on two normal distributions.
  
  The two standard diviations are assumed to be equal if their difference is
  smaller than \code{1e-08}.
}
\value{Numeric vector with sensitivity and specificity.
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
## (standardized) mean difference to sensitivity/specificity
## equal variances
library(ggplot2)
delta <- seq(from = 0.0, to = 6, by = 0.05)
res <- sapply(delta, md2sens, sd1 = 1, sd2 = 1)
DF <- data.frame(SMD = delta, sensitivity = res[1,], 
                 specificity = res[2,])
ggplot(DF, aes(x = SMD, y = sensitivity)) +
  geom_line() + ylim(0.5, 1.0) + xlab("(standardized) mean difference") +
  ylab("sensitivity = specificity")

## unequal variances
delta <- seq(from = 0.0, to = 6, by = 0.05)
res <- sapply(delta, md2sens, sd1 = 1, sd2 = 2)
DF <- data.frame(MD = delta, performance = c(res[1,], res[2,]),
                 measure = c(rep("sensitivity", length(delta)),
                             rep("specificity", length(delta))))
ggplot(DF, aes(x = MD, y = performance, color = measure)) +
  geom_line() + ylim(0, 1.0) + xlab("mean difference") +
  scale_color_manual(values = c("darkblue", "darkred")) +
  ggtitle("SD1 = 1, SD2 = 2")
}
\keyword{univar}
