\name{md2zfactor}
\alias{md2zfactor}
\title{ Compute z-Factor based on MD. }
\description{
  The function calculates the z-factor for a given mean 
  difference (MD) and specified standard deviations of the two groups. 
}
\usage{
md2zfactor(delta, sd1 = 1, sd2 = 1)
}
\arguments{
  \item{delta}{ (absolute) mean difference of the two groups. }
  \item{sd1}{ standard deviation of first group. }
  \item{sd2}{ standard deviation of second group. }
}
\details{
  This function can be used to demonstrate how (standardized) effect sizes 
  translate into z-factor, based on two normal distributions.
}
\value{Numeric vector with sensitivity and specificity.
}
\references{
  J.H. Zhang, T.D.Y. Chung and K.R. Oldenburg (1999). A simple statistical 
  parameter for use in evaluation and validation of high throughput screening assays. 
  \emph{Journal of Biomolecular Screening} \bold{4} (2): 67-73.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\examples{
## (standardized) mean difference to sensitivity/specificity
## equal variances
library(ggplot2)
delta <- seq(from = 2, to = 18, by = 0.05)
res <- sapply(delta, md2zfactor)
DF <- data.frame(SMD = delta, zfactor = res)
ggplot(DF, aes(x = SMD, y = zfactor)) +
  geom_line() + xlab("(standardized) mean difference") +
  ylab("z-factor") + ggtitle("SD1 = SD2 = 1") + 
  geom_hline(yintercept = 1, linetype = "dotted")

## unequal variances
delta <- seq(from = 2.5, to = 20, by = 0.05)
res <- sapply(delta, md2zfactor, sd1 = 1, sd2 = 2)
DF <- data.frame(MD = delta, zfactor = res)
ggplot(DF, aes(x = MD, y = zfactor)) +
  geom_line() + xlab("mean difference") +
  ylab("z-factor") + ggtitle("SD1 = 1, SD2 = 2") +
  geom_hline(yintercept = 1, linetype = "dotted")
}
\keyword{univar}
