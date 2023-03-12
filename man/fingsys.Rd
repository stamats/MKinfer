\name{fingsys}
\alias{fingsys}
\docType{data}
\title{ Systolic Blood Pressure Dataset used in Bland and Altman (1995) }
\description{
  This dataset was used in Band and Altman (1995) to demonstrate why plotting 
  difference against standard method is misleading.
}
\usage{data(fingsys)}
\format{
  A data.frame with two colums
  \describe{
    \item{\code{armsys}}{ systolic blood pressure (mm Hg) measured by a standard arm cuff. }
    \item{\code{fingsys}}{ systolic blood pressure (mm Hg) measured by a finger monitor. }
  }
  For more details see Close et al. (1986) as well as Bland and Altman (1995).
}
\details{
  The dataset is a random subset of 200 observations from a larger dataset
  of Close et al. (1986) that was used in Bland and Altman (1995) to demonstrate
  why plotting difference against standard method is misleading.
}
\source{
  The data set was obtained from
  \url{https://www-users.york.ac.uk/~mb55/datasets/fingsys.dct}
}
\references{
  Close A, Hamilton G, Muriss S (1986). Finger systolic pressure: its use in screening 
  for hypertension and monitoring. \emph{Brit Med J} 293:775-778.
  
  Bland JM, Altman DG. (1995) Comparing methods of measurement: why plotting difference against 
  standard method is misleading. \emph{Lancet} 346, 1085-7. 
}
\examples{
  data(fingsys)
  str(fingsys)
  head(fingsys)
  tail(fingsys)
}
\keyword{datasets}
