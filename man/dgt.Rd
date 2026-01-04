\name{dgt}
\alias{dgt}
\alias{pgt}
\alias{qgt}
\title{A Generalized Central t Distribution}
\description{
  Density, distribution function and quantile function for the generalized 
  central t distribution defined in Section 5.3 of Xiao (2018)).
}
\usage{
dgt(x, n1, n2, v1tov2, log = FALSE, MAX = 1e4)
pgt(q, n1, n2, v1tov2, lower.tail = TRUE, log.p = FALSE, 
    rel.tol = .Machine$double.eps^0.25, parallel = FALSE, cl = NULL)
qgt(p, n1, n2, v1tov2, lower.tail = TRUE, log.p = FALSE, 
    tol = .Machine$double.eps^0.5, parallel = FALSE, cl = NULL)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n1}{sample size of first group.}
  \item{n2}{sample size of second group.}
  \item{v1tov2}{ratio of variances (var1/var2).}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{MAX}{maximum number of iterations to determine the hypergeometric function; see 
    function \code{\link[hypergeo]{hypergeo}}.}
  \item{rel.tol}{relative accuracy of \code{\link{integrate}}.}
  \item{tol}{the desired accuracy (convergence tolerance) of \code{\link{uniroot}}.}
  \item{parallel}{logical; use parellel computing via function \code{\link{parRapply}}.}
  \item{cl}{a cluster object, created by package \code{\link[parallel:makeCluster]{parallel}} 
    or by package \code{\link[snow:makeCluster]{snow}}. If \code{NULL}, a new cluster 
    will be generated with function \code{\link[parallel:makeCluster]{makeCluster}} 
    using the maximum number of available cores - 1.}
}
\value{
  \code{dgt} gives the density,
  \code{pgt} gives the distribution function, and
  \code{qgt} gives the quantile function
}
\details{
  It is possible to provide arguments of length larger than 1 to the functions.
  
  The density of the distribution can be found on page 58 of Xiao (2018). If the
  hypergeometric function does not converge, one can try to increase the value 
  of \code{MAX}.
  
  The distribution function is determined by numerical intergration via function
  \code{\link{integrate}}.
  
  The quantile function is determined by applying function \code{\link{uniroot}}.
}
\references{
  Y. Xiao (2018). On the Solution of a Generalized Behrens-Fisher Problem. 
  \emph{Far East Journal of Theoretical Statistics}, \bold{54} (1), 21-140.
}
\seealso{\code{\link{dt}}}
\examples{
  ## equal variances
  dgt(0, n1 = 5, n2 = 10, v1tov2 = 1)
  dt(0, df = 5+10-2)
  
  pgt(2, n1 = 5, n2 = 10, v1tov2 = 1)
  pt(2, df = 5+10-2)
  
  qgt(0.975, n1 = 5, n2 = 10, v1tov2 = 1)
  qt(0.975, df = 5+10-2)
  
  library(ggplot2)
  ggplot(data = data.frame(x = c(-4, 4)), aes(x)) +
      stat_function(fun = function(x) dgt(x, n1 = 5, n2 = 10, v1tov2 = 1), 
                    n = 101, col = "darkblue", lwd = 1) + 
      stat_function(fun = function(x) dt(x, df = 5+10-2), 
                    n = 101, col = "darkgreen", lwd = 1) + ylab("density") +
      ggtitle("Generalized Central (blue) t vs. Student t Distribution (green)")
                    
  
  ## unequal variances
  nu <- (1^2/5+5^2/10)^2/(1^4/(5^2*(5-1)) + 5^4/(10^2*(10-1)))
  dgt(0, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2)
  dt(0, df = nu)
  
  pgt(2, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2)
  pt(2, df = nu)
  
  qgt(0.975, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2)
  qt(0.975, df = nu)
  
  ggplot(data = data.frame(x = c(-4, 4)), aes(x)) +
      stat_function(fun = function(x) dgt(x, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2), 
                    n = 101, col = "darkblue", lwd = 1) + 
      stat_function(fun = function(x) dt(x, df = nu), 
                    n = 101, col = "darkred", linetype = "dashed", lwd = 1) + 
      ylab("density") +
      ggtitle("Generalized Central t (blue) vs. Welch t Distribution (red)")
  
  \dontrun{
    ## parallel computing
    pgt(seq(from = -2, to = 2, length = 30), n1 = 5, n2 = 10, v1tov2 = 1^2/5^2, 
        parallel = TRUE)
    qgt(seq(from = 0.1, to = 0.9, length = 30), n1 = 5, n2 = 10, v1tov2 = 1^2/5^2, 
  +     parallel = TRUE)
  }    
}
\keyword{distribution}
