\name{volcano}
\alias{volcano}
\alias{volcano.default}
\title{ Volcano Plots }
\description{
  Produce volcano plot(s) of the given effect size and p values. 
}
\usage{
volcano(effect, pval, \dots)

\S3method{volcano}{default}(effect, pval, effect0 = 0, sig.level = 0.05, 
                            effect.low = NULL, effect.high = NULL,
                            color.low = "#4575B4", color.high = "#D73027",
                            xlab = "effect size", ylab = "-log10(p value)",
                            title = "Volcano Plot", alpha = 1, shape = 19, \dots)
}
\arguments{
  \item{effect}{ numeric, measure of effect size.}
  \item{pval}{ numeric, (adjusted) p values.}
  \item{effect0}{ single numeric, value for no effect.}
  \item{sig.level}{ single numeric, significance level.}
  \item{effect.low}{ \code{NULL} or single numeric, boundary for low effect sizes.}
  \item{effect.high}{ \code{NULL} or single numeric, boundary for low effect sizes.}
  \item{color.low}{ color used if effect size smaller than \code{effect.low} and 
  (adjusted) p value smaller than \code{sig.level}.}
  \item{color.high}{ color used if effect size larger than \code{effect.high} and 
  (adjusted) p value smaller than \code{sig.level}.}
  \item{xlab}{ label of x-axis.}
  \item{ylab}{ label of y-axis.}
  \item{title}{ title of plot.}
  \item{alpha}{ bleding factor (default: no blending.}
  \item{shape}{ point shape used.}
  \item{\dots}{ further arguments that may be passed through).}
}
\details{
  The plot generates a \code{ggplot2} object that is shown. 

  Missing values are handled by the \code{ggplot2} functions.
}
\value{
  Object of class \code{gg} and \code{ggplot}.
}
\references{
  Wikipedia contributors, \emph{Volcano plot (statistics)}, Wikipedia, The Free Encyclopedia, 
  \url{https://en.wikipedia.org/w/index.php?title=Volcano_plot_(statistics)&oldid=900217316} 
  (accessed December 25, 2019). 

  For more sophisticated and flexible volcano plots see for instance:
  Blighe K, Rana S, Lewis M (2019). EnhancedVolcano: Publication-ready volcano 
  plots with enhanced colouring and labeling. R/Bioconductor package.
  \url{https://github.com/kevinblighe/EnhancedVolcano}.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
\examples{
## Generate some data
x <- matrix(rnorm(1000, mean = 10), nrow = 10)
g1 <- rep("control", 10)
y1 <- matrix(rnorm(500, mean = 11.75), nrow = 10)
y2 <- matrix(rnorm(500, mean = 9.75, sd = 3), nrow = 10)
g2 <- rep("treatment", 10)
group <- factor(c(g1, g2))
Data <- rbind(x, cbind(y1, y2))
pvals <- apply(Data, 2, function(x, group) hsu.t.test(x ~ group)$p.value,
               group = group)
## compute log-fold change
logfc <- function(x, group){
  res <- tapply(x, group, mean)
  log2(res[1]/res[2])
}
lfcs <- apply(Data, 2, logfc, group = group)

volcano(lfcs, pvals, xlab = "log-fold change")
volcano(lfcs, pvals, effect.low = -0.25, effect.high = 0.25, 
        xlab = "log-fold change")
volcano(lfcs, p.adjust(pvals, method = "fdr"), 
        effect.low = -0.25, effect.high = 0.25, 
        xlab = "log-fold change", ylab = "-log10(adj. p value)")
volcano(2^lfcs, pvals, effect0 = 1, effect.low = 1/2^0.25, effect.high = 2^0.25,
        xlab = "mean difference")
}
\keyword{hplot}
