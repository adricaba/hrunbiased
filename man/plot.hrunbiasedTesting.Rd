\name{plot.hrunbiasedTesting}
\alias{plot.hrunbiasedTesting}

\title{  Plotting Cox proportional hazard model outcome for gene signatures
  testing }
\description{
  Plot function for hrunbiasedTesting objects
}
\usage{
\method{plot}{hrunbiasedTesting}(x, legend.out = TRUE, show.pvalues = TRUE,
          show.stat = TRUE,  density.tg.minpoints = 10, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{hrunbiasedTesting}}
  \item{legend.out}{If TRUE, an automatic legend is shown in the plot}
  \item{show.pvalues}{ If TRUE, an automatic legend with the testing
     p-values is created}
  \item{show.stat}{If TRUE, the hazard ratio statistic is shown,
    if FALSE, the log hazard ratio is shown}
  \item{density.tg.minpoints}{minimum number of genes to show a
	density estimation of the signature hazard ratios genewise}
  \item{\dots}{
    Arguments passed to or from other methods to the low level
  }
}
  \value{
   plot of object of class \code{\link{hrunbiasedTesting}}.
  }

\author{ Adria Caballe Mestres }

\references{
Caballe Mestres A, Berenguer Llergo A and Stephan-Otto Attolini C.
Adjusting for systematic technical biases in risk assessment of gene signatures
in transcriptomic cancer cohorts. bioRxiv (2018).
}
\seealso{\code{\link{hrunbiasedTesting}} and
   \code{\link{hrunbiasedDiagnostic}}}
\examples{
eh <- ExperimentHub()
nda.brca <- query(eh, "mcsurvdata")[["EH1497"]]
nda.gse1456 <- nda.brca[,nda.brca$dataset=="GSE1456"]
data(mammaprint.sig)
mammaprint.sig <- intersect(mammaprint.sig, rownames(nda.gse1456))
correction.type <-  "Gsignature"
tout.GS <- hrunbiasedTesting(nda.gse1456,
    gene.signature = mammaprint.sig, test.type = c("asymptotic",
    "random.signatures", "random.genes"), correction.type = "Gsignature",
    mc.cores=1, nrs = 100, executation.info = FALSE)
plot(tout.GS)
}

