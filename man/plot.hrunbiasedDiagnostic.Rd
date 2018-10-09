\name{plot.hrunbiasedDiagnostic}
\alias{plot.hrunbiasedDiagnostic}

\title{    Cox model diagnostic plots using gene random signatures}
\description{
  Plot function for hrunbiasedDiagnostic objects
}
\usage{
\method{plot}{hrunbiasedDiagnostic}(x, diagnostic.plot, id.size = 1,
   show.stat = TRUE, legend.out = TRUE, perms.to.show = 1:10,
   seq.alpha = seq(0,.5, length.out=50), seq.alpha.pow = c(0.01,0.05),
   power.whplots = NULL, col1 = "gray", col2 = "orange",
   pointsq = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{hrunbiasedDensityRS},
  \code{hrunbiasedPermutations}, \code{hrunbiasedSimulations} or
  \code{hrunbiasedPCpower} }
  \item{diagnostic.plot}{At least one of \code{"density.rs"},
  \code{"density.rs.size"}, \code{"geneMean.geneSign"} if \code{x} is
  of class
     \code{hrunbiasedDensityRS}, \code{"perm.violin"},
   \code{"perm.GSvsEvents"}, \code{"perms.corr.GS"} if \code{x} is of
   class \code{hrunbiasedPermutations},
      \code{"simulations"} if \code{x} is of class
    \code{hrunbiasedSimulations}, and \code{"positive.cont.power"},
    \code{"bias.power"}, \code{"corNCsig.power"}
     if \code{x} is of class \code{hrunbiasedPCpower},
    \code{"PCAcorrelation"}, if \code{x} is of class 
    \code{hrunbiasedPCAcorrelation}}
  \item{id.size}{For class \code{hrunbiasedDensityRS} related
  plots. Show random signature distribution for signature size
  \code{sigSize[id.size]}}
  \item{show.stat}{For \code{"density.rs"} plots. If TRUE, hazard ratio
  statistic is shown, if FALSE, log hazard ratio is shown}
  \item{legend.out}{If TRUE, an automatic legend is shown in the plot}
  \item{perms.to.show}{For \code{"permutations"} diagnostic
  plots. Number of instances to show (keep it low for
  interpretability)}
  \item{seq.alpha}{For \code{"positive.cont.power"} diagnostic
  plots. Sequence of rejecting levels to be used to plot power curves
  }
  \item{seq.alpha.pow}{For \code{"corNCsig.power"} or
  \code{"bias.power"} diagnostic plots. Sequence of rejecting levels
  to be used to summarize the power using several negative control
  levels}
  \item{power.whplots}{For class \code{hrunbiasedPCpower} related
  plots. Show power plots for contaminated positive control signature
  defined by \code{prop.pos.cont[power.whplots]}}
  \item{col1}{For \code{"simulations"} diagnostic plot. Color for
  non-adjusted random signature distributions }
  \item{col2}{For \code{"simulations"} diagnostic plot. Color for
  GS-adjusted random signature distributions }
  \item{pointsq}{If \code{TRUE}, population logHR for target signature
  is shown}
  \item{\dots}{
    Arguments passed to or from other methods to the low level.
  }
}
\details{

  Following details section in  \code{\link{hrunbiasedDiagnostic}},
  module (i) contains \code{"density.rs"}: hazard ratio
  t-statistic density plot obtained from random signatures;
  \code{"density.rs.size"}: density plot for hazard ratio t-statistic as
  function of signature size; \code{"geneMean.geneSign"}: average log HR
  genewise in a signature vs log HR in the signature. Module (ii) contains
  \code{"perm.violin"}: hazard ratio violin plots obtained from random
  signatures for several instances with a shuffle in the time-to-even
  outcome; \code{"perm.GSvsEvents"}: global signature boxplots that
  distinguish between event and not event; \code{"perms.corr.GS"}:
  relationship between GS event and not event average difference and
  observed average lHR. Module (iii) contains  \code{"positive.cont.power"}:
  power curves using positive control random signatures;
  \code{"bias.power"}: two y-axis plot showing average values for the hazard
  ratio t-statistic and power given several charactarizations of negative
  controls; \code{"corNCsig.power"}: two y-axis plot showing correlation to
  positive controls and power of several charactarizations of negative
  controls; Module (iv) contains \code{"simulations"}: random signatures
  hazard ratio distributions for simulated time-to-event data.
  }
  \value{
    plot of object of class \code{\link{hrunbiasedDiagnostic}}.
  }
\author{ Adria Caballe Mestres }

\references{
Caballe Mestres A, Berenguer Llergo A and Stephan-Otto Attolini C. 
Adjusting for systematic technical biases in risk assessment of gene signatures 
in transcriptomic cancer cohorts. bioRxiv (2018). 
}

\seealso{ \code{\link{hrunbiasedDiagnostic}}}
\examples{
data(out.brca)
plot(out.brca[["metabric"]], diagnostic.plot = "density.rs", main = "metabric",
     xlab = "lhr.stat", ylab="", lwd=3, legend.out = FALSE, id.size = 3 )

}

%  \item{globalSignatureEvents}{ Whether to boxplot global signature values vs
%  recurrence information (shuffled)}

