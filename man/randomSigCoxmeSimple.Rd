\name{randomSigCoxmeSimple}
\alias{randomSigCoxmeSimple}

\title{  Cox model using gene random signatures }
\description{
  Cox proportional hazard model with random and fixed effects using gene
  signatures
}
\usage{
randomSigCoxmeSimple(nda, randsign,  sigSize = NULL, mc.cores = 1,
    signature.method = "zscore", needScaling = FALSE, nite = 500, GS = NULL,
    adjusted.var.random = NA, adjusted.var.fixed = NA, adj.var.GS = "GSadj",
    executation.info = TRUE, perm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nda}{ \code{\link{ExpressionSet}}  object. \code{pData(nda)}
  must contain survival information in attributes \code{evn} (numeric
  vector with event description - 0 no event 1 event) and \code{tev}
  (time to event)  }
  \item{randsign}{NULL or list of gene signatures (names must match with
  \code{rownames(nda)}) }
  \item{sigSize}{If \code{nda = NULL}, random signatures of size
  \code{sigSize} are generated  }
  \item{mc.cores}{The number of cores to use for parallel executions}
  \item{signature.method}{Signature score method, to choose from
  \code{"zscore"}, \code{"gsva"} or \code{"plage"} }
  \item{needScaling}{If TRUE, the expression matrix in \code{nda} is
  scaled  }
  \item{nite}{Number of signatures (only if \code{randsign = NULL}) }
  \item{GS}{Global signature used as an adjusting variable as specified by
            attribute \code{adj.var.GS}. If NULL, no correction is applied}
  \item{adjusted.var.random}{Variables that enter in the model as random
  effect }
  \item{adjusted.var.fixed}{ Confounding variables that enter in the
  model as fixed effect }
  \item{adj.var.GS}{To choose among \code{'GSadj'} (default), where GS
  enters in the Cox model as confounding variable, \code{'GScor'}, where GS
  is deducted a priori from the expression matrix, and \code{'GSlmcor'},
  where gene signatures are adjusted a priori by GS using a linear model  }
  \item{executation.info}{ If \code{TRUE} information regarding
    executation stages is shown }
  \item{perm}{If \code{TRUE}, a permutation of the sample labels in the
  expression matrix is applied }

}
\details{
  An adjustment for global signature is provided to correct for
  technical (or global) artefacts that can bias the results of the
  model (see reference)
}
\value{
   list with outcome of cox models for tested gene signatures
   (contains lHR coefficient, HR coefficient, standard.error, statistic
   and p-value)
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
rs.s <- randomSigCoxmeSimple(nda.gse1456, randsign = NULL,
    sigSize = 50, mc.cores = 1, signature.method = "zscore", needScaling = TRUE,
    nite = 100, GS = NULL, executation.info = FALSE)
}


