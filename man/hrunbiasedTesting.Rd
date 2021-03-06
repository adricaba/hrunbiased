\name{hrunbiasedTesting}
\alias{hrunbiasedTesting}
\alias{print.hrunbiasedTesting}

\title{  Cox proportional hazard model for gene signatures testing }
\description{
  Cox model for a target signature of genes allowing several correction
  types and testing procedures
}
\usage{
hrunbiasedTesting(nda, gene.signature, alternative ="two.sided",
           test.type = c("asymptotic", "random.signatures", "random.genes"),
           correction.type = c("NCsignatures", "Gsignature", "RUV", "none"),
           score = NULL, seqquant = 0.25,  nrs = 5000,
           comp.RUV = c(1,2), adj.var.GS = "GSadj",
           adjusted.var.random = NA, adjusted.var.fixed = NA, mc.cores = 1,
           signature.method = "zscore", maxGenes = 5000, executation.info=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nda}{ \code{\link{ExpressionSet}}  object. \code{pData(nda)}
  must contain survival information in attributes \code{evn} (numeric
  vector with event description - 0 no event 1 event) and \code{tev}
  (time to event) }
  \item{gene.signature}{Vector with the names of the gene signature to
  be tested (must match with \code{rownames(nda)})}
  \item{alternative}{ Character string specifying the alternative
  hypothesis, must be one of \code{"two.sided"} (default),
  \code{"greater"} or \code{"less"}}
  \item{test.type}{ Character string specifying the reference
  distribution to be used, must include at least one of
  \code{"asymptotic"} or \code{"random.signatures"} (all three by
  default)}
  \item{correction.type}{ Standard Cox model can be used when
\code{correction.type = "none"}, otherwise correction for global
signature \code{"Gsignature"}, correction for negative control
signatures \code{"NCsignatures"} or factor analysis \code{"RUV"}
 are available }
  \item{score}{Vector of the same length as number of rows in \code{nda}
  with scores to be used for  \code{"NCsignatures"} and \code{"RUV"}
  corrections }
  \item{seqquant}{Value between zero and one  used to set negative
  controls}
  \item{nrs}{ Number of random signatures to be used when
  \code{test.type = "random.signatures"}}
  \item{comp.RUV}{ Factors used in factor analysis signature correction
  when \code{correction.type = "RUV"}   }
  \item{adj.var.GS}{Only if correction.type is not "none" and
  siganture.method is "zscore". Transforming signature scores by
  subtracting the correction variable scores "GScor", residuals
  of linear model using the correction variable as covariate "GSlmcor",
  or using the correction variable in the Cox model "GSadj" (default)}
  \item{adjusted.var.random}{Variables that enter in the model as random
  effect }
  \item{adjusted.var.fixed}{ Confounding variables that enter in the
  model as fixed effect }
  \item{mc.cores}{The number of cores to use for parallel executions}
  \item{signature.method}{Signature score method, to choose from
  \code{"zscore"}, \code{"gsva"} or \code{"plage"} }
  \item{maxGenes}{ Number of genes to be used to determine reference
  distribution if \code{test.type = "random.genes"}}
    \item{executation.info}{ If \code{TRUE} information regarding
    executation stages is shown }
    \item{\dots}{
    Arguments passed to or from other methods to the low level.
  }
}
\details{
  Testing statistical significance of the relationship between a gene
  signature and recurrence can be done by considering both asymptotic
  distributions or random signature based distributions. The
  interpretation of the two is different: asymptotic testing finds
  evidence against the null hypothesis  HR = 1;  random signatures
  quantifies how special is the relation of a gene signature to
  recurrence in comparison to any possible random set of genes.

  A useful correction we consider is the adjusment by a global signature
  which uses the information of all genes and that may point towards
  technical / sampling directions. This global signature can be
  influenced by both technical and biological effects. In order to
  reduce the biological impact, adjusting variables are determined by
  negative control signatures based on a score ordered from low to high.
  For example, we suggest to use a measure of dispersion so that only
  low-variable genes are used in an  adjusting signature. Other scores
  based on cross-study reproducibility or any sort of biological prior
  information could also be employed. See
  \code{\link{hrunbiasedDiagnostic}} for more details in selecting the
  correction type.

  For \code{seqquant = 0}, then the correction type coincides with
  \code{correction.type = "none"} and for \code{seqquant = 1} then  the
  correction type coincides with  \code{correction.type = "Gsignature"}.
 }
\value{
   An object of class hrunbiasedTesting containing the following items
   \item{"hr.ts"}{Cox model outcome of target signature}
   \item{"hr.rs"}{Cox model outcome of random signatures}
   \item{"hr.tg"}{Cox model outcomes of genes composing the target signature}
   \item{"hr.rg"}{Cox model outcomes of random genes}
   \item{"pvals.rs"}{P-value found using as reference distribution the
  statistics obtained for random signatures }
   \item{"test.type"}{ Reference distribution used for hypothesis testing}
   \item{"gene.signature"}{ Gene signature used for hypothesis testing}
}
\author{ Adria Caballe Mestres }

\references{
Caballe Mestres A, Berenguer Llergo A and Stephan-Otto Attolini C.
Adjusting for systematic technical biases in risk assessment of gene signatures
in transcriptomic cancer cohorts. bioRxiv (2018).

Gagnon-Bartsch J.  and Speed T. Using control genes to correct for unwanted
variation in microarray data. Biostatistics 13, 539-552 (2012).

Hanzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for
microarray and RNA-Seq data. BMC Bioinformatics 7 (2013).
}

\seealso{\code{\link{plot.hrunbiasedTesting}},
  \code{\link{randomSigCoxmeSimple}} and
  \code{\link{hrunbiasedDiagnostic}} }
\examples{
eh <- ExperimentHub()
nda.brca <- query(eh, "mcsurvdata")[["EH1497"]]
nda.gse1456 <- nda.brca[,nda.brca$dataset=="GSE1456"]
data(mammaprint.sig)
mammaprint.sig <- intersect(mammaprint.sig, rownames(nda.gse1456))
correction.type <-  "Gsignature"
tout.GS <- hrunbiasedTesting(nda.gse1456,
    gene.signature = mammaprint.sig,
    test.type = c("asymptotic", "random.signatures", "random.genes"),
    correction.type = "Gsignature", mc.cores=1,
    nrs = 100,  executation.info = FALSE)
tout.GS

}

