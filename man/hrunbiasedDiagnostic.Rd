\name{hrunbiasedDiagnostic}
\alias{hrunbiasedDiagnostic}
\alias{hrunbiasedDensityRS}
\alias{hrunbiasedPermutations}
\alias{hrunbiasedPCpower}
\alias{hrunbiasedSimulations}
\alias{hrunbiasedPCAcorrelation}
\alias{print.hrunbiasedDiagnostic}

\title{  Cox model diagnostics using gene random signatures }
\description{
  Contains several diagnostic functions for adjusted Cox proportional hazard
  models using gene signatures
}
\usage{

hrunbiasedDiagnostic(nda, correction.type = "NCsignatures",
        score, signature.method = "zscore",
        seqquant = c(0,seq(0.05,1,length.out=6)), sigSize = 50,
        diagnostic.type = c("density.rs", "permutations", "simulations",
        "positive.cont.power", "PCAcorrelation"), ...)

hrunbiasedDensityRS(nda, correction.type = "NCsignatures", score,
        seqquant = c(0,seq(0.05,1,length.out=6)), sigSize = 50,
        nrs = 5000, comp.RUV = c(1,2), adj.var.GS = "GSadj",
        adjusted.var.random = NA, adjusted.var.fixed= NA,
        executation.info = TRUE, mc.cores = 1,  GSs = NULL, rg.cox = TRUE,
        id.size = 1, basic.cont = TRUE, signature.method = "zscore", ...)

hrunbiasedPermutations(nda, correction.type = "NCsignatures",
       score, seqquant = c(0,seq(0.05,1,length.out=6)), sigSize = 50,
        nrs = 5000,  comp.RUV = c(1,2), nperm = 100,
       adj.var.GS = "GSadj", permutation.type = "none", GSassociation = TRUE,
       executation.info = TRUE, mc.cores = 1, GSs = NULL,  basic.cont = TRUE,
       GS = NULL, signature.method = "zscore",  ...)

hrunbiasedPCpower(nda, positive.controls, correction.type = "NCsignatures", score,
     seqquant = c(0,seq(0.05,1,length.out=6)), alternative = "two.sided",
     bootstrapping = TRUE,  prop.pos.cont = 1, sigSize = 50,
     id.size = 1, nrs = 5000, nrpcs = 500, comp.RUV = c(1,2),
     null.dist = "asymptotic",  adj.var.GS = "GSadj", adjusted.var.random = NA,
     adjusted.var.fixed= NA, executation.info = TRUE, mc.cores = 1, GSs = NULL,
     cdRS = NULL, basic.cont = TRUE, ...)

hrunbiasedSimulations(nda, positive.controls, GScoef = 0,
     tscoef = 0, ninst.rs = 10, ninst.ts = 500, nrs = 500, lambdas = c(0.1,0.2),
     gammas = c(1.5,2),  mc.cores = 1,  stat = TRUE,
     executation.info = TRUE, GSs = NULL, signature.method = "zscore",
     basic.cont = TRUE, sigSize = 50, correction.type = "none", ...)

hrunbiasedPCAcorrelation(nda, nrs = 100, sigSize = 50, mc.cores = 1,
      comp.PCA = c(1,2), GS = NULL, executation.info = TRUE,
     signature.method = "zscore", ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nda}{ \code{\link{ExpressionSet}}  object. \code{pData(nda)}
 must contain survival information in attributes \code{evn} (numeric
 vector with event description - 0 no event 1 event) and \code{tev}
 (time to event)}
 \item{correction.type}{ Standard Cox model can be used when
 \code{correction.type = "none"}, otherwise correction for global
 signature \code{"Gsignature"}, correction for negative control
 signatures \code{"NCsignatures"} (default) or factor analysis on
 negative control genes \code{"RUV"} are available }
  \item{score}{Vector of the same length as number of rows in \code{nda}
 with scores to be used for  \code{"NCsignatures"} and \code{"RUV"}
 corrections (see details)}
  \item{signature.method}{Signature score method: to choose from
  \code{"zscore"}, \code{"gsva"} or \code{"plage"} }
  \item{seqquant}{Numeric vector with values between zero and one  used
 to set negative controls}
  \item{diagnostic.type}{Diagnostic to be performed: at least one among
 \code{"density.rs"},  \code{"permutations"},
 \code{"positive.cont.power"}, \code{"simulations"} or
 \code{"PCAcorrelation"} (see details)}
  \item{sigSize}{Numeric vector of at least one element with signature
 sizes to be considered (a single signature size is used for
 \code{diagnostic.type = "positive.cont.power"} following attribute
 \code{id.size})}
  \item{nrs}{ Number of random signature instances used for approximated
 null distribution}

  \item{comp.RUV}{ Factors used in factor analysis signature correction
  when \code{correction.type = "RUV"}   }
  \item{GSs}{List with correction variables. If NULL, GSs is found by
        the approach specified in \code{correction.type} }

  \item{adj.var.GS}{Only if correction.type is not "none" and
  siganture.method is "zscore". Transforming signature scores by
  subtracting the correction variable scores "GScor", residuals
  of linear model using the correction variable as covariate "GSlmcor",
  or using the correction variable in the Cox model "GSadj" (default)}
  \item{adjusted.var.random}{Variables that enter in the Cox model as
  random effect }
  \item{adjusted.var.fixed}{ Confounding variables that enter in the Cox
  model as fixed effect }
  \item{executation.info}{ If \code{TRUE} information regarding
  executation stages is shown }
  \item{mc.cores}{The number of cores to use for parallel executions}
  \item{rg.cox}{If TRUE, survival effects at gene level are computed}

  \item{alternative}{Character string specifying the alternative
        hypothesis, must be one of \code{"two.sided"} (default),
        \code{"greater"} or \code{"less"}}
  \item{nperm}{ Number of outcome shuffle instances to be done }
  \item{permutation.type}{ Either \code{"none"}, \code{"Gsignature"} or
        \code{"norm01"} (see details) }
  \item{GSassociation}{ GS association with survival is computed to show
        global tendency}
  \item{GS}{List with the global signature (for internal use)}
  \item{positive.controls}{ Character string with gene signature to be
        used as positive control (must match \code{rownames(nda)})}
  \item{bootstrapping}{If \code{TRUE}, random positive control
        signatures of size \code{sigSize} are sampled with replacement}
  \item{prop.pos.cont}{ Numeric vector of at least one element that
        indicates the amount (proportion between 0 and 1) of purity of
        the positive control signatures to be tested (one by default,
        zero would correspond to random signatures, and in betweeen
        there would be contaminated positive control signatures)}
  \item{nrpcs}{ Number of random positive control signature instances to
        be used to approximate power}
  \item{id.size}{Signature size used \code{sigSize[id.size]} for some
        random signature based diagnostics}
  \item{null.dist}{ Either \code{"asymptotic"} (default) or
        \code{"random.signatures"}}
  \item{cdRS}{A \code{hrunbiasedDensityRS} object (to avoid computation
        burden) or \code{NULL}}
  \item{basic.cont}{Basic controls are performed, these are done no
        matter the \code{basic.cont} specification if
        hrunbiasedDiagnostic is called}
  \item{GScoef}{The log HR effect of global signature (common tendency)
        for simulations}
  \item{tscoef}{The log HR effect of target signature (unique tendency)
        for simulations}
  \item{ninst.rs}{Instances of random signature distributions for
        simulated data analysis (must be larger than zero)}
  \item{ninst.ts}{Instances of target signature estimated effect for
        simulated data analysis (can be zero)}
  \item{lambdas}{Weibull distribution shape parameters (see
        \code{\link{simsurv}}), first element for times to event distribution
        and second element for follow-up times distribution}
  \item{gammas}{Weibull distribution scale parameters (see
        \code{\link{simsurv}}), first element for times to event distribution
        and second element for follow-up times distribution}
  \item{stat}{If \code{TRUE}, simulated data analysis is performed on
    the HR statistic instead of the logHR coefficient}
  \item{comp.PCA}{ Numeric(2) vector with the specific components to
    compute correlation between random signatures and the projections}
  \item{\dots}{
    Arguments passed to or from other methods to the low level.
  }
}
\details{
  Several diagnostic plots via \code{\link{plot.hrunbiasedDiagnostic}}
        are available which can be useful to decide  whether an
        adjustment in estimating hazard ratios for gene signatures is
        needed and, in some extent, which class of correction could be
        appropriate. The diagnostics implemented in this function can be
        grouped in five modules: (i) density random signatures; (ii)
        permutated outcome;  (iii) positive control power; (iv)
        simulated survival data; and (v) random signatures correlation
        with PC. Detailed information about the plots in each of the
        modules is provided in \code{\link{plot.hrunbiasedDiagnostic}}.

  Correction types include adjustments by global signature
  (\code{"Gsignature"}),  negative control signatures
  (\code{"NCsignatures"}), or negative control factors (\code{"RUV"}).
  Global signature is a summary across all standardized genes in the
  data whereas negative control adjustments are found by summarising the
  information  of only  a few particular genes.    All genes are given a
  score via  attribute \code{score} (i.e., gene variability, cross-study
  reproducibility,...).
  Selected negative control genes are those whose score is smaller than
  the percentile indicated in \code{seqquant}. Several
  values for \code{seqquant} can be given so the diagnostic can be done
  over different specifications of negative control genes.

  Permutations based diagnostic is done to show global biases in random
  signature distributions due to both sampling and gene-to-gene
  dependence structure.
  The \code{permutation.type} equal to "none" or "Gsignature" have the
  same interpretation as the one for the correction type. Besides, in
  "norm01" we also  consider a shuffling of the samples independently
  for all genes in the data so gene-to-gene relationships are broken.
}
\value{
  If \code{diagnostic.type = "density.rs"} or hrunbiasedDensityRS function
  is called, return an object of class  \code{hrunbiasedDensityRS}  with
  attributes
  \item{"hr.rs"}{Cox model outcome of random signatures}
  \item{"geneMeanHR.signatureHR"}{ Average log HR genewise in random
  signature and log HR in the same random signatures }

  If \code{diagnostic.type = "permutations"} or hrunbiasedPermutations function
  is called, return an object of class  \code{hrunbiasedPermutations}  with
  attributes
   \item{"perms"}{Time-to-event data permuted indexes}
   \item{"perm.out"}{Outcome of permutations study}

  If \code{diagnostic.type = "positive.cont.power"} or hrunbiasedPCpower
  function is called, return an object of class  \code{hrunbiasedPCpower} with
  attributes
   \item{"hr.pc"}{Cox model outcome of positive control signature}
   \item{"pvals.pc"}{ P-values found for positive control signatures}

  If \code{diagnostic.type = "simulations"} or hrunbiasedSimulations function
  is called, return an object of class  \code{hrunbiasedSimulations}  with
  attributes
   \item{"simu.info"}{Simulated data analysis output for diagnostic plot}

  If \code{diagnostic.type = "PCAcorrelations"} or hrunbiasedPCAcorrelation
  function is called, return an object of class \code{hrunbiasedPCAcorrelation}
  with attributes
   \item{"rs.cor"}{Correlation between random signatures and the
  projections of the comp.PCA components}
   \item{"rs.gs.cor"}{Correlation between random signatures adjusted by
  GS and the projections of the comp.PCA components}
   \item{"gs.cor"}{Correlation between Global signature and the
    projections of the comp.PCA components}

  For all called functions or \code{diagnostic.type} options, the
  following items are also included
   \item{"GS"}{Global signature}
   \item{"seqquant"}{Sequence seqquant used}
   \item{"sigSize"}{Signature size}
   \item{"correction.type"}{Correction approach used}
   \item{"diagnostic.type"}{Diagnostic type used}

 }

\author{ Adria Caballe Mestres }

\references{
Caballe Mestres A, Berenguer Llergo A and Stephan-Otto Attolini C.
Adjusting for systematic technical biases in risk assessment of gene signatures
in transcriptomic cancer cohorts. bioRxiv (2018).

Gagnon-Bartsch J.  and Speed T. Using control genes to correct for
unwanted variation in microarray data. Biostatistics 13, 539-552 (2012).

Hanzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for
microarray and RNA-Seq data. BMC Bioinformatics 7 (2013).
}

\seealso{\code{\link{plot.hrunbiasedDiagnostic}} and
  \code{\link{hrunbiasedTesting}}}
\examples{
eh <- ExperimentHub()
nda.brca <- query(eh, "mcsurvdata")[["EH1497"]]
nda.gse1456 <- nda.brca[,nda.brca$dataset=="GSE1456"]

out.brca.GSE1456 <- hrunbiasedDiagnostic(nda.gse1456,
 correction.type = "none", mc.cores = 1, diagnostic.type = "density.rs",
 sigSize = 50, nrs = 100, id.size = 1, executation.info = FALSE)
}


