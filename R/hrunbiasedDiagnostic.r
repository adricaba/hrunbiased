## function hrunbiasedDiagnostic
hrunbiasedDiagnostic <- function(nda, correction.type = "NCsignatures",
                  score, signature.method = "zscore",
                  seqquant = c(0,seq(0.05,1,length.out=6)), sigSize = 50,
                  diagnostic.type = c("density.rs", "permutations",
                  "simulations", "positive.cont.power",
                  "PCAcorrelation"), ...){
  ## controls
  correction.type <- basic.controls(nda = nda,
         correction.type = correction.type, sigSize = sigSize)

  diagnostic.type <- diagnostic.controls(diagnostic.type)
  signature.method <- signmethod.control(signature.method)
   if(!all(diagnostic.type %in% c("simulations","PCAcorrelation"))){
    if(missing(score) & correction.type %in% c("none", "Gsignature"))
          score = NULL
    if(missing(score) & !correction.type %in% c("none", "Gsignature"))
          stop("score must be specified")
  }

  ## initialization
  out <- list()

  class.out <- character()
  densRS.info <- NULL

  eps <- 10e-8
  if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

  ## diagnostics
  if("simulations" %in% diagnostic.type){
      out$simu.info <- hrunbiasedSimulations(nda,
                                    signature.method = signature.method,
                                    basic.cont = FALSE,...)$simu.info
      class.out <- c(class.out, class(out$simu.info)[1])
  }

  out$GS <- GSnegControls(nda = nda, correction.type =  "Gsignature",
              score = NULL, seqquant = seqquant,  ...)
  if(correction.type =="Gsignature" | all(diagnostic.type %in%
         c("simulations","PCAcorrelation")))
     GSs <- out$GS
  else
     GSs <- GSnegControls(nda = nda, correction.type =  correction.type,
                score = score, seqquant = seqquant,  ...)
  if("permutations" %in% diagnostic.type){
    perm.info <- hrunbiasedPermutations(nda = nda,
                   correction.type = correction.type, score = score,
                   seqquant = seqquant, GSs = GSs, sigSize = sigSize,
                   basic.cont = FALSE, GS = out$GS,
                   signature.method = signature.method, ...)
    out$perm.out <- perm.info$perm.out
    out$perms <- perm.info$perms
    class.out <- c(class.out, class(perm.info)[1])
  }
  if("density.rs" %in% diagnostic.type){
    densRS.info <- hrunbiasedDensityRS(nda, correction.type = correction.type,
                       score = score, seqquant = seqquant, GSs = GSs,
                       sigSize = sigSize,  basic.cont = FALSE,
                       signature.method = signature.method, ...)
    out$hr.rs <- densRS.info$hr.rs
    out$randsign <- densRS.info$randsign
    out$geneMeanHR.signatureHR <- densRS.info$geneMeanHR.signatureHR
    class.out <- c(class.out, class(densRS.info)[1])
  }
  if("positive.cont.power" %in% diagnostic.type){
   pcPower.info <- hrunbiasedPCpower(nda, correction.type = correction.type,
                      score = score, seqquant = seqquant, GSs = GSs,
                      cdRS = densRS.info, sigSize = sigSize,
                      basic.cont = FALSE, signature.method = signature.method, ...)

   out$hr.pc <- pcPower.info$hr.pc
   out$pvals.pc <- pcPower.info$pvals.pc
   class.out <- c(class.out, class(pcPower.info)[1])

   attr(out,"null.dist") <-  attr(pcPower.info, "null.dist")
   if(correction.type %in% c("NCsignatures","RUV") &  length(seqquant) > 1){
          attr(out, "cor.nc.rs") <- attr(pcPower.info, "cor.nc.rs")
          attr(out, "means.rs") <- attr(pcPower.info, "means.rs")
    }

  }
  if("PCAcorrelation" %in% diagnostic.type){
   PCAcorr.info <- hrunbiasedPCAcorrelation(nda, sigSize = sigSize,
                           GS = out$GS, signature.method = signature.method, ...)
   out$rs.cor <-  PCAcorr.info$rs.cor
   out$rs.gs.cor <-  PCAcorr.info$rs.gs.cor
   out$gs.cor <-  PCAcorr.info$gs.cor
   attr(out,"comp.PCA") <- attr(PCAcorr.info, "comp.PCA")
   class.out <- c(class.out, class(PCAcorr.info)[1])
  }

  out$sigSize <- sigSize
  out$correction.type <- correction.type
  out$diagnostic.type <- diagnostic.type
  out$seqquant <- seqquant
  class(out) <- c(class.out,"hrunbiasedDiagnostic")
  attr(out,"ev.info") <- nda$evn

  return(out)
}

##### hrunbiasedDensityRS
hrunbiasedDensityRS <- function(nda, correction.type = "NCsignatures", score,
          seqquant = c(0,seq(0.05,1,length.out=6)), sigSize = 50,
          nrs = 5000, comp.RUV = c(1,2), adj.var.GS = "GSadj",
          adjusted.var.random = NA, adjusted.var.fixed= NA,
          executation.info = TRUE, mc.cores = 1,  GSs = NULL, rg.cox = TRUE,
          id.size = 1, basic.cont = TRUE, signature.method = "zscore", ...){
 # controls
 if(basic.cont){
   correction.type <- basic.controls(nda = nda,
          correction.type = correction.type, sigSize = sigSize)

   if(missing(score) & correction.type %in% c("none", "Gsignature"))
          score = NULL
   if(missing(score) & !correction.type %in% c("none", "Gsignature"))
       stop("score must be specified")
    signature.method <- signmethod.control(signature.method)
   if(correction.type %in% c("NCsignatures","Gsignature") & signature.method=="zscore")
       adj.var.GS <- adjvarGS.control(adj.var.GS)
 }
 mod.controls(nda, adjusted.var.random, adjusted.var.fixed, GSs)
 if(id.size > length(sigSize))
     stop("id.size must be between 1 and length(sigSize)")

 # initialization
 rs.HR <- NULL
 geneMeanHR.signatureHR <- NULL

 eps <- 10e-8
 if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

 if(is.null(GSs))
  GSs <- GSnegControls(nda, correction.type =  correction.type, score=score,
                 seqquant = seqquant,  comp.RUV = comp.RUV, ...)

 # diagnostics
 randsign <- list()
 if(executation.info) cat("\n random signatures cox models \n")
 lz <- length(sigSize)
 rs.HR.size <- list()
 for(K in seq_len(lz)){
  randsign[[K]] <- lapply(seq_len(nrs), function(j) gs <- sample(rownames(exprs(nda)),
                       sigSize[K]))
  rs.HR.size[[K]] <- rsCoxmeByCorrections(nda,
        correction.type = correction.type, score=score, seqquant=seqquant,
        randsign=randsign[[K]], adj.var.GS = adj.var.GS,
        adjusted.var.random = adjusted.var.random, signature.method = signature.method,
        adjusted.var.fixed = adjusted.var.fixed, mc.cores = mc.cores,
        GSs = GSs, executation.info=executation.info, ...)
 }

 if(rg.cox){
  if(executation.info) cat("\n random gene cox models \n")
  geneMeanHR.signatureHR <- geneMeanHRvsSignatureHR(nda,
        correction.type = correction.type, score = score, seqquant = seqquant,
        randsign = randsign[[id.size]], adj.var.GS = adj.var.GS,
        adjusted.var.random = adjusted.var.random,
        adjusted.var.fixed = adjusted.var.fixed,  mc.cores = mc.cores,
        GSs = GSs, rs.HR =rs.HR.size[[id.size]],
        executation.info =  executation.info)
 }

 out <- list(hr.rs = rs.HR.size, randsign = randsign, sigSize = sigSize,
             geneMeanHR.signatureHR = geneMeanHR.signatureHR )
 class(out) <- c("hrunbiasedDensityRS", "hrunbiasedDiagnostic")
 return(out)
}


##### hrunbiasedPermutations
hrunbiasedPermutations <- function(nda, correction.type = "NCsignatures",
      score, seqquant = c(0,seq(0.05,1,length.out=6)), sigSize = 50,
      nrs = 5000,  comp.RUV = c(1,2), nperm = 100,
      adj.var.GS = "GSadj", permutation.type = "none", GSassociation = TRUE,
      executation.info = TRUE, mc.cores = 1, GSs = NULL,  basic.cont = TRUE,
      GS = NULL, signature.method = "zscore",  ...){

 # controls
 if(basic.cont){
   correction.type <- basic.controls(nda = nda,
          correction.type = correction.type, sigSize = sigSize)

    if(missing(score) & correction.type %in% c("none", "Gsignature"))
          score = NULL
    if(missing(score) & !correction.type %in% c("none", "Gsignature"))
          stop("score must be specified")
    signature.method <- signmethod.control(signature.method)
    if(correction.type %in% c("NCsignatures","Gsignature") & signature.method=="zscore")
       adj.var.GS <- adjvarGS.control(adj.var.GS)

 }
 permutation.type <- perm.controls(permutation.type)

 # initialization
 eps <- 10e-8
 if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

 if(is.null(GSs))
  GSs <- GSnegControls(nda, correction.type = correction.type, score = score,
           seqquant = seqquant,  comp.RUV = comp.RUV, ...)

 if(is.null(GS))
     GS  <- GSnegControls(nda = nda, correction.type =  "Gsignature",
                          score=NULL, seqquant = seqquant,  ...)
 # diagnostics
 if(executation.info) cat("\n permutations diagnostic \n")
 ds1 <-  dim(nda)[2]
 perms <- vapply(seq_len(nperm), function(i) sample(seq_len(ds1)), integer(ds1))

 perm.out <- permutationsOut(nda, sigSize = sigSize, perms = perms,
               type = permutation.type, nite = nrs,
               mc.cores = mc.cores, GSassociation = GSassociation,
               GS = GSs[[1]], adj.var.GS = adj.var.GS,
               executation.info=executation.info, signature.method=signature.method)

 out <- list(perm.out =  perm.out, perms = perms, sigSize = sigSize, GSs = GS)
 class(out) <- c("hrunbiasedPermutations", "hrunbiasedDiagnostic")

 attr(out, "ev.info") <- nda$evn
 return(out)
}


##### hrunbiasedPCpower
hrunbiasedPCpower <- function(nda, positive.controls,
     correction.type = "NCsignatures", score,
     seqquant = c(0,seq(0.05,1,length.out=6)),
     alternative = "two.sided", bootstrapping = TRUE,  prop.pos.cont = 1,
     sigSize = 50,  id.size = 1, nrs = 5000, nrpcs = 500,
     comp.RUV = c(1,2), null.dist = "asymptotic",  adj.var.GS = "GSadj",
     adjusted.var.random = NA, adjusted.var.fixed= NA, executation.info = TRUE,
     mc.cores = 1, GSs = NULL, cdRS = NULL, basic.cont = TRUE, ...){

  # controls
  if(basic.cont){
   correction.type <- basic.controls(nda = nda,
          correction.type = correction.type, sigSize = sigSize)

    if(missing(score) & correction.type %in% c("none", "Gsignature"))
          score = NULL
    if(missing(score) & !correction.type %in% c("none", "Gsignature"))
          stop("score must be specified")
    signature.method <- signmethod.control(signature.method)
   if(correction.type %in% c("NCsignatures","Gsignature") & signature.method=="zscore")
       adj.var.GS <- adjvarGS.control(adj.var.GS)
  }
  if(missing(positive.controls))
    stop("positive.controls must be specified")
  if(any(!positive.controls %in% rownames(nda)))
    stop("positive.controls must be in rownames(nda)")

  mod.controls(nda, adjusted.var.random, adjusted.var.fixed, GSs)
  bc2 <- pow.controls(alternative, prop.pos.cont, null.dist, id.size, sigSize)
  alternative <- bc2$alternative
  null.dist <- bc2$null.dist

  # initialization
  pc.HR <- pc.pow <- list()
  eps <- 10e-8
  if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

  if(is.null(GSs))
    GSs <- GSnegControls(nda, correction.type =  correction.type,
             score=score, seqquant = seqquant,  comp.RUV = comp.RUV, ...)

  # diagnostics
  if(is.null(cdRS)){
   sigSize <- sigSize[id.size]
   cdRS <- hrunbiasedDensityRS(nda = nda, correction.type = correction.type,
            score = score, seqquant = seqquant, sigSize = sigSize,
            nrs = nrs, comp.RUV = comp.RUV, adj.var.GS = adj.var.GS,
            adjusted.var.random = adjusted.var.random,
            geneMean.geneSign = FALSE, adjusted.var.fixed = adjusted.var.fixed,
            executation.info = executation.info, mc.cores = mc.cores, ...)
    rs.HR    <- cdRS$hr.rs[[1]]
    randsign <- cdRS$randsign[[1]]
    sigSize  <- cdRS$sigSize[1]
   }
   else{
    rs.HR    <- cdRS$hr.rs[[id.size]]
    randsign <- cdRS$randsign[[id.size]]
    sigSize  <- cdRS$sigSize[id.size]
   }

   if(executation.info) cat("\n positive controls \n")
   positive.controls <- intersect(positive.controls, rownames(exprs(nda)))
   genes.not.poscontrol <- rownames(exprs(nda))[!rownames(exprs(nda))%in%
                               positive.controls]

   for(ppc.i in seq_len(length(prop.pos.cont))){
    ppc.x  <- prop.pos.cont[ppc.i]
    randsign.s <- lapply(seq_len(nrpcs), function(j)
        gs <- c(sample(positive.controls, round(sigSize*ppc.x),
                replace = bootstrapping), sample(genes.not.poscontrol,
                round(sigSize*(1-ppc.x)))))

        pc.HR[[ppc.i]] <- rsCoxmeByCorrections(nda, correction.type, score,
                    seqquant, randsign.s, adj.var.GS,
                    adjusted.var.random = adjusted.var.random,
                    adjusted.var.fixed = adjusted.var.fixed,
                    mc.cores = mc.cores, GSs = GSs,
                    executation.info = executation.info, ...)
       pc.pow[[ppc.i]] <- posControlsPower(rs.HR, pc.HR[[ppc.i]],
                  alternative = alternative, correction.type = correction.type,
                  seqquant = seqquant, test = null.dist)
   }
   names(pc.HR) <- names(pc.pow) <- paste0("prop.pos ", prop.pos.cont)
   out <- list(hr.rs = rs.HR, hr.pc = pc.HR, pvals.pc = pc.pow)
   class(out) <- c("hrunbiasedPCpower","hrunbiasedDiagnostic")

   attr(out,"null.dist") <-  null.dist
   if(correction.type %in% c("NCsignatures","RUV") &  length(seqquant) > 1){
      attr(out, "cor.nc.rs") <- corNCSvsRS(apply(t(exprs(nda)[rownames(nda)%in%
            positive.controls,]),1, eval("mean"), na.rm = TRUE), GSs, seqquant)
      attr(out, "means.rs") <-  vapply(rs.HR, function(x) mean(x[,4]), double(1))
   }
   return(out)
}


##### hrunbiasedSimulations
hrunbiasedSimulations <- function(nda, positive.controls, GScoef = 0,
     tscoef = 0, ninst.rs = 10, ninst.ts = 500, nrs = 500, lambdas = c(0.1,0.2),
     gammas = c(1.5,2),  mc.cores = 1, stat = TRUE,
     executation.info = TRUE, GSs = NULL, signature.method = "zscore",
     basic.cont = TRUE, sigSize = 50, correction.type = "none", ...){
 # controls
 if(missing(positive.controls))
    stop("positive.controls must be specified")
 if(any(!positive.controls %in% rownames(nda)))
    stop("positive.controls must be in rownames(nda)")
 signature.method <- signmethod.control2(signature.method)

 # initialization
 eps <- 10e-8
 if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

 randsign  <- lapply(seq_len(nrs), function(j) gs <- sample(rownames(exprs(nda)),sigSize))
 randsign.tar <- lapply(seq_len(nrs), function(j) gs <- sample(positive.controls,sigSize))

 ts <- scale(apply(exprs(nda)[positive.controls,],2,mean))
 if(is.null(GSs)) GS <- (GSnegControls(nda = nda,
        correction.type =  "Gsignature", score = NULL, seqquant = 0, ...)[[1]])
 else GS <- (GSs[[1]])

 if(signature.method == "zscore"){
     if(correction.type=="GSignature")
         rs <- scale(vapply(randsign, function(o) apply(exprs(nda)[o,],2,mean) - GS, double(length(GS))))
     else
         rs <- scale(vapply(randsign, function(o) apply(exprs(nda)[o,],2,mean),double(dim(nda)[2])))
     ts2 <- ts
     if(correction.type=="GSignature")
         rs.tar <- scale(vapply(randsign.tar, function(o) apply(exprs(nda)[o,],2,mean) - GS,double(length(GS))))
     else
         rs.tar <- scale(vapply(randsign.tar, function(o) apply(exprs(nda)[o,],2,mean),double(dim(nda)[2])))
 }
 if(signature.method =="gsva"){
     gsva.aux <- gsva(nda, randsign, method="gsva", parallel.sz =mc.cores)
     if(is(gsva.aux,"list")) gsva.aux <- exprs(gsva.aux$es.obs)
     rs <- scale(t(exprs(gsva.aux)))

     gsva.aux <- gsva(nda, list(positive.controls), method="gsva", parallel.sz =mc.cores)
     if(is(gsva.aux,"list")) gsva.aux <- exprs(gsva.aux$es.obs)
     ts2 <- scale(as.numeric(exprs(gsva.aux)))

     gsva.aux <- gsva(nda, randsign.tar, method="gsva", parallel.sz =mc.cores)
     if(is(gsva.aux,"list")) gsva.aux <- exprs(gsva.aux$es.obs)
     rs.tar <- scale(t(exprs(gsva.aux)))
 }

 GS <- scale(GS)
 covs <- as.data.frame(rs)
 covs2 <- data.frame(GS = GS, ts = ts)
 covs <- data.frame(covs, ts2)

 # diagnostic
 simu.info <- coefsSimulations(ts = ts2, GS = GS, covs = covs, covs2 = covs2,
              tscoef = tscoef, GScoef = GScoef, lambdas = lambdas,
              gammas = gammas, ninst.rs = ninst.rs, ninst.ts = ninst.ts,
              mc.cores = mc.cores, executation.info = executation.info,
              stat = stat, rs.tar = rs.tar)

 simu.info$stat <- stat
 simu.info$show.target <- ninst.ts>0
 out <- list(simu.info = simu.info)
 class(out) <- c("hrunbiasedSimulations","hrunbiasedDiagnostic")
 return(out)
}

##### hrunbiasedPCAcorrelation
hrunbiasedPCAcorrelation <- function(nda, nrs = 100, sigSize = 50,
                     mc.cores = 1, comp.PCA = c(1,2),
                     GS = NULL, executation.info = TRUE, signature.method = "zscore", ...){

  signature.method <- signmethod.control2(signature.method)
  # initialization
  eps <- 10e-8
  if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

  randsign <- list()
  for(K in seq_len(length(sigSize))){
      randsign[[K]] <- lapply(seq_len(nrs), function(j)
          gs <- sample(rownames(exprs(nda)), sigSize[K]))
  }
  xs <- t(exprs(nda))

  # diagnostic
  out <- pcaCorrFun(xs = xs, GS = GS, sigSize = sigSize,
          executation.info = executation.info, signature.method = signature.method,
          randsign = randsign, comp.PCA = comp.PCA, mc.cores = mc.cores)

  class(out) <- c("hrunbiasedPCAcorrelation", "hrunbiasedDiagnostic")
  attr(out, "comp.PCA") <- comp.PCA
  return(out)
}

##### print.hrunbiasedDiagnostic
print.hrunbiasedDiagnostic <- function(x,  ...){
  cat(paste("Diagnostics for Gene signature Cox proportional",
              "hazard models:"),  "\n")
  if(inherits(x, "hrunbiasedDensityRS"))
  {
    cat("  ",paste("density.rs"))
    cat("\t\t", paste0("plots available: density.rs, density.rs.size ",
       "and geneMean.geneSign"), "\n")
  }
  if(inherits(x, "hrunbiasedPermutations"))
  {
    cat("  ",paste("permutations"))
    cat("\t\t", paste0("plots available: perm.violin, perm.GSvsEvents ",
        "and perms.corr.GS"), "\n")
  }
  if(inherits(x, "hrunbiasedSimulations"))
  {
   cat("  ",paste("simulations"))
   cat("\t\t", paste0("plots available: simulations"), "\n")
  }
  if(inherits(x, "hrunbiasedPCpower"))
  {
   cat("  ",paste("power.pos.cont"))
   cat("\t", paste0("plots available: positive.cont.power, ",
       "corNCsig.power and bias.power"), "\n")
  }
  if(inherits(x, "hrunbiasedPCAcorrelation"))
  {
    cat("  ",paste("PCAcorrelation "))
    cat("\t", paste0("plots available: PCAcorrelation"), "\n")
  }
}


##### plot.hrunbiasedDiagnostic
plot.hrunbiasedDiagnostic <- function(x, diagnostic.plot, id.size = 1,
          show.stat = TRUE, legend.out = TRUE, perms.to.show = 1:10,
          seq.alpha = seq(0,.5, length.out=50), seq.alpha.pow = c(0.01,0.05),
          power.whplots = NULL, col1 = "gray", col2 = "orange",
          pointsq = FALSE, ...){

  if(inherits(x, "hrunbiasedDensityRS"))
  {
   if("density.rs" %in% diagnostic.plot)
      plotrsCoxmeByCorrections(x$hr.rs[[id.size]],
             correction.type = x$correction.type, show.stat = show.stat,
             legend.out = legend.out, seqquant = x$seqquant, ...)
   if("density.rs.size" %in% diagnostic.plot)
      rsCoxmeBySigSizes(x$hr.rs, x$sigSize, x$seqquant,
             correction.type = x$correction.type,  show.stat = show.stat,
             legend.out = legend.out, ...)
   if("geneMean.geneSign" %in% diagnostic.plot)
            plotGeneMeanHRvsSignatureHR(x$geneMeanHR.signatureHR)
  }
  if(inherits(x, "hrunbiasedPermutations"))
  {
   if(any(c("perm.violin", "perm.GSvsEvents", "perms.corr.GS")%in%
          diagnostic.plot))
       plotPermutationsOut(x$perm.out, plot.type = diagnostic.plot,
           perms.to.show = perms.to.show, legend.out = legend.out,
           perms = x$perms, GS = x$GS[[1]], ev.info = attr(x, "ev.info"),
           sigSize = x$sigSize, ...)
  }
  if(inherits(x, "hrunbiasedSimulations"))
  {
    if("simulations"%in%diagnostic.plot)
       plotcorrectionsSimulations(x$simu.info, col1 = col1,
            col2 = col2, legend.out = legend.out, ...)
  }
  if(inherits(x, "hrunbiasedPCpower"))
  {
    if("positive.cont.power" %in% diagnostic.plot)
       plotposControlsPower(x$pvals.pc, x$correction.type,
           seqquant = x$seqquant, seq.alpha = seq.alpha,
           test = attr(x, "null.dist"), legend.out = legend.out,
           power.whplots =  power.whplots,  ...)
    if("corNCsig.power" %in% diagnostic.plot)
       plotpowAndBias(means = attr(x, "cor.nc.rs") , pc.pow = x$pvals.pc,
           seq.alpha =  seq.alpha.pow, seqquant =  x$seqquant,
           test = attr(x, "null.dist"), legend.out = legend.out,
           meansPlot =FALSE,  power.whplots =  power.whplots, ...)
    if("bias.power" %in% diagnostic.plot)
       plotpowAndBias(means = attr(x, "means.rs"), pc.pow = x$pvals.pc,
           seq.alpha =  seq.alpha.pow, seqquant =  x$seqquant,
           test = attr(x, "null.dist"), legend.out = legend.out,
           meansPlot = TRUE, power.whplots =  power.whplots, ...)

  }
  if(inherits(x, "hrunbiasedPCAcorrelation"))
  {
    if("PCAcorrelation"%in%diagnostic.plot)
       hrunbiasedPCAcorrelationplot(x, col1 = col1, col2 = col2,
                       id.size = id.size, legend.out = legend.out, ...)

  }
}

