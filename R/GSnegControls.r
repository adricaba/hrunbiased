## function GSnegControls
GSnegControls <- function(nda, correction.type, score, seqquant,
                          comp.RUV = c(1,2), covar.RUV = NULL, ...){
 if(correction.type %in%  c("Gsignature", "none")) ln <- 1
 else  ln <- length(seqquant)

 GS.all <- lapply(seq_len(ln), function(K){
   if(correction.type %in% c("none") |
     (!(correction.type %in% c("Gsignature","SVA")) & seqquant[K] == 0)) GS <- NULL
   else{
    if(correction.type %in% c("NCsignatures"))
    {
     genes.GS  <- which(score <= quantile(score, seqquant[K]))
     GS        <- apply(exprs(nda)[genes.GS,], 2, mean, na.rm=TRUE)
    }
    if(correction.type %in% "RUV")
    {
     genes.GS  <- which(score <= quantile(score, seqquant[K]))
     corr.ruv  <- correction.by.ruv(t(exprs(nda)), pData(nda),
                       covar = covar.RUV, ncomp = comp.RUV, nmin = 200,
                       which.cases = genes.GS, steps = 2)
     GS  <- corr.ruv$W
    }
    if(correction.type %in% "SVA")
    {
        GS  <- correction.by.sva(t(exprs(nda)), pData(nda),
                       covar = covar.RUV, ncomp = comp.RUV[1])
    }
    if(correction.type %in% c("Gsignature"))
     GS        <- apply(exprs(nda), 2, mean, na.rm=TRUE)
  }
  GS
 })
 return(GS.all)
}

require(dSVA)
correction.by.sva <- function(Y, X, covar = NULL, ncomp = 0){
    if(!is.null(covar)){
        form1 <- formula(paste0("Y[,1]~", paste(covar, collapse ="+")))
        modm <- model.matrix(form1, data = X)[,-1]
        svafit <- dSVA(Y, modm, ncomp = ncomp)
        coefs <- svafit$Z
    }
    else{
        coefs <- prcomp(Y)$x[,seq_len(max(ncomp,1))]
    }
    return(coefs)
}



