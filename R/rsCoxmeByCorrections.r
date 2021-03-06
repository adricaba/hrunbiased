
rsCoxmeByCorrections <- function(nda, correction.type, score, seqquant,
      randsign, adj.var.GS = "GSadj", comp.RUV = c(1,2), mc.cores = 1,
      adjusted.var.random = NA, adjusted.var.fixed= NA, GSs = NULL,
      executation.info=TRUE, signature.method = "zscore", ...){

  if(correction.type %in%  c("Gsignature", "none", "SVA")) ln <- 1
  else  ln <- length(seqquant)
  if(is.null(GSs))
    GSs <- GSnegControls(nda=nda, correction.type = correction.type,
           score=score, seqquant=seqquant,   comp.RUV = comp.RUV , ...)

  rcs.nda.cor.1 <- lapply(seq_len(ln), function(K){
    rcs <- randomSigCoxmeSimple(nda, randsign = randsign,  sigSize = NULL,
            mc.cores =  mc.cores, signature.method = signature.method, needScaling = FALSE,
            nite=500, GS = GSs[[K]],
            adjusted.var.random =  adjusted.var.random,
            adj.var.GS = adj.var.GS,  adjusted.var.fixed = adjusted.var.fixed,
            executation.info=executation.info)
    do.call(rbind,rcs)
  })

 return(rcs.nda.cor.1)
}


plotrsCoxmeByCorrections <- function(rcs.nda.cor.1,
         correction.type = "NCsignatures", show.stat = TRUE, legend.out= TRUE,
         seqquant = NA, ...){
     whstat <- ifelse(show.stat, 4, 1)
     if(legend.out) par(mar=c(5.1, 5.1, 4.1, 7.5), xpd= FALSE)
     else  par(xpd= FALSE)

     LIMS <- vapply(rcs.nda.cor.1, function(x){
         aux <- density(x[,whstat])
         c(min(aux$x),max(aux$x), max(aux$y))
     }, double(3))
     ln <- length(rcs.nda.cor.1)
     for( i in seq_len(ln)){
         if(i==1)
            plot(density(rcs.nda.cor.1[[1]][,whstat]), xlim=c(min(LIMS[1,]),
                         max(LIMS[2,])), ylim = c(0,max(LIMS[3,])), ...)
         else
            lines(density(rcs.nda.cor.1[[i]][,whstat]), col=grey.colors(ln)[i],
                  lty=i, ...)
     }
     coord <- par("usr")
     lines(c(0,0),coord[3:4], col=3,lty=2)

     if(legend.out){
      par(xpd=TRUE)
      if(correction.type %in% c("NCsignatures", "Gsignature", "RUV"))
       legend(x = coord[2] * 1.05, y = coord[4],
              legend = paste0(round(seqquant,3)), title = "neg.cont",
              col = grey.colors(ln), lty=seq_len(ln), lwd = 2)
       par(xpd=FALSE)
      }
     #par(mfrow=c(1,1))
 }
