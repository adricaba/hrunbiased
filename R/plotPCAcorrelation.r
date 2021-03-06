## Function hrunbiasedPCAcorrelationplot
pcaCorrFun <- function(xs, GS = NULL, sigSize, randsign, comp.PCA,
                 mc.cores=1,executation.info = TRUE, signature.method = "zscore"){

  pca1 <- prcomp(xs)
  if(is.null(GS)) GS <- apply(xs,1,mean)
  else GS <- GS[[1]]

  if(executation.info){
      cat("\n PCA correlations \n")
      pb <- txtProgressBar(min = 0,
                           max = length(sigSize) * 2, style = 3)
  }
  rs.cor <- mclapply(seq_len(length(sigSize)), function(ii){
    if(executation.info) setTxtProgressBar(pb, ii)
    y <- randsign[[ii]]
    if(length(y[[1]])>1){
       rs <- sapply(y, function(w)
             apply(xs[,w],1,mean))
    }
    else
       rs <- sapply(y, function(w) xs[,w])

      cor(rs,pca1$x[,comp.PCA])
  },mc.cores = mc.cores)

  if(signature.method=="zscore"){
     rs.gs.cor <- mclapply(seq_len(length(sigSize)), function(ii){
      if(executation.info)
          setTxtProgressBar(pb, ii +length(sigSize))
      y <- randsign[[ii]]
      if(length(y[[1]]) > 1){
        rs <- sapply(y, function(w)
             lm(apply(xs[,w],1,mean) ~ GS)$resid)
      }
      else
        rs <- sapply(y, function(w) lm(xs[,w] ~ GS)$resid)

      cor(rs,pca1$x[,comp.PCA])
    },mc.cores = mc.cores)
  }
  if(signature.method == "gsva")
  {
          rs.gs.cor <- lapply(seq_len(length(randsign)), function(j){
              gsva.aux <- gsva(xs, randsign[[j]], method="gsva")
               if(is(gsva.aux,"list")) gsva.aux <- exprs(gsva.aux$es.obs)
              cor(t(exprs(gsva.aux)), pca1$x[,comp.PCA])
          })
  }
  gs.cor <- cor(GS, pca1$x[,comp.PCA])

  out <- list(rs.cor = rs.cor, rs.gs.cor = rs.gs.cor,
              gs.cor = gs.cor)
  return(out)
}


hrunbiasedPCAcorrelationplot <- function(x, col1 = "gray", col2 = "orange",
                       id.size = 1, legend.out = TRUE,...){
    plot(0,0,xlim=c(-1,1), ylim=c(-1,1),  pch="",
         xlab= paste0("PC", attr(x,"comp.PCA")[1]),
         ylab=paste0("PC", attr(x,"comp.PCA")[2]),...)
    arrows(0, 0, x$rs.cor[[id.size]][,1], x$rs.cor[[id.size]][,2],
           col= col1, length = 0.1, lty=2)
    points(x$rs.cor[[id.size]][,1], x$rs.cor[[id.size]][,2],
           col=col1, pch = 16)

    arrows(0, 0, x$rs.gs.cor[[id.size]][,1], x$rs.gs.cor[[id.size]][,2],
           col= col2, length = 0.1, lty=2)
    points(x$rs.gs.cor[[id.size]][,1], x$rs.gs.cor[[id.size]][,2],
           col=col2, pch = 15)

    arrows(0, 0, x$gs.cor[1], x$gs.cor[2], col = 1, length = 0.1,
           lty=1, lwd = 2)
    points(x$gs.cor[1], x$gs.cor[2], col = 1,  pch = 17)

    if(legend.out)
        legend("topleft", legend = c("non.adj.sign.","adj.sign.", "GS"),
           col=c(col1,col2,"black"), lty=c(2,3,1), lwd=c(1,1,2),
           cex=1.3, pch = c(16, 15,17), bty = "n")
}


