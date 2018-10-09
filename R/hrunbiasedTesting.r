## Function hrunbiasedTesting
hrunbiasedTesting <- function(nda, gene.signature, alternative ="two.sided",
           test.type = c("asymptotic", "random.signatures", "random.genes"),
           correction.type = c("NCsignatures", "Gsignature", "RUV", "none"),
           score = NULL, seqquant = 0.25,  nrs = 5000,
           comp.RUV = c(1,2), adj.var.GS = "GSadj",
           adjusted.var.random = NA, adjusted.var.fixed = NA, mc.cores = 1,
           signature.method = "zscore", maxGenes = 5000, executation.info=TRUE, ...){

  eps <- 10e-8
  if(abs((sum((exprs(nda)*exprs(nda))) - prod(dim(nda)) + dim(nda)[1])) > eps)
    exprs(nda) <- t(scale(t(exprs(nda))))

  if(length(seqquant) != 1 & correction.type %in% c("NCsignatures", "RUV"))
  {
        warning("only the first element in seqquant is going to be used")
        seqquant <- seqquant[1]
  }
  GS <- GSnegControls(nda = nda, correction.type = correction.type,
           score = score, seqquant = seqquant, comp.RUV = comp.RUV, ...)[[1]]
  signature.method <- signmethod.control(signature.method)

  pvals.rs <- NA
  rsc.gs <- randomSigCoxmeSimple(nda, list(gene.signature),  sigSize = NULL,
              mc.cores = mc.cores, signature.method = signature.method,
              needScaling = FALSE, nite = 500, GS = GS,
              adjusted.var.random = adjusted.var.random,
              adjusted.var.fixed = adjusted.var.fixed, adj.var.GS = adj.var.GS,
              executation.info=executation.info)
  rsc.gs <- do.call(rbind, rsc.gs)
  rsc <- NULL
  rsc.by.gene.all <- NULL
  rsc.by.gene <- NULL
  if("random.signatures" %in% test.type)
  {
     sigSize <- length(gene.signature)
     randsign <- lapply(seq_len(nrs), function(j)
                   gs <- sample(rownames(exprs(nda)), sigSize))
     rsc <- do.call(rbind, randomSigCoxmeSimple(nda, randsign,  sigSize = NULL,
           mc.cores = mc.cores, signature.method = signature.method, needScaling = FALSE,
           nite = 500, GS = GS, adjusted.var.random = adjusted.var.random,
           adjusted.var.fixed = adjusted.var.fixed, adj.var.GS = adj.var.GS,
           executation.info=executation.info))

     if(alternative == "two.sided")
       pvals.rs <-  vapply( rsc.gs[,4], function(x) 2 * min(mean(x < rsc[,4]),
                      mean(x > rsc[,4])), double(1))
     if(alternative == "greater")
       pvals.rs <- vapply( rsc.gs[,4], function(x)  1 - mean(x > rsc[,4]), double(1))
     if(alternative == "less")
       pvals.rs <- vapply( rsc.gs[,4], function(x)  1 - mean(x < rsc[,4]), double(1))
  }
  rsc.by.gene <- do.call(rbind, randomSigCoxmeSimple(nda,
        as.list(gene.signature), sigSize = NULL, mc.cores = mc.cores,
        signature.method = signature.method, needScaling = FALSE, nite = 500, GS = GS,
        adjusted.var.random = adjusted.var.random,
        adjusted.var.fixed =  adjusted.var.fixed, adj.var.GS = adj.var.GS,
        executation.info = executation.info))
  rownames( rsc.by.gene) <- gene.signature

  if("random.genes" %in%  test.type){
     urs <- rownames(nda)#unique(unlist(randsign))
     urs <- sample(urs, min(length(urs), maxGenes))
     rsc.by.gene.all <- do.call(rbind, randomSigCoxmeSimple(nda, as.list(urs),
          sigSize = NULL, mc.cores = mc.cores, signature.method = signature.method,
          needScaling = FALSE, nite = 500, GS = GS,
          adjusted.var.random = adjusted.var.random,
          adjusted.var.fixed = adjusted.var.fixed, adj.var.GS = adj.var.GS,
          executation.info = executation.info))
  }

  out <- list(hr.ts = rsc.gs, hr.rs = rsc, hr.tg = rsc.by.gene,
            hr.rg = rsc.by.gene.all, pvals.rs = pvals.rs,
            test.type = test.type, gene.signature = gene.signature)
  class(out) <- "hrunbiasedTesting"
  return(out)
}

print.hrunbiasedTesting <- function(x, ...){
    cat("\t",paste("Cox model results for target signature"),  "\n")
    res <- cbind(x$hr.ts, x$pvals.rs)
    colnames(res)[5:6] <- c("pval.asy", "pval.rs")
    rownames(res) <- "ts"
    print(res)

    cat("\n\n")
    cat("\t", paste("random signature summary lHR"), "\n")
    print(summary(x$hr.rs[,1]))
    cat("\t",paste("random signature summary stat"), "\n")
    print(summary(x$hr.rs[,4]))
}

plot.hrunbiasedTesting <- function(x, legend.out = TRUE, show.pvalues = TRUE,
          show.stat = TRUE,  density.tg.minpoints = 10, ...){
  if(legend.out) par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=FALSE)
  else par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=FALSE)

  rgtest <- "random.genes" %in%  x$test.type
  rstest <- "random.signatures" %in%  x$test.type

  indic <- ifelse(show.stat,4,1)
  if(length(x$hr.tg[,indic]) >= density.tg.minpoints) ds2 <- density(x$hr.tg[,indic])
  else ds2 <- list(x=0,y=0)
  lhn <- ifelse(show.stat,"stat","lHR")

  if(rstest) ds1 <- density(x$hr.rs[,indic])
  else ds1 <- list(x=0,y=0)

  if(rgtest) ds4 <- density(x$hr.rg[,indic])
  else ds1 <- list(x=0,y=0)

  ylimm <- max(c(ds1$y, ds2$y, ds4$y, 0.2))
  rus <- runif(length(x$gene.signature), 0, ylimm/8)
  plot(x$hr.tg[,indic], rus, xlim = c(min(c(ds1$x,ds2$x,ds4$x,x$hr.ts[,indic])),
                max(c(ds1$x,ds2$x, ds4$x,x$hr.ts[,indic]))),
        ylim = c(0,ylimm ),  pch="*", col=2, type="p",...)
  lines(rep(x$hr.ts[,indic],2), c(0,ylimm/8), col=3, lwd=2, lty=1)

  if(rstest) lines(ds1,  lwd=2, lty=3)
  if(rgtest)  lines(ds4, col=4, lwd=2, lty=4)
  if(length(x$hr.tg[,indic])>=density.tg.minpoints) lines(ds2, col=2, lwd=2, lty=2)

  if(legend.out){
      par(xpd=TRUE)
      pos.lev <- c(".rs", ".rgene", ".ind.gene", ".ts")
      ltys <- c(3,4, ifelse(length(x$hr.tg[,indic])>density.tg.minpoints, 2, NA))
      cols <- c(1,4,2,3)
      pchs <- c(NA,NA,"*",NA)
      lev.out <- pos.lev[c(rstest,rgtest,TRUE,TRUE)]
      cols.out <- cols[c(rstest,rgtest,TRUE,TRUE)]
      ltys.out <- ltys[c(rstest,rgtest,TRUE,TRUE)]
      pchs.out <- pchs[c(rstest,rgtest,TRUE,TRUE)]

      if(!rgtest)  legend("topright",
         legend = paste0(lhn, lev.out), col=cols.out,
         lty= ltys.out, lwd = 2, pch=c(NA,"*",NA))
  else{
     coord <- par("usr")
     legend(x = coord[2] * 1.05, y = coord[4],
       legend = paste0(lhn, c(".rs", ".rgene", ".ind.gene", ".ts")),
       col = c(1,4,2,3), lty= c(3,4, ifelse(length(x$hr.tg[,indic])>density.tg.minpoints, 2, NA),1),
       lwd=c(2,2,2,2), pch = pchs.out, title = paste0(lhn, " types"), cex = 1)
    }
  }
  if(show.pvalues){
      if(rstest){
          leg <- paste0(c("pval.a = ", "pval.rs = "),
            format.pval(c(x$hr.ts[,5], x$pvals.rs),1))
          pchs <- c(20,21)
      }
      else{
          leg <- paste0(c("pval.a = "),
            format.pval(c(x$hr.ts[,5]),1))
          pchs <- c(20)
      }
      legend("topleft", legend= leg, pch= c(20),
      box.col="white", bty="n")
  }
}





