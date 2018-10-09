## posControlsPower function
posControlsPower <- function(rs.HR, pc.HR, alternative = "two.sided",
       seqquant =  seq(0.05,1,length.out=6), correction.type,
       test = c("asymptotic","random.signatures"), ...){
   ln <- length(rs.HR)
   pvals.nda.cor.1.rs <- pvals.nda.cor.1.asy <- NULL

   if("random.signatures" %in% test)
   {
     if(alternative == "two.sided")
        pvals.nda.cor.1.rs <- sapply(seq_len(ln), function(K)
           sapply(pc.HR[[K]][,4], function(x) 2 * min(mean(x < rs.HR[[K]][,4]),
                   mean(x > rs.HR[[K]][,4]))))
     if(alternative == "greater")
         pvals.nda.cor.1.rs <- sapply(seq_len(ln), function(K)
            sapply(pc.HR[[K]][,4], function(x) 1 - mean(x > rs.HR[[K]][,4])))
     if(alternative == "less")
         pvals.nda.cor.1.rs <- sapply(seq_len(ln), function(K)
            sapply(pc.HR[[K]][,4], function(x) 1 - mean(x < rs.HR[[K]][,4])))
   }
   if("asymptotic" %in% test)
    pvals.nda.cor.1.asy <- sapply(seq_len(ln), function(K) pc.HR[[K]][,5])

   return(list(pvals.nda.cor.1.asy = pvals.nda.cor.1.asy,
               pvals.nda.cor.1.rs =  pvals.nda.cor.1.rs))
}

plotposControlsPower <- function(hr.pc,  correction.type,
   seqquant =  seq(0.05,1,length.out=6), seq.alpha = seq(0,.5, length.out=50),
   test = c("asymptotic", "random.signatures"), legend.out = TRUE,
   power.whplots = NULL, ...){

  if(is.null( power.whplots)) power.whplots <- seq_len(length(hr.pc))
   if("random.signatures" %in% test)
   {
    ln <- dim(hr.pc[[1]][[2]])[2]
    if(legend.out) par(mar=c(5.1, 5.1, 4.1, 9.5), xpd=TRUE)
    else par(mar=c(5.1, 5.1, 4.1, 5.1), xpd=TRUE)

    for(ss in power.whplots){
        plot(seq.alpha, sapply(seq.alpha, function(o)
             mean(hr.pc[[ss]]$pvals.nda.cor.1.rs[,1] < o)), type="l",
             col = grey.colors(ln)[1], lty = 1,
             ylab = "proportion rejections pc.rs",  xlab = "alpha",
             ylim=c(0,1), ...)
        for(K in seq_len(ln))
          lines(seq.alpha, sapply(seq.alpha, function(o)
            mean(hr.pc[[ss]]$pvals.nda.cor.1.rs[,K] < o)),
             col = grey.colors(ln)[K], lty=K, ...)
        if(correction.type %in% c("NCsignatures", "Gsignature", "RUV")
             & legend.out )
        {
           coord <- par("usr")
           legend(x = coord[2] * 1.05, y = coord[4],
             legend = paste0("neg.cont.",seqquant), col = grey.colors(ln),
             lty = seq_len(ln), lwd = 2)
        }
        lines(c(min(seq.alpha), max(seq.alpha)), c(min(seq.alpha),
                max(seq.alpha)), col = 5,lty = 2)
     }
    }
   if("asymptotic" %in% test)
   {
    ln <- dim(hr.pc[[1]][[1]])[2]
    if(legend.out) par(mar=c(5.1, 5.1, 4.1, 9.5), xpd=TRUE)
    else par(mar=c(5.1, 5.1, 4.1, 5.1), xpd=TRUE)
    for(ss in power.whplots){
      plot(seq.alpha, sapply(seq.alpha, function(o)
        mean(hr.pc[[ss]]$pvals.nda.cor.1.asy[,1] < o)),  type="l",
        col = grey.colors(ln)[1], lty=1,
        ylab = "proportion rejections pc.asy", xlab = "alpha",
        ylim=c(0,1), ...)
      for(K in seq_len(ln)) lines(seq.alpha, sapply(seq.alpha, function(o)
          mean(hr.pc[[ss]]$pvals.nda.cor.1.asy[,K] < o)),
          col=grey.colors(ln)[K], lty=K, ...)
      if(correction.type %in% c("NCsignatures", "Gsignature", "RUV") &
          legend.out)
      {
                coord <- par("usr")
                legend(x = coord[2] * 1.05, y = coord[4],
                     legend = paste0(round(seqquant,3)), title = "neg.cont",
                     col = grey.colors(ln), lty=seq_len(ln), lwd=2)
      }
      lines(c(min(seq.alpha),max(seq.alpha)), c(min(seq.alpha),max(seq.alpha)),
            col=5,lty=2)
     }
   }
}
