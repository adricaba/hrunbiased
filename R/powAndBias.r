## function corNCSvsRS
corNCSvsRS <- function(Gsmr, GSs, seqquant){
           d1 <- unlist(lapply(GSs,function(x) dim(x)[2]))[1]
           if(is.null(d1)) d1 <- 1
           if(d1==1){
               means  <- as.numeric(cor(do.call(cbind,GSs), Gsmr))
               if(seqquant[1] == 0) means <- c(0, means)
           }
           else{
              means <- do.call(cbind,
                        lapply(GSs, function(x) if(!is.null(x)) cor(x, Gsmr)))
               if(seqquant[1] == 0) means <- cbind(rep(0,d1), means)
           }
       return(means)
}


plotpowAndBias <- function(means,  pc.pow, seq.alpha = 0.05,
                           seqquant =  seq(0.05,1,length.out=6),
                           test = "asymptotic", legend.out = TRUE,
                           meansPlot = TRUE, power.whplots = NULL, ...){

 if(legend.out) par(mar=c(5.1, 5.1, 4.1, 9.9), xpd=TRUE)
 else par(mar=c(5.1, 5.1, 4.1, 5.1), xpd=TRUE)
 if(is.null(power.whplots)) power.whplots <- 1:length(pc.pow)

 for(ss in power.whplots){
   if(test == "asymptotic")
     POWERS <- sapply(seq.alpha, function(o)
                 apply(pc.pow[[ss]][["pvals.nda.cor.1.asy"]], 2,
                 function(x) mean(x < o)))
   if(test =="random.signatures")
     POWERS <- sapply(seq.alpha, function(o)
                 apply(pc.pow[[ss]][["pvals.nda.cor.1.rs"]], 2,
                 function(x) mean(x < o)))

   ylims <- c(-max(abs(means)),max(abs(means)))
   if(is.matrix(means)){
     for(i in 1:dim(means)[1]){
        if(i ==1)
          plot(seqquant, means[1,], pch = 16, axes=FALSE,  ylim =ylims,
               xlab="", ylab="", type="b",col="red", ...)
         else
          lines(seqquant, means[i,], pch=16+i, lty= i,type="b",col="red")
     }
   }
   else
     plot(seqquant, means, pch=16, axes=FALSE,  ylim =ylims, xlab="",
          ylab="", type="b",col="red", ...)
   lines(seqquant,rep(0,length(seqquant)), col=3,lty=3)
   axis(2, col="red",col.axis="red",las=1,...)
   if(meansPlot)
     mtext("log HR stat mean",side=2,line=2.5, col="red", ...)
   else
     mtext("correlation (nc.sig, pc.sig)",side=2,line=2.5, col="red", ...)
   box()

   par(new=TRUE)

   plot(seqquant, POWERS[,1], ylim=c(0,1), pch=15,  xlab="",
         ylab="", axes=FALSE, type="b", col="grey", lty=2, ...)
   for(io in 1:length(seq.alpha))
       lines(seqquant, POWERS[,io], pch=c(15-io+1),   xlab="", ylab="",
             xaxt='n', yaxt="n", type="b", col="grey", lty=io+1)
   coord <- par("usr")
   if(legend.out)
       legend(x = coord[2] * 1.15, y = coord[4],
         legend = c(paste0(round(seq.alpha,3))), title = "rej.levels",
         col="grey", lty = c(2:(length(seq.alpha)+1)),
         pch= c(15-1:length(seq.alpha)+1))

    ## a little farther out (line=4) to make room for labels
    axis(4,  col="grey",col.axis="grey",las=1,...)
    mtext("proportion rejections pc.rs",side=4,col="grey",line=2.5, ...)

    axis(1, col="black",las=1,...)
    mtext("quantile neg.controls",side=1,col="black",line=2.5, ...)
}
      #      par(mfrow=c(1,1))
 }


