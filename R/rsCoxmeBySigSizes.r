## Function rsCoxmeBySigSizes
rsCoxmeBySigSizes <- function(rs.HR.size, sigSize, seqquant, correction.type,
                       show.stat = TRUE, legend.out = TRUE, ...){
    ln <- length(rs.HR.size[[1]])
    whstat <- ifelse(show.stat, 4, 1)

      op <- par(oma = c(0,0,0,1) + 0.1,
         mar = c(1,1,.5,2) + 0.1, xpd = FALSE)
    if(legend.out){
      layout(rbind(t( vapply(2*c(0:(ln-1)), function(l) l+c(1,2), double(2))),
        rep(ln*2+1,2)), widths=c(.5,3),
        heights = c(rep(1,ln), min(1,ln * 0.1)))
    }
    else
    {
     layout(rbind(t(vapply(2*c(0:(ln-1)), function(l) l+c(1,2), double(2) )),
        rep(ln*2+1,2)), widths=c(.5,3),
        heights = c(rep(1,ln), min(1,ln * 0.05)))
    }
    xlims1 <- c(min( unlist(lapply(rs.HR.size, function(x) min(unlist(lapply(x,
             function(y) min(y[,whstat]))))))), max( unlist(lapply(rs.HR.size,
             function(x) max(unlist(lapply(x,function(y) max(y[,whstat]))))))))

    lz <- length(rs.HR.size)
    greys <- gray.colors(lz)
    for(k in seq_len(ln)){
        plot(0,0, axes = FALSE, frame.plot = FALSE, main="", xlim= xlims1*1.2,
            bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xlab="",ylab="")
        if(ln>1){
            text(x = 0, y = 0, paste("q_",round(seqquant[k],3)), cex = 1,
            col = "black", srt=90)
        }
        aux <- lapply(rs.HR.size, function(x) x[[k]])
        ylims <- c(0,max(sapply(aux, function(x) max(density(x[,whstat])$y))))
      if(k==ln){
       plot(density(aux[[1]][,whstat]), col=greys[1], lty=1, ylim=ylims,
             axes=TRUE,xlab="dsd", frame.plot=FALSE,  xlim= xlims1*1.2,
             main = " ", ...)
      }

     else{
       plot(density(aux[[1]][,whstat]), col=greys[1], lty=1, ylim=ylims,
             axes=TRUE,xaxt="n",xlab="", frame.plot=FALSE,  xlim= xlims1*1.2,
             main = " ", ...)
      }
      for(i in seq_len(length(aux)))  lines(density(aux[[i]][,whstat]),
            col=greys[i], lty=i, ...)
        lines(c(0,0), ylims, col="grey")
    }
    if(legend.out){
        plot(0,0, axes=FALSE,frame.plot=FALSE, main="", xlim = xlims1*1.2,
             bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n', xlab="", ylab = "")
        legend("center", legend=paste("size = ", sigSize), col=greys, lty=seq_len(lz),
               lwd=1, bty="n", ncol = ifelse(lz>5,  round(lz/2), lz), cex=1)
    }
    par(mfrow=c(1,1))
}
