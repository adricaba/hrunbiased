############ function plotcorrectionsSimulations
coefsSimulations <- function(ts, GS, covs, covs2, tscoef, GScoef, stat,
   lambdas, gammas, ninst.rs, ninst.ts, mc.cores = 1, executation.info = TRUE, rs.tar){
 ds.info <- ds.info.gs <- ds.ts <- ds.ts.gs  <- list()

 # diagnostics
 if(executation.info) cat("\n simulations diagnostic \n")
 if(executation.info )
   pb <- txtProgressBar(min = 0, max = ninst.rs, style = 3)

 si <- ifelse(stat, 4, 1)
 for(KAS in seq_len(ninst.rs)){
  if(executation.info )  setTxtProgressBar(pb, KAS)
  ssd <- simsurv(lambdas = lambdas[1], gammas = gammas[1], x = covs2,
                 betas = c(ts = tscoef, GS = GScoef))

  tfollow <- rweibull(dim(covs2)[1], gammas[2], lambdas[2]^(-1/gammas[2]))
  evn <- ifelse(ssd[,2] - tfollow<0, 1, 0)
  tev <-  evn * ssd[,2] + (1-evn) * tfollow
  ssd$eventtime <- tev
  ssd$status <- evn

  ds <- data.frame(evn = ssd$status, tev = ssd$eventtime)
  ds$ts <- ts
  ds$GS <- GS

  rs.hr <-   do.call(cbind,mclapply(seq_len(dim(covs)[2]), function(i){
    gs <- scale(covs[,i])
    m <- coxph(formula(paste0("Surv(tev, evn) ~  gs")), data=ds)
    summary(m)["coefficients"][[1]]["gs",]
  },  mc.cores= mc.cores))
  ds.info[[KAS]] <- rs.hr[si,seq_len(dim(covs)[2]-1)] #density(rs.hr[si,1:(dim(covs)[2]-1)])

  rs.hr.gs <-   do.call(cbind,mclapply(seq_len(dim(covs)[2]), function(i){
    gs <- scale(covs[,i])
    m <- coxph(formula(paste0("Surv(tev, evn) ~ GS + gs")), data=ds)
    summary(m)["coefficients"][[1]]["gs",]
  },mc.cores=mc.cores))
  ds.info.gs[[KAS]] <- rs.hr.gs[si,seq_len(dim(covs)[2]-1)]#density(rs.hr.gs[si,1:(dim(covs)[2]-1)])

  rs.hr <-   do.call(cbind,mclapply(seq_len(dim(rs.tar)[2]), function(i){
    gs <- scale(rs.tar[,i])
    m <- coxph(formula(paste0("Surv(tev, evn) ~  gs")), data=ds)
    summary(m)["coefficients"][[1]]["gs",]
  },  mc.cores= mc.cores))
  ds.ts[[KAS]] <- rs.hr[si,seq_len(dim(covs)[2]-1)] #density(rs.hr[si,1:(dim(covs)[2]-1)])

  rs.hr.gs <-   do.call(cbind,mclapply(seq_len(dim(rs.tar)[2]), function(i){
    gs <- scale(rs.tar[,i])
    m <- coxph(formula(paste0("Surv(tev, evn) ~ GS + gs")), data=ds)
    summary(m)["coefficients"][[1]]["gs",]
  },mc.cores=mc.cores))
  ds.ts.gs[[KAS]] <- rs.hr.gs[si,seq_len(dim(covs)[2]-1)]#density(rs.hr.gs[si,1:(dim(covs)[2]-1)])

 }
 out <- list(ds.info = ds.info, ds.info.gs = ds.info.gs,
             ds.ts = ds.ts, ds.ts.gs = ds.ts.gs, tscoef = tscoef)
 return(out)

}

plotcorrectionsSimulations <- function(simu.info,
                    col1 = "gray", col2 = "orange", legend.out = TRUE,  ...){

     xlims.a <- max(max(abs(unlist(simu.info$ds.info))),max(abs(unlist(simu.info$ds.info.gs))))
     xlims <- c(-xlims.a,xlims.a)
     ylims.a <- max(max(unlist(vapply(simu.info$ds.info, function(y) max(density(y)$y), double(1)) )),
                    max(unlist(vapply(simu.info$ds.info.gs, function(y) max(density(y)$y), double(1)))))
     ylims <- c(0,ylims.a)
     plot(density(simu.info$ds.info[[1]]), col = col1, ylim = ylims, xlim = xlims ,...)
     for(i in 2:length(simu.info$ds.info))
         lines(density(simu.info$ds.info[[i]]), col = col1, lty=i, ...)

     for(i in seq_len(length(simu.info$ds.info.gs)))
         lines(density(simu.info$ds.info.gs[[i]]), col = col2, lty=i, ...)
#     curve(dnorm(x,0,1),from=-4,to=4,add=TRUE, col=2,lwd=2)
     abline(v=0,col=3,lty=2)
     if(legend.out) legend("topleft", c("rs (no adj)","rs (GS adj)"),
                            col=c(col1,col2), lwd = 2, lty=c(1,2))

}
#   if(simu.info$stat) curve(dnorm(x,0,1),
#        from=-4, to=4, add=TRUE, col="blue", lty=6, lwd=2)
