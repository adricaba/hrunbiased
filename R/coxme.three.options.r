## function coxme.three.options
coxme.three.options <- function(nda, perms, randsign,
  type = c("none","Gsignature", "norm01"),  GS = NULL, mc.cores=1,
  executation.info =TRUE, adj.var.GS = "GSadj",
  signature.method = "zscore"){

    if(!is(randsign,"matrix")) randsign <- t(as.matrix(randsign))
    nite <- dim(randsign)[2]
    nperms <- dim(perms)[2]
    enda2 <- exprs(nda)
    ds <- pData(nda)

    if(type[1] == "norm01"){
        enda2 <- t(apply(enda2,1, sample))
        exprs(nda) <- enda2
    }
    if(type[1] == "Gsignature")
        Gsmr.GS <- GS
    else
        Gsmr.GS <- NULL
    if(signature.method=="gsva"|signature.method=="plage")
            {
                randsign <- lapply(seq_len(nite), function(i) randsign[,i])
                if(length(randsign[[1]])>1){
                    gsva.aux <- gsva(nda, randsign, method=signature.method, parallel.sz = mc.cores )
                    if(is(gsva.aux,"list")) gsva.aux <- exprs(gsva.aux$es.obs)
                    if(is(gsva.aux,"ExpressionSet"))  gsva.aux <- exprs(gsva.aux)
                    if(signature.method=="plage"){
                        cors <- as.numeric(sign(cor(t(gsva.aux), apply(enda2,2,mean))))
                        gsva.aux <- (gsva.aux) * cors
                    }
                }
                else{
                    gsva.aux <- enda2[as.character(randsign),]
                }
            }


    new.coefs.NC <- array(0, dim = c(nperms, nite, 5))
    if(executation.info)
       pb <- txtProgressBar(min = 0, max = nperms, style = 3)
    for (ss in seq_len(nperms)){
        if(executation.info) setTxtProgressBar(pb, ss)
        perm <- perms[,ss]
        new.coefs <- mclapply(seq_len(nite), function(k){
            if(signature.method=="zscore"){
                gs.1 <- randsign[,k]
                if(length(gs.1)==1) Gsmr <- enda2[rownames(nda)%in%gs.1,]
                else Gsmr <- apply(t(enda2[rownames(nda)%in%gs.1,]),1, mean,
                               na.rm = TRUE)
            }
            if(signature.method=="gsva"|signature.method=="plage"){
              Gsmr <- gsva.aux[k,]
            }
            if(adj.var.GS == "GScor"  & !is.null(Gsmr.GS))
            {
                Gsmr <- Gsmr - Gsmr.GS
                Gsmr.GS <- NULL
            }
            ds$evn <- as.numeric(ds$evn)[perm]
            ds$tev <- as.numeric(ds$tev)[perm]
            if(type[1] =="Gsignature" & !is.null(Gsmr.GS) & adj.var.GS == "GSlmcor")
                ds$gs <- scale(lm(Gsmr ~ Gsmr.GS)$resid)
            else
                ds$gs <- scale(Gsmr)
            if(type[1] =="Gsignature" & !is.null(Gsmr.GS) & adj.var.GS == "GSadj")
                m <- coxph(Surv(tev, evn) ~  gs + Gsmr.GS, data=ds)
            else
                m <- coxph(Surv(tev, evn) ~  gs, data=ds)
            summary(m)[ "coefficients"][[1]]["gs",]
        }, mc.cores=mc.cores)
        new.coefs.NC[ss,,] <- do.call(rbind,new.coefs)
    }
    return(new.coefs.NC)
}
