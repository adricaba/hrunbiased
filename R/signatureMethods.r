## Function signatureMethods
signatureMethods <- function(nda,  signature.method = "zscore", randsign,  mc.cores = 1, needScaling = FALSE){

        if(needScaling) enda2 <- t(scale(t(exprs(nda))))
        else enda2 <- exprs(nda)
        exprs(nda) <- enda2

        if(signature.method=="gsva"|signature.method=="plage")
            {
                if(length(randsign[[1]])>1){
                    gsva.aux <- gsva(nda, randsign, method=signature.method, parallel.sz = mc.cores )
                    if(is(gsva.aux,"list")) out <- exprs(gsva.aux$es.obs)
                    if(is(gsva.aux,"ExpressionSet")) out <- exprs(gsva.aux)
                    if(signature.method=="plage"){
                        cors <- as.numeric(sign(cor(t(out), apply(enda2,2,mean))))
                        out <- (out) * cors
                    }
                }
                else{
                   out <- enda2[as.character(randsign),]
                }
            }
        if(signature.method=="zscore")
         {
             if(length(randsign[[1]])>1)
                 out <- do.call(rbind,mclapply(randsign, function(o) apply(enda2[o,],2,mean), mc.cores= mc.cores))
             else
                 out <- enda2[as.character(randsign),]
         }
       return(out)
    }







