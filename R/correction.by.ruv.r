## Function correction.by.ruv
correction.by.ruv <- function(Y, X, covar = NULL, ncomp, nmin = 200,
                       which.cases = NA, steps = 2){
    if(is.na(which.cases)[1]){
            low.var <- apply(Y, 2, mad)
            ctl <- rep(FALSE, dim(Y)[2])
            ctl[order(low.var)[seq_len(nmin)]] <-  TRUE
    }
    else{
            ctl <- rep(FALSE, dim(Y)[2])
            ctl[which.cases] <-  TRUE
    }
    ncomp.v <- ncomp
    ncomp <- max(ncomp.v)

    if(!is.null(covar)){
      form1 <- formula(paste0("Y[,1]~", paste(covar, collapse ="+")))
      modm <- model.matrix(form1, data=X)[,-1]

      if(steps == 2) svafit <- RUV2(Y, modm, ctl = ctl, k=ncomp)
      if(steps == 4) svafit <- RUV4(Y, modm, ctl = ctl, k=ncomp)
    }
    else
      svafit <- fa.neg.cont(Y, ctl = ctl, k=ncomp)

    return(list(W=svafit$W))
}

fa.neg.cont <- function(Y,ctl, k){
  Yc = Y[, ctl]
  W = svd(Yc %*% t(Yc))$u[, seq_len(k), drop = FALSE]
  out <- list(W=W)
  return(out)
}

