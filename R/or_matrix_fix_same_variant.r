OR.matrix.same.variant <- function (n.variants, OR.del, OR.pro = 1/OR.del, p.causal, prob.pro, maf, maf.threshold){
    if (length(OR.del) != length(OR.pro))
        stop("Dimensions mismatch")
    #Causal variants only among variants with maf < threshold
    w <- which(maf <= maf.threshold)
    if(length(w)==0) stop(paste0("No variants with MAF lower than ", maf.threshold))
    nb.causal <- length(w)*p.causal
    #NULL vector if !condition
    v.causal <- if(nb.causal>0) sample(w, nb.causal) 
    v.protect <- if(prob.pro>0) sample(v.causal, nb.causal*prob.pro)

    OR.tot <- matrix(rep(1, n.variants*nrow(OR.del)), nrow=nrow(OR.del))    

    #If NULL vector: OR.tot won't be change
    OR.tot[,v.causal] <- OR.del[,v.causal]
    OR.tot[,v.protect] <- OR.pro[,v.protect]
    return(list(OR = OR.tot, causal = v.causal))
}
