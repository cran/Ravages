#Function for the simulations
#Sample causal variants, compute haplotypes burden
#and thresholds corresponding to the desired prevalence
#p.causal=P(causal variant) ; p.protect=P(protective variant | causal variant)

rbm.haplos.thresholds <- function(haplos, weights = c("SKAT", "constant"), max.maf.causal = 0.01, 
             p.causal = 0.5, p.protect = 0, h2, prev, normal.approx = TRUE, size, replicates, rep.by.causal, verbose = TRUE) {

  weights <- match.arg(weights)
  
  if(any(colMeans(haplos)==0 | colMeans(haplos)==1)){
     warning("Some variants have a maf equal to 0 and won't be kept")
     haplos <- haplos[,which(colMeans(haplos)>0 & colMeans(haplos)<1)]
  }
  
  if (weights=="SKAT") {
        weights <- -0.4*log10(colMeans(haplos))
  }else{
    if (weights=="constant") {
          weights <- rep(1, ncol(haplos))
    }
  }

  
  if( (replicates %% rep.by.causal) != 0 ) 
    stop("replicates should be a multiple of rep.by.causal")

  if(length(h2) != length(prev) | length(h2) != length(size) | length(prev) != length(size)) stop("h2 and prev should have same size as size")
  if(length(max.maf.causal)==1){
    max.maf.causal <- rep(max.maf.causal, length(h2))
  }else{
    if(length(max.maf.causal) != length(h2)) stop("max.maf.causal should be of size 1 or of same size as h2 and prev")
  }
  if(length(p.causal)==1){
    p.causal <- rep(p.causal, length(h2))
  }else{
    if(length(p.causal) != length(h2)) stop("p.causal should be of size 1 or of same size as h2 and prev")
  }
  if(length(p.protect)==1){
    p.protect <- rep(p.protect, length(h2))
  }else{
    if(length(p.protect) != length(h2)) stop("p.protect should be of size 1 or of same size as h2 and prev")
  }
  
  if(verbose) cat(length(h2), " groups of individuals will be simulated \n")
    
  
  x <- new.bed.matrix(nb_inds=sum(size), nb_snps=ncol(haplos)*replicates)
  #For causal variants 
  x@snps$Causal <- ""
 
  for(i in 1:length(seq(1,replicates, by=rep.by.causal))){
    burdens.causal <- mapply(get.causal.burden, weights=list(weights), haplos = list(haplos), max.maf.causal, p.causal, p.protect, h2, SIMPLIFY=FALSE)
    burdens <- lapply(burdens.causal, function(z) z$burdens) 
  
    ##Computation of thresholds to respect 'prev' 
    ##normal.approx = TRUE : ok if small h2
    if(normal.approx) {
      s <- qnorm(prev, lower.tail = FALSE)
    } else {  
      BRD <- mapply(sample, x = burdens, size = 1e6, replace = TRUE) + mapply(sample, x = burdens, size = 1e6, replace = TRUE) + mapply(rnorm, n = 1e6, sd = sqrt(1-h2))
      s <- sapply(1:length(prev), function(z) quantile(BRD[,z], 1 - prev[z]))
    }

    .Call('rbm_haplos_thresholds_filling', PACKAGE = 'Ravages', x@bed,  haplos, burdens, sd = sqrt(1-h2), thr1 = s, thr2 = rep(Inf, length(size)), size, repNumber = i-1, reps=rep.by.causal)
    #add flag on causal variants in cases
    for(j in which(prev<1)){
      x@snps[as.vector(outer(((i-1)*rep.by.causal):(i*rep.by.causal-1)*ncol(haplos), burdens.causal[[j]]$causal, FUN = "+")), "Causal"] <- paste0(x@snps[as.vector(outer(((i-1)*rep.by.causal):(i*rep.by.causal-1)*ncol(haplos), burdens.causal[[j]]$causal, FUN = "+")), "Causal"], "Gpe", j)
    }
  }

  x@ped$pheno <- factor(rep.int( 1:length(size) - 1, size))
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = ncol(haplos)) )
  x@snps$Causal[which(x@snps$Causal=="")] <- "NonCausal"
 
  if( is.null(colnames(haplos)) )
    x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  else
    x@snps$id <- paste( x@snps$genomic.region, colnames(haplos), sep = "_")

  x <- set.stats(x, verbose = FALSE)
  x
}



get.causal.burden <- function(weights, haplos, maf.threshold, p.causal, p.protect, h2){
  weights[colMeans(haplos) == 0 | colMeans(haplos) > maf.threshold ] <- 0
  w <- which(weights > 0) # parmi les SNPs potentiellement causaux
  
  if(length(w) == 0) 
    stop("There are not enough positively weighted variants")

  nb.causal <- round(p.causal*length(w))
  nb.pro <- round(nb.causal * p.protect)
  nb.del <- nb.causal - nb.pro
  
  w.causal <- sample(w, nb.causal)
  directions <- numeric(ncol(haplos))
  directions[w.causal] <- c(rep(-1,nb.pro), rep(1,nb.del))

  ##Computation of burden (genetic component of liability)
  ##Based on the vector 'weights'
  we <- directions * weights
  burdens <- haplos %*% we
  #Standardisation of burdens
  burdens <- burdens - mean(burdens)

  ##Respect h2: adjust burdens' variance of the two haplotypes 
  tau <- as.numeric(2*var(burdens))
  burdens <- burdens / sqrt(tau) * sqrt(h2)
  #To generate the liability, a gaussian variable of standard deviation sqrt(1-h2) should be added to the burdens
  return(list(burdens = burdens, causal = w.causal))
}
