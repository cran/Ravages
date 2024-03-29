SKAT.theoretical <- function(x, NullObject, genomic.region = x@snps$genomic.region, 
                     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                     estimation.pvalue = "kurtosis", cores = 10, debug = FALSE){


  #Check between number of individuals
  if(nrow(x) != length(NullObject$group)) stop("Different number of individuals in 'x' and 'NullObject'")

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  #Replace (if needed...) genomic region in x@snps -- so it will be modified by select.snps below
  x@snps$genomic.region <- genomic.region
  
  if(any(table(genomic.region)==1)) 
    stop("All 'genomic.region' sould contain at least 2 variants, please use 'filter.rare.variants()' to filter the bed matrix")

  x@snps$weights <- weights
  x <- select.snps(x, x@snps$maf <= maf.threshold & x@snps$maf > 0)
  
  group <- NullObject$group
  
  lev <- levels(group) 
  #Number of individuals by group
  n <- as.vector(table(group))
  
  #Matrix pi
  Pi <- NullObject$Pi.data
  
  #Matrix of Pi in all groups
  P1 <- NullObject$P1

  #Matrix of (YY - Pi^) with YY indicatrice variable in each group 
  YY <- sapply(lev, function(l) as.numeric(group == l))
  ymp = YY -  Pi
  
  #P-value for all regions
  res.allregions <- do.call(rbind, mclapply(levels(genomic.region), function(REG) get.parameters.pvalue.theoretical(x, RR = REG, P1 = P1, ymp = ymp, n = n, estimation.pvalue = estimation.pvalue), mc.cores=cores))
  res.final <- as.data.frame(res.allregions)
  colnames(res.final) <- colnames(res.allregions)
  rownames(res.final) <- levels(genomic.region)
  if(debug)
    res.final
  else
    res.final[, c("stat", "p.value")]
} 
  

  
