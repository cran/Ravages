burden.mlogit <- function(x, NullObject, genomic.region = x@snps$genomic.region, burden, maf.threshold = 0.5, get.effect.size = FALSE, alpha = 0.05, cores = 10){
  
  if (is.numeric(burden)) {
    if (!is.matrix(burden)) {
      stop("Score is not a matrix")
    }
    else {
      if (is.null(colnames(burden))) {
        colnames(burden) <- make.names(1:ncol(burden))
      }
    }
    #Check between number of individuals
    if(nrow(burden) != length(NullObject$group)) stop("Different number of individuals in 'burden' and 'NullObject'")
    score <- burden

  }else {
    if (!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
    genomic.region <- droplevels(genomic.region)
    if (missing(x)) stop("a bed.matrix 'x' is needed to compute the score")
  
    #Check between number of individuals
    if(nrow(x) != length(NullObject$group)) stop("Different number of individuals in 'x' and 'NullObject'")

    if (burden == "CAST") {
      score <- CAST(x, genomic.region, maf.threshold)
    }else if (burden == "WSS") {
      score <- WSS(x, genomic.region)
    }else {
      stop("'burden' should be \"CAST\", \"WSS\", or a matrix of pre-computed burdens")
    }
  }
  score <- as.data.frame(score)
  old.names <- colnames(score)
  names(score) <- make.names(names(score))
 
  alt.levels <- levels(NullObject$group)[levels(NullObject$group) != NullObject$ref.level]
  data.reg <- cbind(NullObject$data, score) ; rownames(data.reg) <- NULL
  data.reg <- dfidx(data.reg, varying = NULL, shape = "wide", choice = "ind.pheno")
  R <- do.call(rbind, mclapply(names(score), function(reg) run.mlogit.withNull(pheno = NullObject$group,  score = score, region = reg, ref.level = NullObject$ref.level, alt.levels = alt.levels, covar.toinclude = NullObject$covar.toinclude, data = data.reg, alpha = alpha, get.effect.size = get.effect.size, H0.LogLik = NullObject$H0.LogLik), mc.cores = cores))
  if (get.effect.size) colnames(R) <- c("p.value", "is.err", paste("OR", alt.levels, sep = "."), paste("l.lower", alt.levels, sep = "."), paste("l.upper", alt.levels, sep = "."))
  else colnames(R) <- c("p.value", "is.err")
  rownames(R) <- old.names
  return(as.data.frame(R))
}



run.mlogit.withNull <- function (pheno, score, region, ref.level, alt.levels, covar.toinclude, data, alpha, get.effect.size, H0.LogLik){
    if (is.null(covar.toinclude)) {
      my.formula <- Formula(as.formula(paste("ind.pheno ~ 0 |", region)))
    }else {
      my.formula <- Formula(as.formula(paste("ind.pheno ~ 0 |", region, " + ", covar.toinclude)))
    }
    fit <- tryCatch(mlogit(my.formula, data = data, reflevel = ref.level), error = identity, warning = identity)
    if (is(fit, "error")) {
        pval <- NA
        is.err <- 1
        OR.values <- data.frame(Estimate = rep(NA, nlevels(pheno) - 1), sd = rep(NA, nlevels(pheno) - 1))
        rownames(OR.values) <- paste(alt.levels, "region", sep = ":")
    }
    else {
      if (is.null(covar.toinclude)) {
        my.model <- summary(fit)
        pval <- as.numeric(my.model$lratio$p.value)
        if (get.effect.size == TRUE) OR.values <- my.model$CoefTable
      }else {
        my.model.H1 <- summary(fit)
        pval <- pchisq(-2 * H0.LogLik + 2 * as.numeric(my.model.H1$logLik), nlevels(pheno) - 1, lower.tail = FALSE)
        if (get.effect.size == TRUE) OR.values <- my.model.H1$CoefTable
      }
      is.err <- 0
    }
    quantile.alpha <- qnorm(alpha/2, lower.tail = FALSE)
    if (get.effect.size){
      OR.values.estimate <- OR.values[paste(region, alt.levels, sep = ":"), 1]
      OR.values.sd <- OR.values[paste(region, alt.levels, sep = ":"), 2]
      results <- c(pval, is.err, exp(OR.values.estimate), exp(OR.values.estimate - quantile.alpha * OR.values.sd), exp(OR.values.estimate + quantile.alpha * OR.values.sd))
    }else{ results <- c(pval, is.err)}
    
    #Cleaning temporary objects
    rm(score) ; rm(data) ; rm(fit) ; gc()
    return(results)
}
