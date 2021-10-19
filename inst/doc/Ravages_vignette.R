## ----message=FALSE, warning=FALSE, echo = FALSE-------------------------------
library("knitr")
require("Ravages")

## -----------------------------------------------------------------------------
# Import data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)
# Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

# Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

# Group variants within know genes by extending their positions
# 500bp upstream and downstream
# (the function uses build 37 unless told otherwise)
x <- set.genomic.region(x, flank.width=500)

## -----------------------------------------------------------------------------
# a quick look at the result
table(x@snps$genomic.region, useNA = "ifany")

# Filter variants with maf in the entire sample lower than 1%
# And keep only genomic region with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.01, min.nb.snps = 10)
table(x1@snps$genomic.region, useNA="ifany")

# run burden test CAST, using the 1000Genome population as "outcome"
# Null model for CAST
x1.H0.burden <- NullObject.parameters(x1@ped$pop, ref.level = "CEU", 
                                      RVAT = "burden", pheno.type = "categorial")
burden(x1, NullObject = x1.H0.burden, burden = "CAST", cores = 1)

# run SKAT, using the 1000Genome population as "outcome"
# COnstruct null model for SKAT, then run test with only a few permutations
x1.H0.SKAT <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorial")
SKAT(x1, x1.H0.SKAT, params.sampling=list(perm.target = 10, perm.max = 500))

# run a similar analysis but using the RAVA-FIRST approach with WSS
# RAVA-FIRST(x1, filter = "whole", maf.threshold = 0.01, min.nb.snps = 10,
#           burden = TRUE, x1.H0.burden, SKAT = F)


## ---- eval = F----------------------------------------------------------------
#  # Example bed matrix with 4 variants
#  x.ex <- as.bed.matrix(x=matrix(0, ncol=4, nrow=10),
#                        bim=data.frame(chr=1:4, id=paste("rs", 1:4, sep=""),
#                                       dist = rep(0,4), pos=c(150,150,200,250),
#                                       A1=rep("A", 4), A2=rep("T", 4)))
#  
#  # Example genes dataframe
#  genes.ex <- data.frame(Chr=c(1,1,3,4), Start=c(10,110,190,220), End=c(170,180,250,260),
#                         Gene_Name=factor(letters[1:4]))
#  
#  # Attribute genomic regions without splitting the variants
#  # attributed to multiple genomic regions
#  x.ex <- set.genomic.region(x.ex, regions = genes.ex, split = FALSE)
#  x.ex@snps$genomic.region
#  
#  # Split genomic regions
#  x.ex.split <- bed.matrix.split.genomic.region(x.ex, split.pattern = ",")
#  x.ex.split@snps$genomic.region

## -----------------------------------------------------------------------------
# Calculation of the genetic score with a maf threshold of 1%
CAST.score <- CAST(x = x1, genomic.region = x1@snps$genomic.region, maf.threshold = 0.01)
head(CAST.score)

## -----------------------------------------------------------------------------
WSS.score <- WSS(x = x1, genomic.region = x1@snps$genomic.region)
head(WSS.score)

## -----------------------------------------------------------------------------
Sum.score <- burden.weighted.matrix(x = x1, weights = rep(1, ncol(x1)))
head(Sum.score)

## ---- eval = F----------------------------------------------------------------
#  # Null model
#  x1.H0 <- NullObject.parameters(x1@ped$pop, ref.level = "CEU",
#                                 RVAT = "burden", pheno.type = "categorial")
#  # WSS
#  burden(x = x1, NullObject = x1.H0, burden ="WSS",
#        alpha=0.05, get.effect.size=TRUE, cores = 1)
#  
#  # Sex + a simulated variable as covariates
#  sex <- x1@ped$sex
#  u <- runif(nrow(x1))
#  covar <- cbind(sex, u)
#  # Null model with the covariate "sex"
#  x1.H0.covar <- NullObject.parameters(x1@ped$pop, ref.level = "CEU",
#                                       RVAT = "burden", pheno.type = "categorial",
#                                       data = covar, formula = ~ sex)
#  
#  # Regression with the covariate "sex" without OR values
#  # Using the score matrix WSS computed previously
#  burden(NullObject = x1.H0.covar, burden=WSS.score, cores = 1)

## ---- eval = F----------------------------------------------------------------
#  # Random continuous phenotype
#  set.seed(1) ; pheno1 <- rnorm(nrow(x1))
#  # Null model
#  x1.H0.continuous <- NullObject.parameters(pheno1, RVAT = "burden",
#                                            pheno.type = "continuous")
#  # Test CAST
#  burden(x1, NullObject = x1.H0.continuous, burden = "CAST", cores = 1)

## -----------------------------------------------------------------------------
# *** Functionally-informed WSS analysis ***
# Attribution of variants to regions and subregions
x2 <- set.genomic.region.subregion(x, regions = genes.b37,
                                  subregions = subregions.LCT)
# Burden test
burden.subscores(x2, x1.H0.burden, cores = 1)

## ---- eval = F----------------------------------------------------------------
#  # Null model
#  x1.null <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorial")
#  # Permutations because no covariates
#  SKAT(x1, x1.null, get.moments = "permutations", debug = TRUE,
#       params.sampling = list(perm.target = 100, perm.max =5e4))
#  # Theoretical on 1 core
#  SKAT(x1, x1.null, get.moments = "theoretical", debug = TRUE, cores = 1)

## ---- eval = F----------------------------------------------------------------
#  # Random continuous phenotype
#  set.seed(1) ; pheno1 <- rnorm(nrow(x1))
#  # Null Model with covariates
#  x1.H0.c <- NullObject.parameters(pheno1, RVAT = "SKAT", pheno.type = "continuous",
#                                   data = covar)
#  # Run SKAT
#  SKAT(x1, x1.H0.c)

## ---- eval = F----------------------------------------------------------------
#  # Annotation of variants with adjusted CADD scores
#  x <- adjustedCADD.annotation(x)
#  # Attribution of CADD regions
#  x.CADDregions <- set.CADDregions(x)

## ---- eval = FALSE------------------------------------------------------------
#  # Keep only CADD regions with 2 variants and variants with a MAF greater than 1%
#  # and with an adjusted CADD score greater than the median
#  x.median <- filter.adjustedCADD(x.CADDregions, maf.threshold = 0.01, min.nb.snps = 2)

## ---- eval = FALSE------------------------------------------------------------
#  # Functionally-informed WSS analysis
#  x.burden <- burden.subscores(x.median, x1.H0.burden, burden.function = WSS)
#  # SKAT analysis
#  x.SKAT <- SKAT(x.median, x1.H0.SKAT)

## ---- eval = F----------------------------------------------------------------
#  # *** RAVA-FIRST analysis with functionally-informed WSS and SKAT
#  #Burden parameters
#  burden.parameters <- list(get.effect.size = T, burden.function = WSS)
#  # Chromosome by chromosome
#  res.bychr <- vector("list", 22)
#  for(chr in 1:22){
#    x <- read.bed.matrix(paste0(path_file, prefix_vcf, chr, ".vcf.gz"))
#    res.bychr[[chr]] <- RAVA.FIRST(x, H0.burden = H0.burden, H0.SKAT = H0.SKAT,
#                                   burden.parameters = burden.parameters)
#  }

## -----------------------------------------------------------------------------
# Selection of each group of individuals
CEU <- select.inds(x1, pop=="CEU")
CEU
FIN <- select.inds(x1, pop=="FIN")
FIN
GBR <- select.inds(x1, pop=="GBR")
GBR

# Combine in one file:
CEU.FIN.GBR <- rbind(CEU, FIN, GBR)
CEU.FIN.GBR

