## ----message=FALSE, warning=FALSE, echo = FALSE-------------------------------
library("knitr")
require("Ravages")

## -----------------------------------------------------------------------------
# GRR calculated using the same formula as in SKAT,
# with values in the second group of cases being twice the values from the first one

GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, n.case.groups = 2,
                      GRR.multiplicative.factor=2, select.gene = "R1")
GRR.del[,1:5]

# Calculation of genotype frequencies in the two groups of cases and the controls group
# The previous GRR matrix is used with a multiplicative model of the disease
# The prevalence in each group of cases is 0.001

geno.freq.groups <- genotypic.freq(genes.maf = Kryukov, select.gene="R1", 
                                   GRR.het = GRR.del, prev = c(0.001, 0.001),
                                   genetic.model = "multiplicative")
str(geno.freq.groups)

# frequencies of the homozygous alternative genotype in the different groups
geno.freq.groups$freq.homo.alt[,1:5]

## -----------------------------------------------------------------------------
#MAF calculation for the first five variants
geno.freq.groups$freq.homo.alt[,1:5] + 0.5*geno.freq.groups$freq.het[,1:5]

## -----------------------------------------------------------------------------
x <- rbm.GRR(genes.maf = Kryukov, size = c(1000, 500, 500), replicates = 5,
             prev = c(0.001, 0.001), GRR.matrix.del = GRR.del, p.causal = 0.5, 
             p.protect = 0, same.variant = FALSE,
             genetic.model = "multiplicative", select.gene = "R1")
x
table(x@ped$pheno)
table(x@snps$genomic.region)

## -----------------------------------------------------------------------------
#Load LCT dataset for haplotype matrix
data(LCT.haplotypes)

#Simulation based on haplotypes frequencies
#for the variants in the LCT gene in the EUR population
LCT.gene.hap <- LCT.hap[which(LCT.sample$super.population=="EUR"),
                        which(LCT.snps$pos>=136545410 & LCT.snps$pos<=136594750)]
#Individuals from EUR
LCT.sample.EUR <- subset(LCT.sample, super.population=="EUR")
#Matrix of haplotype frequencies
LCT.freqs <- sapply(unique(LCT.sample.EUR$population), function(z) 
                    ifelse(LCT.sample.EUR$population==z, 
                           1/table(LCT.sample.EUR$population)[z], 0))
#Simulation of genetic data for five groups of 50 individuals
x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size=rep(50,5), replicates=5)

#Simulation of 100 controls, and two groups of 50 cases with 30 causal variants
#and with the second group having half h2 and twice the prevalence
#compared to the first one
#5 replicates are performed and causal variants are sampled once
x <- rbm.haplos.thresholds(haplos=LCT.gene.hap, p.causal = 0.3, h2=c(0.01,0.01,0.02),
                           prev=c(1,0.01,0.005), size=c(100, 50, 50), replicates=5,
                           rep.by.causal = 5)

## -----------------------------------------------------------------------------
#CAST theoretical power
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, n.case.groups = 2,
                      GRR.multiplicative.factor=2, select.gene = "R1")
rbm.GRR.power(genes.maf = Kryukov, size = c(1000, 500, 500),
              prev = c(0.001, 0.001), GRR.matrix.del = GRR.del,
              p.causal = 0.5, p.protect = 0, select.gene="R1",
              same.variant = FALSE, genetic.model = "multiplicative",
                        power.type="theoretical")

#CAST power based on simulations on haplotypes (20 replicates)
rbm.haplos.power(haplos=LCT.gene.hap, freqs=LCT.freqs, size=rep(50,5), 
                 replicates=20, rep.by.causal = 10,  RVAT = "CAST", cores = 1)

