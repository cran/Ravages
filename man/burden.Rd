\name{burden}
\alias{burden}
\title{
  Linear, logistic or multinomial regression on a genetic score
}
\description{
  Performs burden tests on categorical or continuous phenotypes
}
\usage{
burden(x, NullObject, genomic.region = x@snps$genomic.region, burden,
       maf.threshold = 0.5, get.effect.size = FALSE, alpha = 0.05, cores = 10,
       verbose = TRUE)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor containing the genomic region of each SNP, \code{x@snps$genomic.region} by default, only needed if \code{burden}="CAST" or \code{burden}="WSS"}
  \item{burden}{"CAST" or "WSS" to directly compute the CAST or the WSS genetic score, or a matrix with one row per individual and one column per \code{genomic.region} if another genetic score is wanted.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score. Set at 0.5 by default}
  \item{get.effect.size}{ TRUE/FALSE: whether to return effect sizes of the tested \code{genomic.region} (OR for categorical phenotypes, betas for continuous phenotypes) }
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{cores}{ How many cores to use, set at 10 by default.}
  \item{verbose}{ Whether to display information about the function actions }
}
\details{
  This function will return results from the regression of the phenotype on the genetic score for each genomic region.

  If only two groups of individuals are present, a classical logistic regression is performed.
  If more than two groups of individuals are present, a non-ordinal multinomial regression is performed,
  comparing each group of individuals to the reference group indicated by the argument \code{ref.level} in \code{NullObject.parameters}.
  The choice of the reference group won't affect the p-values, but only the Odds Ratios.
  In both types of regression, the p-value is estimated using the Likelihood Ratio test and the function \code{burden.mlogit}.  

  If the phenotype is continuous, a linear regression is performed using the function \code{burden.continuous}.  
  
  The type of phenotype is determined from \code{NullObject$pheno.type}.  

  If another genetic score than CAST or WSS is wanted, a matrix with one row per individual and one column per \code{genomic.region}
  containing this score should be given to \code{burden}. In this situation, no bed matrix \code{x} is needed.
}
\value{
  A dataframe with one row per genomic region and at least two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{NullObject$pheno.type = "categorical"} and \code{get.OR.value=TRUE}, additional columns are present:
  \item{OR/beta}{ The OR/beta value(s) associated to the regression. For categorical phenotypes, if there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR/beta}
  \item{l.upper}{ The upper bound of the confidence interval of each OR/beta}
}


\seealso{ \code{\link{NullObject.parameters}}, \code{\link{burden.continuous}}, \code{\link{burden.mlogit}},  \code{\link{CAST}}, \code{\link{WSS}}, \code{\link{burden.weighted.matrix}} }
\references{ Bocher O, et al. DOI: 10.1002/gepi.22210. \emph{Rare variant association testing for multicategory phenotype.} Genet.Epidemiol. 2019;43:646–656.}
\examples{
#Import data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

#Group variants within known genes
x <- set.genomic.region(x)

#Filter of rare variants: only non-monomorphic variants with
#a MAF lower than 2.5%
#keeping only genomic regions with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

#run null model, using the 1000Genome population as "outcome"
x1.H0 <- NullObject.parameters(pheno = x1@ped$pop, ref.level = "CEU",
                               RVAT = "burden", pheno.type = "categorical")

#run burden test WSS
burden(x1, NullObject = x1.H0, burden = "WSS", get.effect.size=TRUE, cores = 1)


}

