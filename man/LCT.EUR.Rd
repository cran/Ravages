\name{LCT.EUR}
\docType{data}
\alias{LCT.EUR.b37}
\alias{LCT.EUR.b37.fam}
\alias{LCT.EUR.b37.bed}
\alias{LCT.EUR.b38}
\alias{LCT.EUR.b38.fam}
\alias{LCT.EUR.b38.bed}

\title{LCT genotypes matrix}

\description{
  These data contain the genotype matrix from the 1000 Genomes EUR ancestry individuals in the Lactase gene (positions retrieved from GeneCards.  
  Information about individuals is present in \code{LCT.EUR.b37/b38.fam}, including the subpopulation, in a format needed to generate a bedmatrix.  
  Two set of files are available: one set in build 37, and one set in build 38.  
  Data have been anonymized by randomly sampling the genotype positions.
}

\format{
  \describe{
    \item{\code{LCT.EUR.b37/b38.bed}}{ Binary genotype file}
    \item{\code{LCT.EUR.b37/b38.bim}}{ Variant information file}
    \item{\code{LCT.EUR.b37/b38.fam}}{ Sample information file}
  }
}

\source{ Data were obtained from the 1000 Genomes Project. }

\references{ McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632 }

\examples{
#Import data in a bed matrix
x <- read.bed.matrix( system.file("extdata", "LCT.EUR.b37.bed", package="Ravages") )

}
