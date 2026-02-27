\name{adjustedCADD.annotation}
\alias{adjustedCADD.annotation}
\title{
  SNVs and Indels annotation with adjusted CADD scores
}
\description{
  Annotate SNVs and Indels with the adjusted CADD scores (CADD PHRED scores for coding, regulatory and intergenic regions)
}
\usage{
adjustedCADD.annotation(x, SNVs.scores = NULL, indels.scores = NULL,
			cores = 10, verbose = T, path.data, 
                        build = c("b37", "b38"))
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{SNVs.scores}{ A dataframe with the columns 'chr', 'pos', 'A1', 'A2' and 'adjCADD' containing the ADJUSTED CADD scores of the SNVs (Optional, useful to gain in computation time if the adjusted CADD scores of variants in the study are available)}
  \item{indels.scores}{ A dataframe with the columns 'chr', 'pos', 'A1', 'A2' and 'PHRED_1.4' (or 'PHRED_1.6' if build 38) containing the CADD PHREDv1.4 (of PHREDv1.6) scores of the indels - Compulsory if indels are present in \code{x}}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
  \item{path.data}{ The repository where data for RAVA-FIRST are or will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/}
  \item{build}{The build of the data, either "b37" or "b38". The CADD Regions in the corresponding build will be considered}
}
\value{
  The bed matrix x with adjusted CADD scores in \code{adjCADD}.
}
\details{
  bedtools need to be installed on the system to run this function if indels are present are if scores are not provided for SNVs.  

  This function calls \code{adjustedCADD.annotation.SNVs} and \code{adjustedCADD.annotation.indels}. See the help of those two functions for more details.  
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{adjustedCADD.annotation.SNVs}, \link{adjustedCADD.annotation.indels}, \link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}} 
}
\examples{
####Build37
#Import 1000Genome data from region around LCT gene
#x37 <- read.bed.matrix( system.file("extdata", "LCT.EUR.b37.bed", package="Ravages") )

#Annotate variants with adjusted CADD score
#x37 <- adjustedCADD.annotation(x37, build = "b37")

####Build38
#Import 1000Genome data from region around LCT gene
#x38 <- read.bed.matrix( system.file("extdata", "LCT.EUR.b38.bed", package="Ravages") )

#Annotate variants with adjusted CADD score
#x38 <- adjustedCADD.annotation(x38, build = "b38")

}
