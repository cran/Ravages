\name{adjustedCADD.annotation.indels}
\alias{adjustedCADD.annotation.indels}
\title{
  Indels annotation with adjusted CADD scores
}
\description{
  Annotate Indels with the adjusted CADD scores (CADD PHRED scores for coding, regulatory and intergenic regions)
}
\usage{
adjustedCADD.annotation.indels(x, variant.scores = NULL, 
			       cores = 10, verbose = T, path.data, 
                               build = c("b37", "b38"))
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{variant.scores}{ A dataframe with the columns 'chr', 'pos', 'A1', 'A2' and 'PHRED_1.4' (or 'PHRED_1.6 if build 38) containing the CADD PHREDv1.4 (or PHREDv1.6) scores of the indels}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
  \item{path.data}{ The repository where data for RAVA-FIRST are or will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/}
  \item{build}{The build of the data, either "b37" or "b38". The CADD Regions in the corresponding build will be considered}
}
\value{
  The bed matrix x with adjusted CADD scores in \code{adjCADD}.
}
\details{
  bedtools need to be installed on the system to run this function.  

  Indels are directly annotated with the adjusted CADD scores in the function using the file "AdjustedCADD_v1.4_202204_indels.tsv.gz" or "AdjustedCADD_v1.6_202602_indels.tsv.gz" downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/ in the repository of the package Ravages.  

  The adjusted CADD scores in "AdjustedCADD_v1.4_202204_indels.tsv.gz" have been computed using a set of 48M indels already annotated in the CADD website. The adjusted CADD scores in "AdjustedCADD_v1.6_202602_indels.tsv.gz" have been computed on the indels observed in GnomAD v1.3.2. If indels not present in this set are to be annotated, they will be given the same adjusted score as the indel with the nearest PHRED score provided in \code{variant.scores} which should contain the chromosome ('chr'), position ('pos'), reference allele ('A1'), alternative allele ('A2') and PHRED CADD scores v1.4 ('PHRED_1.4').  

  Those adjusted scores are used in the \code{RAVA.FIRST()} pipeline to filter rare variants.  
  
  As this function can take time when a large number of SNVs are present, it is recommended to use this function chromosome by chromosome for large datasets or to fitler the bed matrix before the annotation.
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{adjustedCADD.annotation}, \link{adjustedCADD.annotation.SNVs}, \link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}} 
}
