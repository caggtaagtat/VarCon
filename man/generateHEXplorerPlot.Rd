\name{generateHEXplorerPlot}
\alias{generateHEXplorerPlot}
\title{Generates plot with HZEI values and splice site strengths from a list holding information about an SNV.}
\usage{
generateHEXplorerPlot(variationInfoList, ntWindow)
}
\arguments{
\item{variationInfoList}{Output from the \code{getSeqInfoFromVariation} function.}
\item{ntWindow}{Numeric value defining the sequence surrounding of interest.}
}
\description{
This function generates a plot depicting the HZEI score changes and changes in the HBS or MaxEntScan score, from a sequence variation.
}
\examples{
#Defining exemplary input data
transcriptTable <- transCoord    # Using pseudo transcript table
transcriptID <- "pseudo_ENST00000650636"     # Using pseudo transcript 
variation <- "c.412C>G/p.(T89M)"
ntWindow <- 20
gene2transcript <- data.frame(gene_name = "Example_gene", 
gene_ID = "pseudo_ENSG00000147099", transcript_ID = "pseudo_ENST00000650636")

results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID, variation, ntWindow=ntWindow, transcriptTable, gene2transcript)

generateHEXplorerPlot(results)

}
\value{
Plot stating the HZEI values per nt and splice site strength with and without the SNV.
}
