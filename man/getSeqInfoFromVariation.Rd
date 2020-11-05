\name{getSeqInfoFromVariation}
\alias{getSeqInfoFromVariation}
\title{Collects information about genomic context of sequence variants.}
\usage{
getSeqInfoFromVariation(referenceDnaStringSet, transcriptID, 
variation, ntWindow=20, transcriptTable,gene2transcript=gene2transcript)
}
\arguments{
\item{referenceDnaStringSet}{DNAStringset from the reference genome fasta file.}
\item{transcriptID}{Ensembl ID of the transcript of interest.}
\item{variation}{A sequence variation either refering to coding sequence or the genomic sequence (c.12A>T, or g.182284A>T).}
\item{ntWindow}{Numeric value defining the sequence surrounding of interest.}
\item{transcriptTable}{Table of transcrits and their exon coordinates and CDS coordinates.}
\item{gene2transcript}{Gene to transcript conversion table with the gene name in the first column and the gene ID in the second and the transcript ID in the third column.}
}
\description{
This function collects information about genomic context of sequence variants.
}
\examples{
#Defining exemplary input data
transcriptTable <- transCoord
transcriptID <- "pseudo_ENST00000650636"
variation <- "c.412C>G/p.(T89M)"
gene2transcript <- data.frame(gene_name = "Example_gene",
gene_ID = "pseudo_ENSG00000147099", transcriptID = "pseudo_ENST00000650636")

results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID,
variation, ntWindow=20, transcriptTable, gene2transcript=gene2transcript)

#Using a predefined gene to transcript conversion
transcriptID <- "Example_gene"
results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID,
variation, ntWindow=20, transcriptTable, gene2transcript=gene2transcript)

}
\value{
List of informations about the entered variation.
}
