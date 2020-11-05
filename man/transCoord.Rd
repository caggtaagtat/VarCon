\docType{data}
\name{transCoord}
\alias{transCoord}
\title{Small table as exemplary transcript table with exon coordinates}
\format{data frame
\describe{
 \item{Gene.stable.ID}{Ensembl gene ID}
 \item{Transcript.stable.ID}{Ensembl Transcript ID}
 \item{Strand}{Strand of the feature}
 \item{Exon.region.start..bp.}{Smalles coordinate of the exon end coordinates of a specific exon}
 \item{Exon.region.end..bp.}{Largest coordinate of the exon end coordinates of a specific exon} 
 \item{cDNA.coding.start}{Start of the coding sequence}
 \item{cDNA.coding.end}{End of the coding sequence}
 \item{CDS.start}{Covered coding nucleotides start}
 \item{CDS.end}{Covered coding nucleotides end}
 \item{Exon.rank.in.transcript}{Rank of the exon within the respective transcript}
 \item{Exon.stable.ID}{Ensembl exon ID}
 \item{Chromosome.scaffold.name}{Name of the chromosome}
 
}}
\usage{
transCoord
}
\description{
Small table as exemplary transcript table with exon coordinates.
}
\examples{
\donttest{
 transCoord
}
}
\keyword{datasets}


