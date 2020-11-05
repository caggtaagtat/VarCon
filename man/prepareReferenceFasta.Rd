\name{prepareReferenceFasta}
\alias{prepareReferenceFasta}
\title{Imports Fasta file from filepath.}
\usage{
prepareReferenceFasta(filepath)
}
\arguments{
\item{filepath}{R conform filepath to the fasta file of the reference genome to use.}
}
\description{
This function imports Fasta file of the reference genome into R enviroment as DNAStringset.
}
\value{
Creates new DNAStringSet from the object stated by the entered filepath.
}
\examples{
 ## Loading exemplary DNAStringSet
 filepath <- system.file("extdata", "fastaEx.fa", package="Biostrings")
 referenceDnaStringSet <- prepareReferenceFasta(filepath)

}
