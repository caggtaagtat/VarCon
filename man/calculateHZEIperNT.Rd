\name{calculateHZEIperNT}
\alias{calculateHZEIperNT}
\title{Generates table with HZEI scores per nucleotide of a sequence.}
\usage{
calculateHZEIperNT(seq)
}
\arguments{
\item{seq}{Nucleotide sequence longer than 11nt and only containing bases "A", "G", "C", "T".}
}
\description{
This function generates a table with HZEI scores per index nucleotide.
}
\value{
Dataframe with HZEI value per index position.
}
\examples{
calculateHZEIperNT("TTCCAAACGAACTTTTGTAGGGA")

}
