\name{getMaxEntInfo}
\alias{getMaxEntInfo}
\title{Generates table with MaxEntScan scores per potential SA position.}
\usage{
getMaxEntInfo(seq)
}
\arguments{
\item{seq}{Nucleotide sequence longer than 22nt and only containing bases "A", "G", "C", "T".}
}
\description{
This function generates a table with MaxEntScan scores per potential SA position.
}
\examples{
getMaxEntInfo("TTCCAAACGAACTTTTGTAGGGA")

}
\value{
Dataframe of potential acceptor index positons and corresponding MaxEntScan scores.
}
