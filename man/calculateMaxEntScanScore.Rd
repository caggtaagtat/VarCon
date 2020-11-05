\name{calculateMaxEntScanScore}
\alias{calculateMaxEntScanScore}
\title{Calculate MaxEntScan score of a splice site sequence}
\usage{
calculateMaxEntScanScore(seqVector, ssType)
}
\arguments{
\item{seqVector}{Character vector of nucleotide sequence of a splice site sequences. SA sequences should be 23nt long (20 intronic, 3 exonic) and splice donor sequences should be 9nt long (3 exonic, 6 intronic) only contain bases "A", "G", "C", "T".}
\item{ssType}{Numeric indicator, if the entred sequence is a splice donor (5) or acceptor (3)}
}
\description{
This function calculates the MaxEntScan score of either splice donor or acceptor sequences.
}
\value{
Character vector of the MaxEntScan scores generated from the entered \code{seqVector}.
}
\examples{
calculateMaxEntScanScore("TTCCAAACGAACTTTTGTAGGGA",3)
calculateMaxEntScanScore("GAGGTAAGT",5)

}
