## Prepare the uploaded reference DNA string set
prepareReferenceFasta <- function(filepath){

  referenceDnaStringSet2 <- readDNAStringSet(filepath, format="fasta",use.names=TRUE)
  ref_names <- as.character(lapply(names(referenceDnaStringSet2),
                                   function(x){ strsplit(x, " ")[[1]][[1]]}))
  names(referenceDnaStringSet2) <- ref_names

  return(referenceDnaStringSet2)


}


