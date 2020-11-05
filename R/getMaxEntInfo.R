## MaxEntScan score calculations per SA site
getMaxEntInfo <- function(seq){

  ## input seq must consist of a c g or t 
  if(!all(strsplit(seq,"")[[1]] %in% c("a","c","g","t","G","C","T","A"))){
    stop('The entered sequence must be a character string of A C G and T.')
  }

  ## make sequence upper case
  um <- toupper(seq)

  ## Create index for MaxEntScan score calcualtions
  durchzahl <- seq_len(nchar(um)-3)
  d <- durchzahl
  durchzahl <- data.frame(durchzahl)

  ## Write sequence for MaxEntScan score calcualtions
  for(i in d){
    durchzahl$seq9[i] <- um
    durchzahl$seq9[i] <- substr(durchzahl$seq9[i],i,(i+22))
  }

  ## Discard sequences which are too short
  durchzahl$t <- nchar(durchzahl$seq9)
  durchzahl <- durchzahl[durchzahl$t == 23,]

  durchzahl$AG_present <- substr(durchzahl$seq9, 19,20)
  durchzahl$MaxEntScanScore <- as.numeric(calculateMaxEntScanScore(durchzahl$seq9, ssType=3))

  ## Return data frame
  return(durchzahl)
}
