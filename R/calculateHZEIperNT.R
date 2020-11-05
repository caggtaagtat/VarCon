## HEXplorer plotting Function
calculateHZEIperNT <- function(seq){

  ## if there is something else than a c g or t in the sequence, give an error
  if(!all(strsplit(seq,"")[[1]] %in% c("a","c","g","t","G","C","T","A"))){
    stop('The entered sequence must be a character string of A C G and T.')
  } 

  ## Doing the hexplorer plotting 
  um <- toupper(seq)

  ## Create Index for HZEI calculations
  durchzahl <- seq_len(nchar(um)-8)

  ## Make dataframe from index
  d <- durchzahl
  durchzahl <- data.frame(durchzahl)

  ## Get 11nt long sequences overlapping by 10 nt
  for(i in d){

    durchzahl$seq9[i] <- um
    durchzahl$seq9[i] <- substr(durchzahl$seq9[i],i,(i+10))
  }

  ## Discard the sequences which are too short
  durchzahl$t <- nchar(durchzahl$seq9)
  durchzahl <- durchzahl[durchzahl$t == 11,]


  ## Create Sequences for HEXplorer Score entry
  durchzahl$hex1 <- lapply(durchzahl$seq9, function(x) substr(x,1,6))
  durchzahl$hex2 <- lapply(durchzahl$seq9, function(x) substr(x,2,7))
  durchzahl$hex3 <- lapply(durchzahl$seq9, function(x) substr(x,3,8))
  durchzahl$hex4 <- lapply(durchzahl$seq9, function(x) substr(x,4,9))
  durchzahl$hex5 <- lapply(durchzahl$seq9, function(x) substr(x,5,10))
  durchzahl$hex6 <- lapply(durchzahl$seq9, function(x) substr(x,6,11))


  ## Insert HEXplorer Z-Scores
  mtc <- match(durchzahl$hex1, hex$seq)
  durchzahl$hex1value <- hex$value[mtc]

  mtc <- match(durchzahl$hex2, hex$seq)
  durchzahl$hex2value <- hex$value[mtc]

  mtc <- match(durchzahl$hex3, hex$seq)
  durchzahl$hex3value <- hex$value[mtc]

  mtc <- match(durchzahl$hex4, hex$seq)
  durchzahl$hex4value <- hex$value[mtc]

  mtc <- match(durchzahl$hex5, hex$seq)
  durchzahl$hex5value <- hex$value[mtc]

  mtc <- match(durchzahl$hex6, hex$seq)
  durchzahl$hex6value <- hex$value[mtc]


  ## Calculate HEXplorer Score per nucleotide, by the mean of all 
  ## flanking Hexamer Z-Scores, and round the result
  durchzahl$endhex <- round((durchzahl$hex1value+durchzahl$hex2value+
                               durchzahl$hex3value+durchzahl$hex4value+
                               durchzahl$hex5value+durchzahl$hex6value)/6, digits=2)

  ## Return right columns
  d <- durchzahl[,-c(3:15)]
  return(d)
}
