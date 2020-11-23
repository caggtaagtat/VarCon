## Generate the text for describing the difference between the Hexplorer scores of both sequences
getSeqInfoFromVariation <- function(referenceDnaStringSet,
                                        transcriptID, variation,
                                        ntWindow=20, transcriptTable,
                                        gene2transcript=gene2transcript){
  
  if(!transcriptID %in% transcriptTable$Transcript.stable.ID){

      ## Use gene-to-transcript conversion table if needed
      if(transcriptID %in% gene2transcript$gene_name){
        idx <- match(transcriptID, gene2transcript$gene_name)
         transcriptID <- gene2transcript$transcriptID[idx]
      }
  
  
      if(transcriptID %in% gene2transcript$geneID){
        idx <- match(transcriptID, gene2transcript$geneID)
        transcriptID <- gene2transcript$transcriptID[idx]
       }
  }

  ## Get Info from variation annotation
  res <- list()
  trans <- as.character(transcriptID)
  res["transcript"] <- trans
  anno <- strsplit(as.character(variation),"/")[[1]][[1]]
  res["funcAnnotation"] <- as.character(variation)
  findArrow <- regexpr(">", anno)[[1]]

  ## Get reference and variation nucleotide at SNV position
  fromCoding <- substr(anno,3, findArrow-2 )
  res["ref_nuc"] <- substr(anno,findArrow-1,findArrow-1 )
  res["alt_nuc"] <- substr(anno,findArrow+1,findArrow+1 )

  ## Calculate genomic coordinate from SNV variation annotation
  ## either from SNV annotaton refferring to the coding or genomic sequence
  if(substr(anno,1, 1 )=="c"){

    fromCodingMinus <- 0
    fromCodingPlus <- 0

    if(regexpr("-", fromCoding)[[1]] != -1){

      fromCodingMinus <- strsplit(fromCoding,"-")[[1]][[2]]
      fromCoding <- strsplit(fromCoding,"-")[[1]][[1]]


    }

    if(regexpr("\\+", fromCoding)[[1]] != -1){

      fromCodingPlus <- strsplit(fromCoding,"\\+")[[1]][[2]]
      fromCoding <- strsplit(fromCoding,"\\+")[[1]][[1]]


    }
    
    logicVector <- transcriptTable$Transcript.stable.ID == trans
    transcriptSubset <- transcriptTable[logicVector, ]
    transcriptSubset$found <- 0
    transcriptSubset[is.na(transcriptSubset)] <- 0

    for(inns in seq_len(nrow(transcriptSubset))){

      if(as.numeric(fromCoding) %in% transcriptSubset$CDS.start[inns]:
         transcriptSubset$CDS.end[inns]) transcriptSubset$found[inns] <- 1

    }

    strand <- transcriptSubset$Strand[1]

    if(strand == 1){
      logicVector <- transcriptSubset$found == 1
      genomicCoord <- as.numeric(fromCoding)-
        transcriptSubset$CDS.start[logicVector]
      genomicCoord <- transcriptSubset$Exon.region.start..bp.[logicVector]+
        genomicCoord
      genomicCoord <- genomicCoord - as.numeric(fromCodingMinus)
      genomicCoord <- genomicCoord + as.numeric(fromCodingPlus)

       if(transcriptSubset$CDS.start[logicVector] == 1){

        coLeng <- transcriptSubset$cDNA.coding.end[logicVector]-
          transcriptSubset$cDNA.coding.start[logicVector]
        exCodingStart <- transcriptSubset$Exon.region.end..bp.[logicVector]-
          transcriptSubset$Exon.region.start..bp.[logicVector]
        exCodingStart <- exCodingStart - coLeng

        genomicCoord <- as.numeric(fromCoding)-
          transcriptSubset$CDS.start[logicVector]
        genomicCoord <- transcriptSubset$Exon.region.start..bp.[logicVector]+
          exCodingStart+genomicCoord
        genomicCoord <- genomicCoord - as.numeric(fromCodingMinus)
        genomicCoord <- genomicCoord + as.numeric(fromCodingPlus)

      }

    }else{
      logicVector <- transcriptSubset$found == 1
      genomicCoord <- as.numeric(fromCoding)-
        transcriptSubset$CDS.start[logicVector]
      genomicCoord <- transcriptSubset$Exon.region.end..bp.[logicVector]-
        genomicCoord
      genomicCoord <- genomicCoord + as.numeric(fromCodingMinus)
      genomicCoord <- genomicCoord - as.numeric(fromCodingPlus)

      if(transcriptSubset$CDS.start[logicVector] == 1){

        coLeng <- transcriptSubset$cDNA.coding.end[logicVector]-
          transcriptSubset$cDNA.coding.start[logicVector]
        exCodingStart <- transcriptSubset$Exon.region.end..bp.[logicVector]-
          transcriptSubset$Exon.region.start..bp.[logicVector]
        exCodingStart <- exCodingStart - coLeng

        genomicCoord <- as.numeric(fromCoding)-
          transcriptSubset$CDS.start[logicVector]
        genomicCoord <- transcriptSubset$Exon.region.end..bp.[logicVector]-
          exCodingStart-genomicCoord
        genomicCoord <- genomicCoord + as.numeric(fromCodingMinus)
        genomicCoord <- genomicCoord - as.numeric(fromCodingPlus)

      }

    }}else if(substr(anno,1, 1 )=="g"){

      genomicCoord <- fromCoding

      if(strand == -1){

        RefAltNuc <- c( res["ref_nuc"], res["alt_nuc"])

        RefAltNuc <- gsub("A","t",RefAltNuc)
        RefAltNuc <- gsub("C","g",RefAltNuc)
        RefAltNuc <- gsub("T","a",RefAltNuc)
        RefAltNuc <- gsub("G","c",RefAltNuc)
        RefAltNuc <- toupper(RefAltNuc)

        #Change nucleotides if gene on minus strand
        res["ref_nuc"] <- RefAltNuc[1]
        res["alt_nuc"] <- RefAltNuc[2]

      }



    }else{
      stop("Entered variation not in reference to Genome or CDS (g. or c.).")
    }

  res["genomicCoordinate"] <- genomicCoord
  
  ## Retrieve strand and chromosome of the transcript 
  logicVector <- transcriptTable$Transcript.stable.ID == trans
  transcriptSubset <- transcriptTable[logicVector, ]
  strand <-     unique(transcriptSubset$Strand)
  chromosome <- unique(transcriptSubset$Chromosome.scaffold.name )

  ## Get the respective genomic sequence
  sub.info <- data.frame(chrom=as.character(chromosome), 
                         start=genomicCoord-ntWindow,
                         end=genomicCoord+ntWindow)
  seqRange <- as.character(getSeq(referenceDnaStringSet, as(sub.info, "GRanges")))

  ## Minus strand sequence correction
  seqRange[strand=="-1"] <- gsub("A","t",seqRange[strand=="-1"])
  seqRange[strand=="-1"] <- gsub("C","g",seqRange[strand=="-1"])
  seqRange[strand=="-1"] <- gsub("T","a",seqRange[strand=="-1"])
  seqRange[strand=="-1"] <- gsub("G","c",seqRange[strand=="-1"])
  seqRange <- toupper(seqRange)
  seqRange[strand=="-1"] <- reverse(seqRange[strand=="-1"])

  ## Save the sequence surrounding
  res["sequence"] <-   seqRange

  ## Substitut nucleotide at the middle
  splSeq <- strsplit(res$sequence, "")[[1]]
  splSeq[ median(seq_len(length(splSeq))) ] <- res$alt_nuc
  altnative_seq <- paste(splSeq, collapse = "")
  res["altSeq"]  <- altnative_seq

  res["genomic_range"] <- paste(c((genomicCoord-ntWindow),
                                  (genomicCoord+ntWindow)),collapse = " ")

  ## Return gathered information
  return(res)
}
