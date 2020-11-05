#Generate HEXplorer plot
generateHEXplorerPlot <- function(variationInfoList, ntWindow=20){

  ## Get HZEI scores per nt of reference sequence
  ## with and without variation
  durchzahl <- calculateHZEIperNT(as.character(variationInfoList$sequence))
  durchzahl2 <- calculateHZEIperNT(as.character(variationInfoList$altSeq))
  durchzahl$Sequence <- "reference"
  durchzahl2$Sequence <- "reference with variation"
  durchzahl <- rbind(durchzahl, durchzahl2)

  ## Get HBond scores per embeded donor sequences and their position
  durchzahl$hbs <- hbg$hbs[match(durchzahl$seq9,hbg$seq)]
  durchzahl$hbs[is.na(durchzahl$hbs)] <- 0
  durchzahl$hbs_pos <- durchzahl$durchzahl-2

  ## Prepare dataframe for plotting
  durchzahl$id <- seq_len(nrow(durchzahl))
  durchzahl$size <- 3
  durchzahl$size[durchzahl$durchzahl == ntWindow-4] <- 5
  durchzahl$index_nt <- substr(durchzahl$seq9,6,6 )
  durchzahl$ff <- 1
  durchzahl$ff[durchzahl$durchzahl == ntWindow-4] <- 2
  durchzahl$hbs[durchzahl$hbs_pos < 1] <- 0
  
  endhex <- 1
  Sequence <- "A"
  
  ## Generate HEXplorer plot
  plotOutput <- ggplot(durchzahl, aes(x = durchzahl, y = endhex ,fill=Sequence)) +
    geom_bar(stat='identity', position = "dodge")+ xlab("Sequence")+
    scale_y_continuous(name="Hexplorer score",
                       breaks=c(seq(-75,0,5),seq(2,34,2)),
                       limits=c(min(durchzahl$endhex)-6,
                                max(c(durchzahl$hbs,durchzahl$endhex))+1) )+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    annotate("text", label = durchzahl$index_nt[durchzahl$Sequence=="reference"] ,
             x= seq_len((nrow(durchzahl))/2),
             y = min(durchzahl$endhex-2), size = durchzahl$size[durchzahl$Sequence=="reference"],
             colour = "skyblue1", fontface = durchzahl$ff[durchzahl$Sequence=="reference"])+
    annotate("text", label = durchzahl$index_nt[durchzahl$Sequence=="reference with variation"] ,
             x= seq_len((nrow(durchzahl))/2), y = min(durchzahl$endhex-4),
             size = durchzahl$size[durchzahl$Sequence=="reference with variation"], colour = "gray30",
             fontface = durchzahl$ff[durchzahl$Sequence=="reference with variation"])+
    scale_fill_manual(values=c("#56B4E9", "#000000"))

  ## Add HBS to plot for reference sequence
  for(i in durchzahl$id[(durchzahl$hbs>0)&(durchzahl$Sequence=="reference")]){

    f <- durchzahl[durchzahl$Sequence=="reference",]
    x <- durchzahl$hbs_pos[i]-0.2
    yend <- durchzahl$hbs[i]
    plotOutput <- plotOutput+
      geom_segment(x = x, y = 0, xend = x, yend = yend,
                   size=1.5,aes(color="HBS reference"))

  }

  ## Add HBS to plot for reference sequence with variation
  for(i in durchzahl$id[(durchzahl$hbs>0)&(durchzahl$Sequence=="reference with variation")]){

    f <- durchzahl[durchzahl$Sequence=="reference with variation",]
    x <- durchzahl$hbs_pos[i]+0.2
    yend <- durchzahl$hbs[i]
    plotOutput <- plotOutput+
      geom_segment(x = x, y = 0, xend = x, yend = yend,
                   size=1.5,aes(color="HBS reference with variation"))

  }

  ## Add HBS labels
  hbsM <- durchzahl[(durchzahl$hbs>0)&(durchzahl$Sequence == "reference with variation"),]
  hbsW <- durchzahl[(durchzahl$hbs>0)&(durchzahl$Sequence == "reference"),]
  plotOutput <- plotOutput + annotate("text", x = hbsW$hbs_pos-.6 ,
                                      y= hbsW$hbs+1 , label = as.character(hbsW$hbs))
  plotOutput <- plotOutput + annotate("text", x = hbsM$hbs_pos+.6 ,
                                      y= hbsM$hbs+1 , label = as.character(hbsM$hbs))

  ## Add MaxEntScan score to plot
  durchzahlMax <- getMaxEntInfo(as.character(variationInfoList$sequence))
  durchzahlMax2 <- getMaxEntInfo(as.character(variationInfoList$altSeq))
  durchzahlMax$Sequence <- "reference"
  durchzahlMax2$Sequence <- "reference with variation"
  durchzahlMax <- rbind(durchzahlMax,durchzahlMax2)
  durchzahl$seq_ID <- paste(durchzahl$seq9, durchzahl$durchzahl)
  durchzahlMax$seq_ID <- paste(substr(durchzahlMax$seq9,1,11), durchzahlMax$durchzahl)

  durchzahl$MaxEntScanScore <- durchzahlMax$MaxEntScanScore[match(durchzahl$seq_ID, 
                                                                  durchzahlMax$seq_ID)]
  durchzahl[is.na(durchzahl)] <- 0
  durchzahl$MaxEnt_pos <- durchzahl$hbs_pos+15

  for(i in durchzahl$durchzahl[(durchzahl$MaxEntScanScore>0)&(durchzahl$Sequence=="reference")]){

    f <- durchzahl[durchzahl$Sequence=="reference",]
    x <- durchzahl$MaxEnt_pos[i]-0.2
    yend <- durchzahl$MaxEntScanScore[i]
    plotOutput <- plotOutput+
      geom_segment(x = x, y = 0, xend = x, yend = yend, 
                   size=1.5,aes(color="MaxEnt score reference"))

  }

  for(i in durchzahl$id[(durchzahl$MaxEntScanScore>0)&
                        (durchzahl$Sequence=="reference with variation")]){

    f <- durchzahl[durchzahl$Sequence=="reference with variation",]
    x <- durchzahl$MaxEnt_pos[i]+0.2
    yend <- durchzahl$MaxEntScanScore[i]
    plotOutput <- plotOutput+
      geom_segment(x = x, y = 0, xend = x, yend = yend,
                   size=1.5,aes(color="MaxEnt score reference with variation"))

  }

  ## Add MaxEnt score lables to plot
  mesW <- durchzahl[(durchzahl$MaxEntScanScore>0)&
                       (durchzahl$Sequence == "reference"),]
  mesM <- durchzahl[(durchzahl$MaxEntScanScore>0)&
                       (durchzahl$Sequence == "reference with variation"),]

  plotOutput <- plotOutput + annotate("text", x = mesW$MaxEnt_pos-0.6 ,
                                      y= mesW$MaxEntScanScore+1 , 
                                      label = as.character(round(mesW$MaxEntScanScore,1)))
  plotOutput <- plotOutput + annotate("text", x = mesM$MaxEnt_pos+0.6 ,
                                      y= mesM$MaxEntScanScore+1 ,
                                      label = as.character(round(mesM$MaxEntScanScore,1)))

  plotOutput <- plotOutput + 
    scale_colour_manual(name='Splite site scores',
                       values=c('HBS reference'='#F0E442',
                              'HBS reference with variation'='#E69F00',
                               'MaxEnt score reference'='#FF6A6A', 
                              'MaxEnt score reference with variation'='#8B3A3A'))

  ## Return plot
  return(plotOutput)
}
