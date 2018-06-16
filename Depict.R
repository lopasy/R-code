# Depict predicted
#{par(mar = c(10,5,3,2))
# Clean data
#tissueenrichment2 = tissueenrichment[order(tissueenrichment$MeSH.second.level.term),]
#stopwords = c("System")
#tissueenrichment2$MeSH.second.level.term = gsub(paste0(stopwords,collapse = "|"), "", tissueenrichment2$MeSH.second.level.term)
#tissueenrichment2$MeSH.second.level.term = as.factor(tissueenrichment2$MeSH.second.level.term)


#plot(tissueenrichment2$MeSH.second.level.term, -log10(tissueenrichment2$Nominal.P.value),
#    main = "Predicted avMSE GWAS", xaxt = "n", xlab = "", ylab = expression('-log'[10] (paste(italic(P)," value"))))
#tlab = 1:length(table(tissueenrichment2$MeSH.second.level.term))
#axis(1, at = tlab, labels = FALSE, tick = F)
#text(x = tlab, y = par()$usr[3]*(par()$usr[4]-par()$usr[3]),
#     labels = unique(tissueenrichment2$MeSH.second.level.term), srt = 45, adj = 1, xpd = TRUE)
#title(xlab = "MeSH  second - level  term", line = 7)
#}


#par(mar = c(8,5,3,2))
tissueenrichment2 = tissueenrichment
stopwords = c("System", "Systems")
tissueenrichment2$MeSH.second.level.term = gsub(paste0(stopwords,collapse = "|"), "", tissueenrichment2$MeSH.second.level.term)
tissueenrichment2$MeSH.second.level.term = as.factor(tissueenrichment2$MeSH.second.level.term)
tissueenrichment2 = tissueenrichment2[order(tissueenrichment2$MeSH.second.level.term),]

function.add.spacers <- function(df) {
  ### USE: UNIVERSIAL
  ### Function to INSERT spacers (dummy rows) in data frame between "groups"
  ### *IMPORTANT*:we are relying the data frame being (correctly) SORTED!!
  df.tissue_enrichment.spaced <- data.frame()
  group.previously <- df$MeSH.second.level.term[1]
  for (i in 1:nrow(df)) {
    ### *OBS*: df MUST be sorted at this point!
    group.now <- df$MeSH.second.level.term[i]
    if ( group.now != group.previously ) { # add "spacer" (blank) rows between "groups"
      df.dummy <- data.frame(Name="dummy",Nominal.P.value=1,False.discovery.rate=">=0.20")
      for (j in 1:n.blank_spacers) {
        # use dplyr::bind_rows() to fill unmatched columns with NA [rbind() will complain]
        df.tissue_enrichment.spaced <- suppressWarnings(dplyr::bind_rows(df.tissue_enrichment.spaced, df.dummy))
      }
      group.previously <- df$MeSH.second.level.term[i]
    }
    df.tissue_enrichment.spaced <- dplyr::bind_rows(df.tissue_enrichment.spaced, df[i,])
    # GTEx --> Name  Nominal.P.value  False.discovery.rate
  }
  return(df.tissue_enrichment.spaced)
}

scaleFUN <- function(x) sprintf("%.1f", x) # change number of decimal places

function.plot.barplot <- function(df.tissue_enrichment, xlabel="MISSING X-LABEL") {
  ### Constructing labels and break positions for x-axis
  # *OBS* df.tissue_enrichment is used for input and NOT the "spaced" data frame
  #df.tissue_enrichment <- as.data.frame(sapply(df.tissue_enrichment, gsub, pattern = "Hemic and Immune s", replacement = "Hemic and Immune", fixed = TRUE))
  df.breaks_and_labels <- df.tissue_enrichment %>% group_by(MeSH.second.level.term) %>% summarize(group.width=n()) # *OBS*: this code relies on the "group" being in the correct order
  df.breaks_and_labels$order.group.numeric <- seq(0,nrow(df.breaks_and_labels)-1) # 0,1,...,3
  df.breaks_and_labels$group.width.cumsum <- with(df.breaks_and_labels, cumsum(group.width)) # cumsum()
  df.breaks_and_labels$group.width.cumsum.shift <- with(df.breaks_and_labels, c(0,group.width.cumsum[-length(group.width.cumsum)])) # *OBS*: "shifting" position by one (pop array). Inserting zero in first position | try diff(x)?
  df.breaks_and_labels$break_position <- with(df.breaks_and_labels, 0.5+(n.blank_spacers*order.group.numeric)+group.width.cumsum.shift+group.width/2) 
  
  ### Adding spacers to data frame | *OBS*: we are relying on CORRECT SORTING of the data frame
  df.tissue_enrichment.spaced <- function.add.spacers(df.tissue_enrichment)
  df.tissue_enrichment.spaced$Order <- 1:nrow(df.tissue_enrichment.spaced) # Add numeric code | maps to x-axis
  df.tissue_enrichment.spaced$FDR.significant <- with(df.tissue_enrichment.spaced, ifelse( (False.discovery.rate=="<0.01"|False.discovery.rate=="<0.05"), TRUE, FALSE))

  p <- ggplot(df.tissue_enrichment.spaced)
  p <- p + geom_bar(aes(x=Order, y=-log10(Nominal.P.value)), stat="identity")
  p <- p + scale_fill_manual(name="", values=c("TRUE"="#F15A22","FALSE"="#606060",guide='legend')) # set colors for FDR significant bars
  p <- p + scale_x_continuous(breaks=df.breaks_and_labels$break_position, labels=df.breaks_and_labels$MeSH.second.level.term) # add x-axis labels
  p <- p + labs(y=expression(-log[10](paste(italic(P)," value")))) # add axis titles
  p <- p + guides(fill=FALSE) # remove legend
  p.theme <- theme_classic() # default is "theme_gray()" "theme_classic()" theme removes background and more [white background, eliminates background, gridlines, and chart border]. Note that it contains black colored axis lines. Therefore it is used as a template for further modificationsp.theme$axis.line.x <- element_blank() # empty list of class "element_blank" and "element" | SAME as theme_bw()$axis.line
  p.theme$axis.line.y <- p.theme$axis.line # copy list
  p.theme$axis.line.y$colour <- 'black' # OBS: sensitive to spelling: "colour" and NOT "color"
  p.theme$axis.ticks.x <- element_blank()
  p.theme$axis.ticks.y <- p.theme$axis.ticks
  p.theme$axis.ticks.y$colour <- 'black' # redundant: theme_bw() has this by default
  p.theme$axis.ticks.length <- unit(0.15, "cm") # adjust distance from axis labels to axis | you may need to play with this a little
  p <- p + scale_y_continuous(expand = c(0, 0),breaks=seq(0,2,0.2), limits = c(0, 2), labels = scaleFUN) # removes expansion of y-axis. Now the y-axis starts at zero!
  p.theme <- p.theme + theme(axis.text.x=element_text(angle=35, hjust=1, size=rel(1.2))) # size=rel(1.15) consider using %+replace% operator
  p <- p + p.theme # saving "p.theme" into plot
  p <- p + labs(x=xlabel)
  return(p) # return ggplot object
}


function.plot.save <- function(p, filename.prefix="MISSING-FILENAME-PREFIX", filename.suffix="MISSING-FILENAME-SUFFIX") {
  ### USE: UNIVERSIAL
  ### Function exports a ggplot object to pdf file
  ### INPUT: 
  # p: a ggplot object
  ### OUTPUT
  # NONE
  # ======================== Save plot =============================== # 
  #p.filename <- sprintf("tissue_plot_%s.pdf", filename.suffix)
  p.filename <- sprintf("tissue_plot_%s_%s.pdf", filename.prefix, filename.suffix)
  suppressWarnings(ggsave(file=p.filename, width=8, height=4, units="in", dpi=300))
  print(sprintf("Saved plot %s", p.filename))
}
