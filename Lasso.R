# Lasso

library(data.table); library(qqman); library(ggplot2)
chr = list.files(pattern = ".lasso", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
colnames(d)[2] = "rs"
lasso = merge(stats, d, by = "rs")



cor.test(lasso$effalt,lasso$EFFECT)
plot(lasso$effalt,lasso$EFFECT)


chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
        "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
for (i in chr){
  assign(paste0(i), lasso[which(lasso$chr == i),])
}

{chr1$seq=1:nrow(chr1)
  chr2$seq=1:nrow(chr2); chr2$seq = chr2$seq + (range(chr1$seq)[2]+100)
  chr3$seq=1:nrow(chr3); chr3$seq = chr3$seq + (range(chr2$seq)[2]+100)
  chr4$seq=1:nrow(chr4); chr4$seq = chr4$seq + (range(chr3$seq)[2]+100)
  chr5$seq=1:nrow(chr5); chr5$seq = chr5$seq + (range(chr4$seq)[2]+100)
  chr6$seq=1:nrow(chr6); chr6$seq = chr6$seq + (range(chr5$seq)[2]+100)
  chr7$seq=1:nrow(chr7); chr7$seq = chr7$seq + (range(chr6$seq)[2]+100)
  chr8$seq=1:nrow(chr8); chr8$seq = chr8$seq + (range(chr7$seq)[2]+100)
  chr9$seq=1:nrow(chr9); chr9$seq = chr9$seq + (range(chr8$seq)[2]+100)
  chr10$seq=1:nrow(chr10); chr10$seq = chr10$seq + (range(chr9$seq)[2]+100)
  chr11$seq=1:nrow(chr11); chr11$seq = chr11$seq + (range(chr10$seq)[2]+100)
  chr12$seq=1:nrow(chr12); chr12$seq = chr12$seq + (range(chr11$seq)[2]+100)
  chr13$seq=1:nrow(chr13); chr13$seq = chr13$seq + (range(chr12$seq)[2]+100)
  chr14$seq=1:nrow(chr14); chr14$seq = chr14$seq + (range(chr13$seq)[2]+100)
  chr15$seq=1:nrow(chr15); chr15$seq = chr15$seq + (range(chr14$seq)[2]+100)
  chr16$seq=1:nrow(chr16); chr16$seq = chr16$seq + (range(chr15$seq)[2]+100)
  chr17$seq=1:nrow(chr17); chr17$seq = chr17$seq + (range(chr16$seq)[2]+100)
  chr18$seq=1:nrow(chr18); chr18$seq = chr18$seq + (range(chr17$seq)[2]+100)
  chr19$seq=1:nrow(chr19); chr19$seq = chr19$seq + (range(chr18$seq)[2]+100)
  chr20$seq=1:nrow(chr20); chr20$seq = chr20$seq + (range(chr19$seq)[2]+100)
  chr21$seq=1:nrow(chr21); chr21$seq = chr21$seq + (range(chr20$seq)[2]+100)
  chr22$seq=1:nrow(chr22); chr22$seq = chr22$seq + (range(chr21$seq)[2]+100)}



chr = list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,
           chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
lasso1 = do.call(rbind,chr)


lasso = lasso[order(lasso$CHR),]
lasso$seq = 1:nrow(lasso)
{gg <- ggplot(lasso1) 
gg <- gg + geom_segment(aes(x=seq, y=EFFECT, xend=seq, yend=0, color=CHR), show.legend=FALSE)
gg <- gg + geom_point(aes(x=seq, y=EFFECT, color=CHR))
gg <- gg + theme_bw()
gg <- gg + theme(strip.background=element_blank())
gg <- gg + theme(strip.text=element_text(hjust=0))
gg <- gg + theme(panel.grid.major.x=element_blank())
gg <- gg + theme(panel.grid.minor.y=element_blank())
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text.x=element_text(size=8))
gg <- gg + theme(axis.text.y=element_text(size=8, vjust=c(0, 0.5, 0.5, 0.5, 1)))
gg <- gg + theme(legend.key=element_blank())
gg <- gg + theme(legend.position="right")
gg <- gg + theme(legend.title = element_text(size = 10,face="bold"))
gg <- gg + geom_text(x = 9100, y = 0.0585, parse = T, label = as.character(expression(paste(N["samples"], "=", 71973, sep = " "))), size = 5)
gg <- gg + geom_text(x = 9100, y = 0.055, label = "MAF > 0.01", size = 5)
gg <- gg + geom_text(x = 9100, y = 0.0515, label = "Total Hits = 8601", size = 5)
gg <- gg + geom_text(x = 9100, y = 0.048, parse = T, label = as.character(expression(paste(h^"2", "NULL%~~%", 0.38, sep = " "))), size = 5)
gg <- gg + geom_text(x = 9100, y = 0.0445, label = "?? = 0.01098, CI(0.010965, 0.010998)", size = 5)
gg <- gg + geom_text(x = 9100, y = 0.0410, parse = T, label = as.character(expression(paste(N["SNP"], "NULL%~~%", "8.5M", sep =" "))), size = 5)
gg <- gg + ggtitle("LASSO effect estimates")
gg <- gg + guides(colour = guide_legend(reverse=T))
gg <- gg + theme(plot.title = element_text(hjust = 0.5, size = 20, face="bold", margin = margin(t = 0, r = 20, b = 20, l = 0)), 
      legend.title = element_text(), 
      axis.text= element_text(size=12), axis.title=element_text(size=16,face="bold"),
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))
gg <- gg + scale_x_continuous(breaks=seq(0,11000,1000))
gg <- gg + scale_y_continuous(breaks=seq(-0.06,0.06,0.01))
gg <- gg + labs(x = "Allignment", y=expression(??["i"]))
gg}




hist(lasso$EFFECT)
hist(lasso$effalt)

hist(lasso$CHR)
table(lasso$CHR)


manhattan(lasso, snp = "rs", bp = "pos", p = "pval")



# Stratify by MAF
lasso$sign[lasso$EFFECT > 0] = "Positive"; lasso$sign[lasso$EFFECT < 0] = "Negative"
lasso$maf = 1- lasso$reffrq


lasso$maf_adjust[lasso$maf <= 0.04] = 0.04



{lab = c(0.01, 0.04, 0.07, 0.1, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34, 0.37, 0.4, 0.43, 0.46, 0.49)
h1 = hist(lasso$maf, plot = FALSE, breaks = 50)
h2 = hist(lasso$maf, plot = FALSE, breaks = 50)
h2$counts = - h2$counts
hmax = max(h1$counts)
hmin = min(h2$counts)
X = c(h1$breaks, h2$breaks)
xmax = max(X)
xmin = min(X)
plot(h1, ylim = c(hmin, hmax), col = "green", xaxt = "n", xlab = "MinorAllele Frequency",
     ylab = "Number of nonzero effects", main = "Distribution of effect sizes", cex.axis = 1.2, 
     cex.lab = 1.5, cex.main = 1.5)
lines(h2, col = "blue")
axis(1, at = lab,labels = lab, las = 1, cex.axis = 1.2)}








write.table(chr1, "chr1.txt", row.names = F, quote = F)
write.table(chr2, "chr2.txt", row.names = F, quote = F)
write.table(chr3, "chr3.txt", row.names = F, quote = F)
write.table(chr4, "chr4.txt", row.names = F, quote = F)
write.table(chr5, "chr5.txt", row.names = F, quote = F)
write.table(chr6, "chr6.txt", row.names = F, quote = F)
write.table(chr7, "chr7.txt", row.names = F, quote = F)
write.table(chr8, "chr8.txt", row.names = F, quote = F)
write.table(chr9, "chr9.txt", row.names = F, quote = F)
write.table(chr10, "chr10.txt", row.names = F, quote = F)
write.table(chr11, "chr11.txt", row.names = F, quote = F)
write.table(chr12, "chr12.txt", row.names = F, quote = F)
write.table(chr13, "chr13.txt", row.names = F, quote = F)
write.table(chr14, "chr14.txt", row.names = F, quote = F)
write.table(chr15, "chr15.txt", row.names = F, quote = F)
write.table(chr16, "chr16.txt", row.names = F, quote = F)
write.table(chr17, "chr17.txt", row.names = F, quote = F)
write.table(chr18, "chr18.txt", row.names = F, quote = F)
write.table(chr19, "chr19.txt", row.names = F, quote = F)
write.table(chr20, "chr20.txt", row.names = F, quote = F)
write.table(chr21, "chr21.txt", row.names = F, quote = F)
write.table(chr22, "chr22.txt", row.names = F, quote = F)











