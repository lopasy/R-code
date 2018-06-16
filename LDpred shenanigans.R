library(data.table); library(qqman)
chr = list.files(pattern = ".txt", full.names = TRUE)
chr = lapply(chr, fread)
chr = do.call(rbind.data.frame, chr)
chr2 = chr[which(chr$V10 < 0.005),]


manhattan(chr2,chr = "V2", bp = "V4",p = "V10",snp = "V3", cex = 0.75)


one=chr[which(chr$V2 == 1),]


a = a[which(a$info > 0.8),]

chr = 1:22
for (i in chr){
  assign(paste0("chr", i), a[which(a$V2 == i),])
}


# There are duplicates. Remove
{chr1 = chr1[!duplicated(chr1$pos),]
  chr2 = chr2[!duplicated(chr2$pos), ]
  chr3 = chr3[!duplicated(chr3$pos), ]
  chr4 = chr4[!duplicated(chr4$pos), ]
  chr5 = chr5[!duplicated(chr5$pos), ]
  chr6 = chr6[!duplicated(chr6$pos), ]
  chr7 = chr7[!duplicated(chr7$pos), ]
  chr8 = chr8[!duplicated(chr8$pos), ]
  chr9 = chr9[!duplicated(chr9$pos), ]
  chr10 = chr10[!duplicated(chr10$pos), ]
  chr11 = chr11[!duplicated(chr11$pos), ]
  chr12 = chr12[!duplicated(chr12$pos), ]
  chr13 = chr13[!duplicated(chr13$pos), ]
  chr14 = chr14[!duplicated(chr14$pos), ]
  chr15 = chr15[!duplicated(chr15$pos), ]
  chr16 = chr16[!duplicated(chr16$pos), ]
  chr17 = chr17[!duplicated(chr17$pos), ]
  chr18 = chr18[!duplicated(chr18$pos), ]
  chr19 = chr19[!duplicated(chr19$pos), ]
  chr20 = chr20[!duplicated(chr20$pos), ]
  chr21 = chr21[!duplicated(chr21$pos), ]
  chr22 = chr22[!duplicated(chr22$pos), ]}

#chr = list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,
#           chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
#stats1 = do.call(rbind,chr)

stats1 = stats1[which(stats1$pval < 0.005),]

manhattan(chr2,chr = "V2", bp = "V4",p = "V10",snp = "V3", cex = 0.75)

# Remove silly columns
{chr1 = chr1[,c(1:7,17,9)]
  chr2 = chr2[,c(1:7,17,9)]
  chr3 = chr3[,c(1:7,17,9)]
  chr4 = chr4[,c(1:7,17,9)]
  chr5 = chr5[,c(1:7,17,9)]
  chr6 = chr6[,c(1:7,17,9)]
  chr7 = chr7[,c(1:7,17,9)]
  chr8 = chr8[,c(1:7,17,9)]
  chr9 = chr9[,c(1:7,17,9)]
  chr10 = chr10[,c(1:7,17,9)]
  chr11 = chr11[,c(1:7,17,9)]
  chr12 = chr12[,c(1:7,17,9)]
  chr13 = chr13[,c(1:7,17,9)]
  chr14 = chr14[,c(1:7,17,9)]
  chr15 = chr15[,c(1:7,17,9)]
  chr16 = chr16[,c(1:7,17,9)]
  chr17 = chr17[,c(1:7,17,9)]
  chr18 = chr18[,c(1:7,17,9)]
  chr19 = chr19[,c(1:7,17,9)]
  chr20 = chr20[,c(1:7,17,9)]
  chr21 = chr21[,c(1:7,17,9)]
  chr22 = chr22[,c(1:7,17,9)]}

#Sort by pos
{chr1 = chr1[order(chr1$pos),]
  chr2 = chr2[order(chr2$pos),]
  chr3 = chr3[order(chr3$pos),]
  chr4 = chr4[order(chr4$pos),]
  chr5 = chr5[order(chr5$pos),]
  chr6 = chr6[order(chr6$pos),]
  chr7 = chr7[order(chr7$pos),]
  chr8 = chr8[order(chr8$pos),]
  chr9 = chr9[order(chr9$pos),]
  chr10 = chr10[order(chr10$pos),]
  chr11 = chr11[order(chr11$pos),]
  chr12 = chr12[order(chr12$pos),]
  chr13 = chr13[order(chr13$pos),]
  chr14 = chr14[order(chr14$pos),]
  chr15 = chr15[order(chr15$pos),]
  chr16 = chr16[order(chr16$pos),]
  chr17 = chr17[order(chr17$pos),]
  chr18 = chr18[order(chr18$pos),]
  chr19 = chr19[order(chr19$pos),]
  chr20 = chr20[order(chr20$pos),]
  chr21 = chr21[order(chr21$pos),]
  chr22 = chr22[order(chr22$pos),]}


{colnames(chr1)[8] = "pval"
  colnames(chr2)[8] = "pval"
  colnames(chr3)[8] = "pval"
  colnames(chr4)[8] = "pval"
  colnames(chr5)[8] = "pval"
  colnames(chr6)[8] = "pval"
  colnames(chr7)[8] = "pval"
  colnames(chr8)[8] = "pval"
  colnames(chr9)[8] = "pval"
  colnames(chr10)[8] = "pval"
  colnames(chr11)[8] = "pval"
  colnames(chr12)[8] = "pval"
  colnames(chr13)[8] = "pval"
  colnames(chr14)[8] = "pval"
  colnames(chr15)[8] = "pval"
  colnames(chr16)[8] = "pval"
  colnames(chr17)[8] = "pval"
  colnames(chr18)[8] = "pval"
  colnames(chr19)[8] = "pval"
  colnames(chr20)[8] = "pval"
  colnames(chr21)[8] = "pval"
  colnames(chr22)[8] = "pval"}




# Create overlaps
{test1 = chr1[which(chr1$pos < 62497273),]
test2 = chr1[which(chr1$pos > 61497273 & chr1$pos < 124994546),]
test3 = chr1[which(chr1$pos > 114994546 & chr1$pos < 200000000),]
test4 = chr1[which(chr1$pos > 199000000 & chr1$pos < 249237748),]

write.table(test1, "chr1_1.txt", row.names = F, quote = F)
write.table(test2, "chr1_2.txt", row.names = F, quote = F)
write.table(test3, "chr1_3.txt", row.names = F, quote = F)
write.table(test4, "chr1_4.txt", row.names = F, quote = F)
################################################################
test1 = chr2[which(chr2$pos < 50885714),]
test2 = chr2[which(chr2$pos > 49885714 & chr2$pos < 121551428),]
test3 = chr2[which(chr2$pos > 120551428 & chr2$pos < 182327142),]
test4 = chr2[which(chr2$pos > 181327142 & chr2$pos < 243092285),]

write.table(test1, "chr2_1.txt", row.names = F, quote = F)
write.table(test2, "chr2_2.txt", row.names = F, quote = F)
write.table(test3, "chr2_3.txt", row.names = F, quote = F)
write.table(test4, "chr2_4.txt", row.names = F, quote = F)
################################################################
test1 = chr3[which(chr3$pos < 49483401),]
test2 = chr3[which(chr3$pos > 48483401 & chr3$pos < 98966802),]
test3 = chr3[which(chr3$pos > 97966802 & chr3$pos < 148450203),]
test4 = chr3[which(chr3$pos > 147450203 & chr3$pos < 197873408),]

write.table(test1, "chr3_1.txt", row.names = F, quote = F)
write.table(test2, "chr3_2.txt", row.names = F, quote = F)
write.table(test3, "chr3_3.txt", row.names = F, quote = F)
write.table(test4, "chr3_4.txt", row.names = F, quote = F)
################################################################
test1 = chr4[which(chr4$pos < 47770824),]
test2 = chr4[which(chr4$pos > 47670824 & chr4$pos < 95541647),]
test3 = chr4[which(chr4$pos > 95441647 & chr4$pos < 142312471),]
test4 = chr4[which(chr4$pos > 142312471 & chr4$pos < 191014509),]

write.table(test1, "chr4_1.txt", row.names = F, quote = F)
write.table(test2, "chr4_2.txt", row.names = F, quote = F)
write.table(test3, "chr4_3.txt", row.names = F, quote = F)
write.table(test4, "chr4_4.txt", row.names = F, quote = F)
################################################################
test1 = chr5[which(chr5$pos < 59000000),]
test2 = chr5[which(chr5$pos > 58000000 & chr5$pos < 120603708),]
test3 = chr5[which(chr5$pos > 119603708 & chr5$pos < 180905563),]

write.table(test1, "chr5_1.txt", row.names = F, quote = F)
write.table(test2, "chr5_2.txt", row.names = F, quote = F)
write.table(test3, "chr5_3.txt", row.names = F, quote = F)
################################################################
test1 = chr6[which(chr6$pos < 52081557),]
test2 = chr6[which(chr6$pos > 51081557 & chr6$pos < 112163114),]
test3 = chr6[which(chr6$pos > 111163114 & chr6$pos < 171244672),]

write.table(test1, "chr6_1.txt", row.names = F, quote = F)
write.table(test2, "chr6_2.txt", row.names = F, quote = F)
write.table(test3, "chr6_3.txt", row.names = F, quote = F)
################################################################
test1 = chr7[which(chr7$pos < 50053064),]
test2 = chr7[which(chr7$pos > 49053064 & chr7$pos < 106106128),]
test3 = chr7[which(chr7$pos > 105106128 & chr7$pos < 159159193),]

write.table(test1, "chr7_1.txt", row.names = F, quote = F)
write.table(test2, "chr7_2.txt", row.names = F, quote = F)
write.table(test3, "chr7_3.txt", row.names = F, quote = F)
################################################################
test1 = chr8[which(chr8$pos < 46783364),]
test2 = chr8[which(chr8$pos > 45783364 & chr8$pos < 97566728),]
test3 = chr8[which(chr8$pos > 96566728 & chr8$pos < 146350093),]

write.table(test1, "chr8_1.txt", row.names = F, quote = F)
write.table(test2, "chr8_2.txt", row.names = F, quote = F)
write.table(test3, "chr8_3.txt", row.names = F, quote = F)
################################################################
test1 = chr9[which(chr9$pos < 47056615),]
test2 = chr9[which(chr9$pos > 46056615 & chr9$pos < 94113230),]
test3 = chr9[which(chr9$pos > 93113230 & chr9$pos < 141169846),]

write.table(test1, "chr9_1.txt", row.names = F, quote = F)
write.table(test2, "chr9_2.txt", row.names = F, quote = F)
write.table(test3, "chr9_3.txt", row.names = F, quote = F)
################################################################
test1 = chr10[which(chr10$pos < 45201572),]
test2 = chr10[which(chr10$pos > 44201572 & chr10$pos < 90403144),]
test3 = chr10[which(chr10$pos > 89403144 & chr10$pos < 135604717),]

write.table(test1, "chr10_1.txt", row.names = F, quote = F)
write.table(test2, "chr10_2.txt", row.names = F, quote = F)
write.table(test3, "chr10_3.txt", row.names = F, quote = F)
################################################################
test1 = chr11[which(chr11$pos < 45042384),]
test2 = chr11[which(chr11$pos > 44042384 & chr11$pos < 90084768),]
test3 = chr11[which(chr11$pos > 89084768 & chr11$pos < 135127152),]

write.table(test1, "chr11_1.txt", row.names = F, quote = F)
write.table(test2, "chr11_2.txt", row.names = F, quote = F)
write.table(test3, "chr11_3.txt", row.names = F, quote = F)
################################################################
test1 = chr12[which(chr12$pos < 44668179),]
test2 = chr12[which(chr12$pos > 43668179 & chr12$pos < 89336358),]
test3 = chr12[which(chr12$pos > 88336358 & chr12$pos <134004537),]

write.table(test1, "chr12_1.txt", row.names = F, quote = F)
write.table(test2, "chr12_2.txt", row.names = F, quote = F)
write.table(test3, "chr12_3.txt", row.names = F, quote = F)
################################################################
test1 = chr13[which(chr13$pos < 67065496),]
test2 = chr13[which(chr13$pos > 66065496 & chr13$pos < 134130993),]

write.table(test1, "chr13_1.txt", row.names = F, quote = F)
write.table(test2, "chr13_2.txt", row.names = F, quote = F)
################################################################
test1 = chr14[which(chr14$pos < 63294825),]
test2 = chr14[which(chr14$pos > 62294825 & chr14$pos <  126589651),]

write.table(test1, "chr14_1.txt", row.names = F, quote = F)
write.table(test2, "chr14_2.txt", row.names = F, quote = F)
################################################################
test1 = chr15[which(chr15$pos < 61275803),]
test2 = chr15[which(chr15$pos > 60275803 & chr15$pos <  122551607),]

write.table(test1, "chr15_1.txt", row.names = F, quote = F)
write.table(test2, "chr15_2.txt", row.names = F, quote = F)
################################################################
test1 = chr16[which(chr16$pos < 45185321),]
test2 = chr16[which(chr16$pos > 44185321 & chr16$pos <  90370643),]

write.table(test1, "chr16_1.txt", row.names = F, quote = F)
write.table(test2, "chr16_2.txt", row.names = F, quote = F)
################################################################
test1 = chr17[which(chr17$pos < 40554920),]
test2 = chr17[which(chr17$pos > 39554920 & chr17$pos <  81109841),]

write.table(test1, "chr17_1.txt", row.names = F, quote = F)
write.table(test2, "chr17_2.txt", row.names = F, quote = F)
################################################################
test1 = chr18[which(chr18$pos < 39024359),]
test2 = chr18[which(chr18$pos > 38024359 & chr18$pos <  78048719),]

write.table(test1, "chr18_1.txt", row.names = F, quote = F)
write.table(test2, "chr18_2.txt", row.names = F, quote = F)
################################################################
write.table(chr19, "chr19.txt", row.names = F, quote = F)
write.table(chr20, "chr20.txt", row.names = F, quote = F)
write.table(chr21, "chr21.txt", row.names = F, quote = F)
write.table(chr22, "chr22.txt", row.names = F, quote = F)}








