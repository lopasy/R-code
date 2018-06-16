library(data.table)
chr = list.files(pattern = ".linear", full.names = TRUE)
chr = lapply(chr, fread)

chr = chr[c(5,9)]



as = 1:32
for (i in as){
  assign(paste0("chr", i), as.data.frame(chr[c(i)]))
}



{chr1 = chr1[which(chr1$TEST == "ADDxUniEdu"),]
chr2 = chr2[which(chr2$TEST == "ADDxUniEdu"),]
chr3 = chr3[which(chr3$TEST == "ADDxUniEdu"),]
chr4 = chr4[which(chr4$TEST == "ADDxUniEdu"),]
chr5 = chr5[which(chr5$TEST == "ADDxUniEdu"),]
chr6 = chr6[which(chr6$TEST == "ADDxUniEdu"),]
chr7 = chr7[which(chr7$TEST == "ADDxUniEdu"),]
chr8 = chr8[which(chr8$TEST == "ADDxUniEdu"),]
chr9 = chr9[which(chr9$TEST == "ADDxUniEdu"),]
chr10 = chr10[which(chr10$TEST == "ADDxUniEdu"),]
chr11 = chr11[which(chr11$TEST == "ADDxUniEdu"),]
chr12 = chr12[which(chr12$TEST == "ADDxUniEdu"),]
chr13 = chr13[which(chr13$TEST == "ADDxUniEdu"),]
chr14 = chr14[which(chr14$TEST == "ADDxUniEdu"),]
chr15 = chr15[which(chr15$TEST == "ADDxUniEdu"),]
chr16 = chr16[which(chr16$TEST == "ADDxUniEdu"),]
chr17 = chr17[which(chr17$TEST == "ADDxUniEdu"),]
chr18 = chr18[which(chr18$TEST == "ADDxUniEdu"),]
chr19 = chr19[which(chr19$TEST == "ADDxUniEdu"),]
chr20 = chr20[which(chr20$TEST == "ADDxUniEdu"),]
chr21 = chr21[which(chr21$TEST == "ADDxUniEdu"),]
chr22 = chr22[which(chr22$TEST == "ADDxUniEdu"),]
chr23 = chr23[which(chr23$TEST == "ADDxUniEdu"),]
chr24 = chr24[which(chr24$TEST == "ADDxUniEdu"),]
chr25 = chr25[which(chr25$TEST == "ADDxUniEdu"),]
chr26 = chr26[which(chr26$TEST == "ADDxUniEdu"),]
chr27 = chr27[which(chr27$TEST == "ADDxUniEdu"),]
chr28 = chr28[which(chr28$TEST == "ADDxUniEdu"),]
chr29 = chr29[which(chr29$TEST == "ADDxUniEdu"),]
chr30 = chr30[which(chr30$TEST == "ADDxUniEdu"),]
chr31 = chr31[which(chr31$TEST == "ADDxUniEdu"),]
chr32 = chr32[which(chr32$TEST == "ADDxUniEdu"),]}


{chr1 = chr1[which(chr1$P < 0.05),]
  chr2 = chr2[which(chr2$P < 0.05),]
  chr3 = chr3[which(chr3$P < 0.05),]
  chr4 = chr4[which(chr4$P < 0.05),]
  chr5 = chr5[which(chr5$P < 0.05),]
  chr6 = chr6[which(chr6$P < 0.05),]
  chr7 = chr7[which(chr7$P < 0.05),]
  chr8 = chr8[which(chr8$P < 0.05),]
  chr9 = chr9[which(chr9$P < 0.05),]
  chr10 = chr10[which(chr10$P < 0.05),]
  chr11 = chr11[which(chr11$P < 0.05),]
  chr12 = chr12[which(chr12$P < 0.05),]
  chr13 = chr13[which(chr13$P < 0.05),]
  chr14 = chr14[which(chr14$P < 0.05),]
  chr15 = chr15[which(chr15$P < 0.05),]
  chr16 = chr16[which(chr16$P < 0.05),]
  chr17 = chr17[which(chr17$P < 0.05),]
  chr18 = chr18[which(chr18$P < 0.05),]
  chr19 = chr19[which(chr19$P < 0.05),]
  chr20 = chr20[which(chr20$P < 0.05),]
  chr21 = chr21[which(chr21$P < 0.05),]
  chr22 = chr22[which(chr22$P < 0.05),]
  chr23 = chr23[which(chr23$P < 0.05),]
  chr24 = chr24[which(chr24$P < 0.05),]
  chr25 = chr25[which(chr25$P < 0.05),]
  chr26 = chr26[which(chr26$P < 0.05),]
  chr27 = chr27[which(chr27$P < 0.05),]
  chr28 = chr28[which(chr28$P < 0.05),]
  chr29 = chr29[which(chr29$P < 0.05),]
  chr30 = chr30[which(chr30$P < 0.05),]
  chr31 = chr31[which(chr31$P < 0.05),]
  chr32 = chr32[which(chr32$P < 0.05),]}

