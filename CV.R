library(data.table)
chr = list.files(pattern = "predicted_3.", full.names = TRUE)
chr = lapply(chr, fread)



as = 1:20
for (i in as){
  assign(paste0("chr", i), as.data.frame(chr[c(i)]))
}


colnames(bim)=c("CHR","SNP","pos","BP","A2","A1")
colnames(chr1)[2]=c("SNP")
chr1=merge(chr1,bim,by="SNP")
colnames(chr2)[2]=c("SNP")
chr2=merge(chr2,bim,by="SNP")
colnames(chr3)[2]=c("SNP")
chr3=merge(chr3,bim,by="SNP")
colnames(chr4)[2]=c("SNP")
chr4=merge(chr4,bim,by="SNP")
colnames(chr5)[2]=c("SNP")
chr5=merge(chr5,bim,by="SNP")
colnames(chr6)[2]=c("SNP")
chr6=merge(chr6,bim,by="SNP")

chr1=chr1[,c(1,2,13,14,6,7,9)]
colnames(chr1)=c("SNP","CHR","A2","A1","NMISS","OR","P")
chr2=chr2[,c(1,2,13,14,6,7,9)]
colnames(chr2)=c("SNP","CHR","A2","A1","NMISS","OR","P")
chr3=chr3[,c(1,2,13,14,6,7,9)]
colnames(chr3)=c("SNP","CHR","A2","A1","NMISS","OR","P")
chr4=chr4[,c(1,2,13,14,6,7,9)]
colnames(chr4)=c("SNP","CHR","A2","A1","NMISS","OR","P")
chr5=chr5[,c(1,2,13,14,6,7,9)]
colnames(chr5)=c("SNP","CHR","A2","A1","NMISS","OR","P")
chr6=chr6[,c(1,2,13,14,6,7,9)]
colnames(chr6)=c("SNP","CHR","A2","A1","NMISS","OR","P")


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
  chr20 = chr20[which(chr20$TEST == "ADDxUniEdu"),]}

{chr = rbind(chr1,chr2)
chr = rbind(chr,chr3)
chr = rbind(chr,chr4)
chr = rbind(chr,chr5)
chr = rbind(chr,chr6)
chr = rbind(chr,chr7)
chr = rbind(chr,chr8)
chr = rbind(chr,chr9)
chr = rbind(chr,chr10)
chr = rbind(chr,chr11)
chr = rbind(chr,chr12)
chr = rbind(chr,chr13)
chr = rbind(chr,chr14)
chr = rbind(chr,chr15)
chr = rbind(chr,chr16)
chr = rbind(chr,chr17)
chr = rbind(chr,chr18)
chr = rbind(chr,chr19)
chr = rbind(chr,chr20)

write.table(chr, "predicted3.txt", row.names = F, quote = F)}


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



write.table(chr1,"gwas1_1.5.txt",quote = F, row.names = F)
write.table(chr2,"gwas2_1.5.txt",quote = F, row.names = F)
write.table(chr3,"gwas3_1.5.txt",quote = F, row.names = F)
write.table(chr4,"gwas4_1.5.txt",quote = F, row.names = F)
write.table(chr5,"gwas5_1.5.txt",quote = F, row.names = F)
write.table(chr6,"gwas6_1.5.txt",quote = F, row.names = F)



chr1=chr1[,c(1,6)]
chr2=chr2[,c(1,6)]
chr3=chr3[,c(1,6)]
chr4=chr4[,c(1,6)]
chr5=chr5[,c(1,6)]
chr6=chr6[,c(1,6)]

colnames(chr1)=c("SNP","P")
colnames(chr2)=c("SNP","P")
colnames(chr3)=c("SNP","P")
colnames(chr4)=c("SNP","P")
colnames(chr5)=c("SNP","P")
colnames(chr6)=c("SNP","P")

write.table(chr1,"meta1_1.5.txt",quote = F, row.names = F)
write.table(chr2,"meta2_1.5.txt",quote = F, row.names = F)
write.table(chr3,"meta3_1.5.txt",quote = F, row.names = F)
write.table(chr4,"meta4_1.5.txt",quote = F, row.names = F)
write.table(chr5,"meta5_1.5.txt",quote = F, row.names = F)
write.table(chr6,"meta6_1.5.txt",quote = F, row.names = F)


