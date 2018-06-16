meta6_4.39 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.000005),]
meta6_4.52 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.000004),]
meta6_4.69 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.00002),]
meta6_5 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.00001),]
meta6_5.39 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.000004),]
meta6_5.52 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.000003),]
meta6_5.69 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.000002),]
meta6_6 = METAANALYSIS1[which(METAANALYSIS1$`P-value` < 0.000001),]




colnames(meta6_4.39) = c("SNP","P");colnames(meta6_4.52) = c("SNP","P");colnames(meta6_4.69) = c("SNP","P");colnames(meta6_5) = c("SNP","P")
colnames(meta6_5.39) = c("SNP","P");colnames(meta6_5.52) = c("SNP","P");colnames(meta6_5.69) = c("SNP","P");colnames(meta6_6) = c("SNP","P")



write.table(meta6_4.39, file = "meta6_4.39.txt", quote = F, row.names = F);write.table(meta6_4.52, file = "meta6_4.52.txt", quote = F, row.names = F)
write.table(meta6_4.69, file = "meta6_4.69.txt", quote = F, row.names = F);write.table(meta6_5, file = "meta6_5.txt", quote = F, row.names = F)
write.table(meta6_5.39, file = "meta6_5.39.txt", quote = F, row.names = F);write.table(meta6_5.52, file = "meta6_5.52.txt", quote = F, row.names = F)
write.table(meta6_5.69, file = "meta6_5.69.txt", quote = F, row.names = F);write.table(meta6_6, file = "meta6_6.txt", quote = F, row.names = F)


rm(meta6_4.39,meta6_4.52,meta6_4.69,meta6_5,meta6_5.39,meta6_5.52,meta6_5.69,meta6_6)






colnames(METAANALYSIS6) = c("SNP","P")
write.table(METAANALYSIS1, file = "METAANALYSIS1.txt", quote = F, row.names = F);write.table(METAANALYSIS2, file = "METAANALYSIS2.txt", quote = F, row.names = F)
write.table(METAANALYSIS3, file = "METAANALYSIS3.txt", quote = F, row.names = F);write.table(METAANALYSIS4, file = "METAANALYSIS4.txt", quote = F, row.names = F)
write.table(METAANALYSIS5, file = "METAANALYSIS5.txt", quote = F, row.names = F);write.table(METAANALYSIS6, file = "METAANALYSIS6.txt", quote = F, row.names = F)
