ukb = merge(ukb, YG, by = "SNPID")

z = rbind[which(rbind$CHR == 12),]

z = ukb[which(ukb$SNP == "rs2008065"),]



ukbA = merge(chr10, chr10S, by = "SNPID")
ukbB = merge(chr11, chr11S, by = "SNPID")
ukbC = merge(chr12, chr12S, by = "SNPID")
rbind = rbind(rbind,ukbC)


write.table(rbind, file = "ukb", row.names = F, quote = F)




final$new = paste(final$X1,final$X2, sep = ":")
final$new = paste(final$new,final$X4,final$X5, sep = "_")
rbind$new = paste(rbind$CHR,rbind$BP,rbind$SNP, sep = "_")

z = final[which(final$X3 != "."),]
zz = final[which(final$X3 == "."),]
z = z[,3]
colnames(z)[1] = "new"
zz = zz[,6]
zzz = rbind(z,zz)
zzzz = rbind(zzz,d)
write.table(zzzz, file = "snps.txt", row.names = F, quote = F, col.names = F)




z = merge(rbind,final, by = "new")




library(data.table)
chr = list.files(pattern = ".bim", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
d = d[,2]
colnames(d)[1] = "new"
write.table(d, file ="dir_gen.txt", quote = F, row.names = F)

