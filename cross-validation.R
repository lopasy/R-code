{a=as.data.frame(sample(predictedMSE_ID$FID,26705,replace = F))
z = predictedMSE_ID[which(!(predictedMSE_ID$FID %in% a$`sample(predictedMSE_ID$FID, 26705, replace = F)`)),]
b=as.data.frame(sample(z$FID,26705,replace = F))
z = z[which(!(z$FID %in% b$`sample(z$FID, 26705, replace = F)`)),]
c=as.data.frame(sample(z$FID,26705,replace = F))
z = z[which(!(z$FID %in% c$`sample(z$FID, 26705, replace = F)`)),]
d=as.data.frame(sample(z$FID,26705,replace = F))
z = z[which(!(z$FID %in% d$`sample(z$FID, 26705, replace = F)`)),]
e=as.data.frame(sample(z$FID,26705,replace = F))
z = z[which(!(z$FID %in% e$`sample(z$FID, 26705, replace = F)`)),]
f = z
rm(z)
a[,2]=a[,1]
b[,2]=b[,1]
c[,2]=c[,1]
d[,2]=d[,1]
e[,2]=e[,1]
write.table(a,file = "gwas1.txt", col.names = F, row.names = F, quote = F)
write.table(b,file = "gwas2.txt", col.names = F, row.names = F, quote = F)
write.table(c,file = "gwas3.txt", col.names = F, row.names = F, quote = F)
write.table(d,file = "gwas4.txt", col.names = F, row.names = F, quote = F)
write.table(e,file = "gwas5.txt", col.names = F, row.names = F, quote = F)
write.table(f,file = "gwas6.txt", col.names = F, row.names = F, quote = F)}
