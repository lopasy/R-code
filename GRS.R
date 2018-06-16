a = merge(recoded, jezs_unrelated, by = 1)
b = a[,c(7:152)];rm(a)
hits2 = hits[which(!(hits$`1000G_SNP` == "rs116226959")),] # rs116226959 rs36067445 rs75770582 rs144472942 rs80253120 rs201775671
hits2 = hits2[which(!(hits2$`1000G_SNP` == "rs201775671")),]
b = b[,-143]
hits3=hits2$ZSCORE
grs = sweep(as.matrix(b, na.rm=T),MARGIN = 2, as.vector(hits2$ZSCORE), `*`)
grs2 = as.data.frame(grs)
grs2$grs = rowSums(grs,na.rm = T)
weighted = (grs2$grs*145)/sum(grs2$grs)

a=merge(ukb,jezs_unrelated, by = 1)
a$grs=weighted
summary(lm(avMSE~Sex+Age+Array+UniEdu+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+UniEdu*grs,data=a))

z = cream_assoc[which(cream_assoc$TEST == "ADD"),]; z=z[-142,]; z=as.data.frame(z$BETA)


test = sweep(as.matrix(b, na.rm=T),MARGIN = 2, as.vector(z$z), `*`)
test = as.data.frame(test)
test$grs = rowSums(test,na.rm = T)
weights = (test$grs*138)/sum(test$grs)
a$grs2=weights
summary(lm(avMSE~Sex+Age+Array+UniEdu+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+UniEdu*grs2,data=a))

a$test = test$grs
summary(lm(avMSE~Sex+Age+Array+UniEdu+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+UniEdu*test,data=a))
