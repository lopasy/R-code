library(flexmix); library(ggplot2)

ukb = ukb[which(ukb$het_miss_out == 0),]
ukb = ukb[which(ukb$in_kinship == 0),]
ukb2 = ukb[!is.na(ukb$AgeSpexWear),]
ukb2 = ukb2[!is.na(ukb2$UniEdu),]

mo1 = FLXMRglm(family = "gaussian")
mo2 = FLXMRglm(family = "gaussian")

flexfit = flexmix(AgeSpexWear ~ 1, data = ukb2, k = 2, model = list(mo1, mo2))

c1 = parameters(flexfit, component=1)[[1]]
c2 = parameters(flexfit, component=2)[[1]]

cat('pred mean:', c1[1], '\n'); cat('pred sd:', c1[2], '\n')
cat('pred mean:', c2[1], '\n'); cat('pred sd:', c2[2], '\n')

plot(flexfit)

plot_mix_comps = function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

lam = table(clusters(flexfit))

ggplot(ukb2) +
  geom_histogram(aes(AgeSpexWear, ..density..), binwidth = 1, colour = "black", fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(c1[1], c1[2], lam[1]/sum(lam)),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(c2[1], c2[2], lam[2]/sum(lam)),
                colour = "blue", lwd = 1.5) +
  ylab("Density")


mod2 = summary(refit(flexfit))

wh_mix <- stepFlexmix(AgeSpexWear ~ 1, data = ukb2, control = list(verbose = 0), k=1:5, nrep=10)
plot(BIC(wh_mix),type='b',ylab='BIC')
points(x = which.min(BIC(wh_mix)),min(BIC(wh_mix)),col='red',pch=20)

wh_best <- getModel(wh_mix,'BIC')
print(wh_best)

##############################################################################################################
#CNV analysis
cnvs = 47:50
results = as.data.frame(matrix(nrow = length(cnvs), ncol = 9))
names(results) = c("CNV", "c1_b", "c1_lci", "c1_uci", "c1_p", "c2_b", "c2_lci", "c2_uci", "c2_p")
for(i in cnvs){
  ukb3 = ukb2[!is.na(ukb2[,i]),]
  flexfit = flexmix(AgeSpexWear ~ Age + Sex + UniEdu + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Geno_array + ukb3[,i], 
                    data = ukb3, k = 1, model = list(mo1, mo2))
  mod2 = summary(refit(flexfit))
  results[i-46,1] = names(ukb3)[i]
  comp1 = summary(mod2)@components@.Data[[1]][1]; comp1 = as.data.frame(comp1)
  comp2 = summary(mod2)@components@.Data[[1]][2]; comp2 = as.data.frame(comp2)
  
  results[i-46,2] = comp1[2,1]
  results[i-46,3] = comp1[2,1] - 1.96*comp1[2,2]
  results[i-46,4] = comp1[2,1] + 1.96*comp1[2,2]
  results[i-46,5] = comp1[2,4]
  
  results[i-46,6] = comp2[2,1]
  results[i-46,7] = comp2[2,1] - 1.96*comp2[2,2]
  results[i-46,8] = comp2[2,1] + 1.96*comp2[2,2]
  results[i-46,9] = comp2[2,4]
}

#matplot(ukb3$Age, fitted(flexfit), pch = 16, type = "p")
#points(ukb3$ukb3[,i], NPreg$AgeSpexWear)
