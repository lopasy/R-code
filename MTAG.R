### MTAG ###
library(data.table)

# Load data
predicted_exposed = fread("F:/Full_ukb/February/predicted_exposed.txt")
predicted_unexposed = fread("F:/Full_ukb/February/predicted_unexposed.txt")
true_all_exposed = fread("F:/Full_ukb/February/all/true_all_exposed.txt")
true_all_unexposed = fread("F:/Full_ukb/February/all/true_all_unexposed.txt")


# To get MAF and ref
predicted_ex = fread("F:/Full_ukb/February/Frequency/predicted_exposed.frq")
predicted_unex = fread("F:/Full_ukb/February/Frequency/predicted_unexposed.frq")
true_ex = fread("F:/Full_ukb/February/Frequency/true_exposed.frq")
true_unex = fread("F:/Full_ukb/February/Frequency/true_unexposed.frq")


predicted_ex = predicted_ex[,c(2,4,5)]
predicted_unex = predicted_unex[,c(2,4,5)]
true_ex = true_ex[,c(2,4,5)]
true_unex = true_unex[,c(2,4,5)]


predicted_exposed = merge(predicted_exposed, predicted_ex, "SNP")
predicted_unexposed = merge(predicted_unexposed, predicted_unex, "SNP")
true_exposed = merge(true_all_exposed, true_ex, "SNP")
true_unexposed = merge(true_all_unexposed, true_unex, "SNP")
rm(predicted_ex, predicted_unex, true_ex, true_unex)


# MTAG summary statistic format
# snpid    chr    bpos    a1    a2    freq    z    pval    n
predicted_exposed = predicted_exposed[,c(1,2,3,4,11,12,8,9,6)]
predicted_unexposed = predicted_unexposed[,c(1,2,3,4,11,12,8,9,6)]
true_exposed = true_exposed[,c(1,2,3,4,11,12,8,9,6)]
true_unexposed = true_unexposed[,c(1,2,3,4,11,12,8,9,6)]

names(predicted_exposed) = c("snpid",	"chr",	"bpos",	"a1",	"a2",	"freq",	"z",	"pval",	"n")
names(predicted_unexposed) = c("snpid",	"chr",	"bpos",	"a1",	"a2",	"freq",	"z",	"pval",	"n")
names(true_exposed) = c("snpid",	"chr",	"bpos",	"a1",	"a2",	"freq",	"z",	"pval",	"n")
names(true_unexposed) = c("snpid",	"chr",	"bpos",	"a1",	"a2",	"freq",	"z",	"pval",	"n")


write.table(predicted_exposed, "predicted_exposed.txt", row.names = F, quote = F)
write.table(predicted_unexposed, "predicted_unexposed.txt", row.names = F, quote = F)
write.table(true_exposed, "true_exposed.txt", row.names = F, quote = F)
write.table(true_unexposed, "true_unexposed.txt", row.names = F, quote = F)












