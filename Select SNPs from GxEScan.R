test <- read_delim("~/Downloads/GxEScan_0_5_0/new.snpinfo", "\t", escape_double = FALSE, trim_ws = TRUE)
testYG <- read_delim("~/Downloads/GxEScan_0_5_0/new_QT_YG.gxeout", "\t", escape_double = FALSE, trim_ws = TRUE)
testVH <- read_delim("~/Downloads/GxEScan_0_5_0/new_QT_VH.gxeout", "\t", escape_double = FALSE, trim_ws = TRUE)
testYGVH <- read_delim("~/Downloads/GxEScan_0_5_0/new_QT_YGVH.gxeout", "\t", escape_double = FALSE, trim_ws = TRUE)


a = merge(test, testYG, by = "SNPID")
a = merge(test, testVH, by = "SNPID")
a = merge(test, testYGVH, by = "SNPID")

b = a[which(a$P < 0.05),]
b = b$SNP

c = a[which(a$P < 0.005),]
c = c$SNP

d = a[which(a$P < 0.0005),]
d = d$SNP

e = a[which(a$P < 0.00005),]
e = e$SNP

f = a[which(a$P < 0.000005),]
f = f$SNP

g = a[which(a$P < 0.0000005),]
g = g$SNP

h = a[which(a$P < 0.00000005),]
h = h$SNP


setwd("~/Downloads/GENS2")
write.table(b, file = "SNP_0.05YGVH.txt", col.names = F, row.names = F, quote = F)
write.table(c, file = "SNP_5e-03YGVH.txt", col.names = F, row.names = F, quote = F)
write.table(d, file = "SNP_5e-04YGVH.txt", col.names = F, row.names = F, quote = F)
write.table(e, file = "SNP_5e-05YGVH.txt", col.names = F, row.names = F, quote = F)
write.table(f, file = "SNP_5e-06YGVH.txt", col.names = F, row.names = F, quote = F)
write.table(g, file = "SNP_5e-07YGVH.txt", col.names = F, row.names = F, quote = F)
write.table(h, file = "SNP_5e-08YGVH.txt", col.names = F, row.names = F, quote = F)
