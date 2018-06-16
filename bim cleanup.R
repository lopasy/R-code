# library(data.table)
chr22 = fread("chr22.plink.bim")

# Make a file with no dots
x = chr22[!(chr22$V2 == "."),]

# Extract the dots
y = chr22[(chr22$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr22.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr21 = fread("chr21.plink.bim")

# Make a file with no dots
x = chr21[!(chr21$V2 == "."),]

# Extract the dots
y = chr21[(chr21$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr21.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr20 = fread("chr20.plink.bim")

# Make a file with no dots
x = chr20[!(chr20$V2 == "."),]

# Extract the dots
y = chr20[(chr20$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr20.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr19 = fread("chr19.plink.bim")

# Make a file with no dots
x = chr19[!(chr19$V2 == "."),]

# Extract the dots
y = chr19[(chr19$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr19.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr18 = fread("chr181.plink.bim")

# Make a file with no dots
x = chr18[!(chr18$V2 == "."),]

# Extract the dots
y = chr18[(chr18$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr18.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr17 = fread("chr17.plink.bim")

# Make a file with no dots
x = chr17[!(chr17$V2 == "."),]

# Extract the dots
y = chr17[(chr17$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr17.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr16 = fread("chr16.plink.bim")

# Make a file with no dots
x = chr16[!(chr16$V2 == "."),]

# Extract the dots
y = chr16[(chr16$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr16.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr15 = fread("chr15.plink.bim")

# Make a file with no dots
x = chr15[!(chr15$V2 == "."),]

# Extract the dots
y = chr15[(chr15$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr15.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr1 = fread("chr1.plink.bim")

# Make a file with no dots
x = chr1[!(chr1$V2 == "."),]

# Extract the dots
y = chr1[(chr1$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr1.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr14 = fread("chr14.plink.bim")

# Make a file with no dots
x = chr14[!(chr14$V2 == "."),]

# Extract the dots
y = chr14[(chr14$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr14.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr13 = fread("chr13.plink.bim")

# Make a file with no dots
x = chr13[!(chr13$V2 == "."),]

# Extract the dots
y = chr13[(chr13$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr13.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr12 = fread("chr12.plink.bim")

# Make a file with no dots
x = chr12[!(chr12$V2 == "."),]

# Extract the dots
y = chr12[(chr12$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr12.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr11 = fread("chr11.plink.bim")

# Make a file with no dots
x = chr11[!(chr11$V2 == "."),]

# Extract the dots
y = chr11[(chr11$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr11.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr10 = fread("chr10.plink.bim")

# Make a file with no dots
x = chr10[!(chr10$V2 == "."),]

# Extract the dots
y = chr10[(chr10$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr10.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr9 = fread("chr9.plink.bim")

# Make a file with no dots
x = chr9[!(chr9$V2 == "."),]

# Extract the dots
y = chr9[(chr9$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr9.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr8 = fread("chr8.plink.bim")

# Make a file with no dots
x = chr8[!(chr8$V2 == "."),]

# Extract the dots
y = chr8[(chr8$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr8.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr7 = fread("chr7.plink.bim")

# Make a file with no dots
x = chr7[!(chr7$V2 == "."),]

# Extract the dots
y = chr7[(chr7$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr7.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr6 = fread("chr6.plink.bim")

# Make a file with no dots
x = chr6[!(chr6$V2 == "."),]

# Extract the dots
y = chr6[(chr6$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr6.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr5 = fread("chr5.plink.bim")

# Make a file with no dots
x = chr5[!(chr5$V2 == "."),]

# Extract the dots
y = chr5[(chr5$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr5.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr4 = fread("chr4.plink.bim")

# Make a file with no dots
x = chr4[!(chr4$V2 == "."),]

# Extract the dots
y = chr4[(chr4$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr4.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr3 = fread("chr3.plink.bim")

# Make a file with no dots
x = chr3[!(chr3$V2 == "."),]

# Extract the dots
y = chr3[(chr3$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr3.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

chr2 = fread("chr2.plink.bim")

# Make a file with no dots
x = chr2[!(chr2$V2 == "."),]

# Extract the dots
y = chr2[(chr2$V2 == "."),]

# Create a column vector for dots replacement
a = paste(y$V1, y$V4, sep = ":")
a = paste(a, y$V6, y$V5, sep = "_")
a = as.data.frame(a)

# Change dot column with chrom:pos_allele
y[,2] = a[,1]
c = rbind(x,y)

write.table(c, file = "chr2.plink.bim", quote = F, col.names = F, row.names = F)

rm(list=ls())

