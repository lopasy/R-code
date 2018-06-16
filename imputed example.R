# library(data.table)
chr11 = fread("imp_maf_chr11.dat", header = T)
b = chr11[which(chr11$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr11snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr11)

# library(data.table)
chr12 = fread("imp_maf_chr12.dat", header = T)
b = chr12[which(chr12$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr12snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr12)

# library(data.table)
chr13 = fread("imp_maf_chr13.dat", header = T)
b = chr13[which(chr13$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr13snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr13)

# library(data.table)
chr14 = fread("imp_maf_chr14.dat", header = T)
b = chr14[which(chr14$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr14snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr14)

chr15 = fread("imp_maf_chr15.dat", header = T)
b = chr15[which(chr15$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr15snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr15)

chr16 = fread("imp_maf_chr16.dat", header = T)
b = chr16[which(chr16$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr16snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr16)

chr17 = fread("imp_maf_chr17.dat", header = T)
b = chr17[which(chr17$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr17snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr17)

chr18 = fread("imp_maf_chr18.dat", header = T)
b = chr18[which(chr18$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr18snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr18)

chr19 = fread("imp_maf_chr19.dat", header = T)
b = chr19[which(chr19$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr19snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr19)

chr20 = fread("imp_maf_chr20.dat", header = T)
b = chr20[which(chr20$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr20snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr20)

chr21 = fread("imp_maf_chr21.dat", header = T)
b = chr21[which(chr21$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr21snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr21)

chr22 = fread("imp_maf_chr22.dat", header = T)
b = chr22[which(chr22$information > 0.5),]

# Make a file with no dots
x = b[!(b$RSID == "."),]
# Create a list of rsIDs
x = as.data.frame(x)
x = x[,2]
x = as.list(x)
# Extract the dots
y = b[(b$RSID == "."),]
# Create a list of dots
y = as.data.frame(y)
y = y[,1]
y = as.list(y)

rm(b)

# Concatenate lists
c = c(x,y)
c = as.data.frame(c)
# Make a transpose
c = t(c)
c = as.data.frame(c)

# Save
write.table(c, file = "chr22snps_to_excude_info_0.5", quote = F, col.names = F, row.names = F)


# P.S. Advisible to rm some files and run gc()
# gc()
rm(chr22)