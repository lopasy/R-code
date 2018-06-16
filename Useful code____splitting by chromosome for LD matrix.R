# Extract SNPs for each CHR
chr = 1:22
for (i in chr){
  assign(paste0("chr", i), EduYearsHigh_avMSE[which(EduYearsHigh_avMSE$CHR == i),])
}

# Make a list
clist = list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,
          chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)


# Check if correct
{list = matrix(ncol = 22)
j = 1
for (i in clist){
  list[, j] = dim(i)[1]
  j = j + 1
}
print(sum(list))
}
  

# Split dataset by p-values and save files
EduYearsHigh_avMSE = merge(snp, EduYearsHigh_avMSE, by = "SNPID")
EduYearsHigh_avMSE = EduYearsHigh_avMSE[,c(3,9)]

plist = c(0.05,0.005,0.0005,0.00005,0.000005,0.0000005,0.00000005)
list = c()
j = 1
for (i in plist){
  list[[j]] = list(EduYearsHigh_avMSE[which(EduYearsHigh_avMSE$P < i),])
  filename <- paste("p_", i, ".txt", sep="")
  write.table(list[[j]], filename, sep = "\t", quote = F, row.names = F)
  j = j + 1
}
  


