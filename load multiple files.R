# Case-control association
cd /scratch/share_PR300/Alfred/Assoc
module load R/3.3.0
R
library(data.table)
chr = list.files(pattern = ".linear", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
d = d[which(d$TEST == "ADD"),]

write.table(d, file ="/scratch/share_PR300/Alfred/Most_Important/maf0.01.txt", quote = F, row.names = F)


# Case-control interaction
cd /scratch/share_PR300/Alfred/GxEScan/Meta/2
module load R/3.3.0
R
library(data.table)
chr = list.files(pattern = ".linear", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
d = d[which(d$TEST == "ADDxUniEdu"),]

write.table(d, file ="/scratch/share_PR300/Alfred/GxEScan/Meta/1/meta2.txt", quote = F, row.names = F)



# Two-stage interaction with multiple covariates
cd /scratch/share_PR300/Alfred/Twostage
module load R/3.3.0
R
library(data.table)
chr = list.files(pattern = ".linear", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
d = d[which(d$TEST == "ADDxUniEdu"),]
d = d[order(d$P),]

write.table(d, file ="/scratch/share_PR300/Alfred/Most_Important/twostage.txt", quote = F, row.names = F)



# cc & case ge
cd /scratch/share_PR300/Alfred/BOLT/pp_1.5
module load R/3.3.0
R
library(data.table)
chr = list.files(pattern = ".CC_GE.gxeout", full.names = TRUE)
#chr=chr[c(1:100)]
chr = lapply(chr, fread)
#d = do.call(rbind.data.frame, chr)
#d = d[,-4]
#d = d[order(d$P),]
#write.table(d, file ="/home/c1669309/test_CC_GE.gxeout", sep = "\t", quote = F, row.names = F)

chr=lapply(chr, function(x) x[,-4])
for (ii in 1:100){
filename <- paste("pp1.5_set1_", ii, "_CC_GE.gxeout", sep="")
write.table(chr[[ii]], filename, sep = "\t", quote = F, row.names = F)
}



chr = list.files(pattern = ".Case_GE.gxeout", full.names = TRUE)
#chr=chr[c(101:176)]
chr = lapply(chr, fread)
#d = do.call(rbind.data.frame, chr)
#d = d[,-6]
#d = d[order(d$P),]
#write.table(d, file ="/home/c1669309/test_Case_GE.gxeout", sep = "\t", quote = F, row.names = F)

chr=lapply(chr, function(x) x[,-6])
for (ii in 1:100){
filename <- paste("pp1.5_set2_", ii, "_Case_GE.gxeout", sep="")
write.table(chr[[ii]], filename, sep = "\t", quote = F, row.names = F)
}


module load R/3.3.0
cd /scratch/share_PR300/Alfred/GRM/Chromosome1
R
library(readr)
chr = list.files(pattern = ".hsq", full.names = TRUE)
chr = lapply(chr, read.delim)
d = do.call(rbind.data.frame, chr)
d = d[which(d$Source == "V(G)/Vp"),]
d$Variance = as.numeric(d$Variance)
sum(d$Variance)

write.table(d, file ="t.hsq", quote = F, row.names = F)

