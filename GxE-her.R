library(data.table)

real = fread("D:/AoO.bim")
bim = fread("D:/AoO_real.bim")
bim$V7 = real$V2
names(bim)[7] = "SNP"


num_chr = 22
interval = 1000000

chr = list()
segments = as.data.frame(matrix(nrow = num_chr, ncol = 1))
for(i in 1:num_chr){
chr[[i]] = bim[which(bim$V1 == i),]
a = round((max(bim[which(bim$V1 == i),4])-min(bim[which(bim$V1 == i),4]))/interval,0)+1
segments[i,1] = a
}

for(chrs in 1:num_chr){
    group1 = as.data.frame(matrix(nrow = 2000, ncol = (segments[chrs,])))
    region = chr[[chrs]]
    range = min(region$V4):(min(region$V4)+interval); range = as.data.frame(range)
    
          for(j in 1:ncol(group1)){
              len = region[which(region$V4 %in% range$range),2]
              
                if(isTRUE(nrow(len) > 1)){
                group1[1:nrow(len),j] = len
                }
              
              range = range + interval
              group2 = group1[,j]
              sav = group2[!is.na(group2)]
              write.table(sav, paste0("D:/Motivation/GxE-her/chr", chrs, "_", j), row.names = F, col.names = F, quote = F)
              
          }
    
}



