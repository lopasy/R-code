######################################
### Add asterisks to certain cells ###
######################################

chr$p_mr = -log10(chr$P_MR)
chr$p_cqr1 = -log10(chr$P_CQR1)
chr$p_cqr2 = -log10(chr$P_CQR2)

chr$p_mr = chr$p_mr/2.932
chr$p_cqr1 = chr$p_cqr1/2.230
chr$p_cqr2 = chr$p_cqr2/1.998

chr$p_mr = formatC(10^(-chr$p_mr), digits = 2, format = "e")
chr$p_cqr1 = formatC(10^(-chr$p_cqr1), digits = 2, format = "e")
chr$p_cqr2 = formatC(10^(-chr$p_cqr2), digits = 2, format = "e")

length(which(as.numeric(chr$p_mr) < 0.05/146))
length(which(as.numeric(chr$p_cqr1) < 0.05/146))
length(which(as.numeric(chr$p_cqr2) < 0.05/146))
             
for(i in 1:nrow(chr)){
  threshold = 0.05/nrow(chr)
  if(as.numeric(chr[i,14]) < threshold){
    chr[i,14] = paste(formatC(chr[i,14], digits = 2, format = "f"), "*", sep="")
  } else{
    chr[i,14] = paste(as.numeric(chr[i,14]))
  }
}

for(i in 1:nrow(chr)){
  threshold = 0.05/nrow(chr)
  if(as.numeric(chr[i,15]) < threshold){
    chr[i,15] = paste(formatC(chr[i,15], digits = 2, format = "f"), "*", sep="")
  } else{
    chr[i,15] = paste(as.numeric(chr[i,15]))
  }
}

for(i in 1:nrow(chr)){
  threshold = 0.05/nrow(chr)
  if(as.numeric(chr[i,16]) < threshold){
    chr[i,16] = paste(formatC(chr[i,16], digits = 2, format = "f"), "*", sep="")
  } else{
    chr[i,16] = paste(as.numeric(chr[i,16]))
  }
}


write.csv(chr,"chr.csv", row.names = F, quote = F)
