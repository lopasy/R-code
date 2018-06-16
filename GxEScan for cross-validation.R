library(data.table)
chr = list.files(pattern = ".snpinfo", full.names = TRUE)
chr = lapply(chr, fread)

gwas = 1:6
for (i in gwas){
  assign(paste0("gwas", i), chr[[i]])
}
rm(chr,i,gwas)

{gwas1_YG = fread("gwas1_QT_YG.gxeout", "\t")
gwas1_GxE = fread("gwas1_QT_GxE.gxeout", "\t")
gwas1_VH = fread("gwas1_QT_VH.gxeout", "\t")
gwas1_YGVH = fread("gwas1_QT_YGVH.gxeout", "\t")
gwas2_YG = fread("gwas2_QT_YG.gxeout", "\t")
gwas2_GxE = fread("gwas2_QT_GxE.gxeout", "\t")
gwas2_VH = fread("gwas2_QT_VH.gxeout", "\t")
gwas2_YGVH = fread("gwas2_QT_YGVH.gxeout", "\t")
gwas3_YG = fread("gwas3_QT_YG.gxeout", "\t")
gwas3_GxE = fread("gwas3_QT_GxE.gxeout", "\t")
gwas3_VH = fread("gwas3_QT_VH.gxeout", "\t")
gwas3_YGVH = fread("gwas3_QT_YGVH.gxeout", "\t")
gwas4_YG = fread("gwas4_QT_YG.gxeout", "\t")
gwas4_GxE = fread("gwas4_QT_GxE.gxeout", "\t")
gwas4_VH = fread("gwas4_QT_VH.gxeout", "\t")
gwas4_YGVH = fread("gwas4_QT_YGVH.gxeout", "\t")
gwas5_YG = fread("gwas5_QT_YG.gxeout", "\t")
gwas5_GxE = fread("gwas5_QT_GxE.gxeout", "\t")
gwas5_VH = fread("gwas5_QT_VH.gxeout", "\t")
gwas5_YGVH = fread("gwas5_QT_YGVH.gxeout", "\t")
gwas6_YG = fread("gwas6_QT_YG.gxeout", "\t")
gwas6_GxE = fread("gwas6_QT_GxE.gxeout", "\t")
gwas6_VH = fread("gwas6_QT_VH.gxeout", "\t")
gwas6_YGVH = fread("gwas6_QT_YGVH.gxeout", "\t")}




{cv1_YG = merge(gwas1, gwas1_YG, "SNPID")
cv1_GxE = merge(gwas1, gwas1_GxE, "SNPID")
cv1_VH = merge(gwas1, gwas1_VH, "SNPID")
cv1_YGVH = merge(gwas1, gwas1_YGVH, "SNPID")

cv2_YG = merge(gwas2, gwas2_YG, "SNPID")
cv2_GxE = merge(gwas2, gwas2_GxE, "SNPID")
cv2_VH = merge(gwas2, gwas2_VH, "SNPID")
cv2_YGVH = merge(gwas2, gwas2_YGVH, "SNPID")

cv3_YG = merge(gwas3, gwas3_YG, "SNPID")
cv3_GxE = merge(gwas3, gwas3_GxE, "SNPID")
cv3_VH = merge(gwas3, gwas3_VH, "SNPID")
cv3_YGVH = merge(gwas3, gwas3_YGVH, "SNPID")

cv4_YG = merge(gwas4, gwas4_YG, "SNPID")
cv4_GxE = merge(gwas4, gwas4_GxE, "SNPID")
cv4_VH = merge(gwas4, gwas4_VH, "SNPID")
cv4_YGVH = merge(gwas4, gwas4_YGVH, "SNPID")

cv5_YG = merge(gwas5, gwas5_YG, "SNPID")
cv5_GxE = merge(gwas5, gwas5_GxE, "SNPID")
cv5_VH = merge(gwas5, gwas5_VH, "SNPID")
cv5_YGVH = merge(gwas5, gwas5_YGVH, "SNPID")

cv6_YG = merge(gwas6, gwas6_YG, "SNPID")
cv6_GxE = merge(gwas6, gwas6_GxE, "SNPID")
cv6_VH = merge(gwas6, gwas6_VH, "SNPID")
cv6_YGVH = merge(gwas6, gwas6_YGVH, "SNPID")}

#rm(list=ls(pattern="gwas"))

{cv1_YG = cbind(cv1_YG, gwas[,6], "SNPID")
  cv1_GxE = cbind(cv1_GxE, gwas[,6], "SNPID")
  cv1_VH = cbind(cv1_VH, gwas[,6], "SNPID")
  cv1_YGVH = cbind(cv1_YGVH, gwas[,6], "SNPID")
  
  cv2_YG = cbind(cv2_YG, gwas[,6], "SNPID")
  cv2_GxE = cbind(cv2_GxE, gwas[,6], "SNPID")
  cv2_VH = cbind(cv2_VH, gwas[,6], "SNPID")
  cv2_YGVH = cbind(cv2_YGVH, gwas[,6], "SNPID")
  
  cv3_YG = cbind(cv3_YG, gwas[,6], "SNPID")
  cv3_GxE = cbind(cv3_GxE, gwas[,6], "SNPID")
  cv3_VH = cbind(cv3_VH, gwas[,6], "SNPID")
  cv3_YGVH = cbind(cv3_YGVH, gwas[,6], "SNPID")
  
  cv4_YG = cbind(cv4_YG, gwas[,6], "SNPID")
  cv4_GxE = cbind(cv4_GxE, gwas[,6], "SNPID")
  cv4_VH = cbind(cv4_VH, gwas[,6], "SNPID")
  cv4_YGVH = cbind(cv4_YGVH, gwas[,6], "SNPID")
  
  cv5_YG = cbind(cv5_YG, gwas[,6], "SNPID")
  cv5_GxE = cbind(cv5_GxE, gwas[,6], "SNPID")
  cv5_VH = cbind(cv5_VH, gwas[,6], "SNPID")
  cv5_YGVH = cbind(cv5_YGVH, gwas[,6], "SNPID")
  
  cv6_YG = cbind(cv6_YG, gwas[,6], "SNPID")
  cv6_GxE = cbind(cv6_GxE, gwas[,6], "SNPID")
  cv6_VH = cbind(cv6_VH, gwas[,6], "SNPID")
  cv6_YGVH = cbind(cv6_YGVH, gwas[,6], "SNPID")}

{cv1_YG = cv1_YG[,c(3,9)]
  cv1_GxE = cv1_GxE[,c(3,9)]
  cv1_VH = cv1_VH[,c(3,7)]
  cv1_YGVH = cv1_YGVH[,c(3,7)]
  
  cv2_YG = cv2_YG[,c(3,9)]
  cv2_GxE = cv2_GxE[,c(3,9)]
  cv2_VH = cv2_VH[,c(3,7)]
  cv2_YGVH = cv2_YGVH[,c(3,7)]
  
  cv3_YG = cv3_YG[,c(3,9)]
  cv3_GxE = cv3_GxE[,c(3,9)]
  cv3_VH = cv3_VH[,c(3,7)]
  cv3_YGVH = cv3_YGVH[,c(3,7)]
  
  cv4_YG = cv4_YG[,c(3,9)]
  cv4_GxE = cv4_GxE[,c(3,9)]
  cv4_VH = cv4_VH[,c(3,7)]
  cv4_YGVH = cv4_YGVH[,c(3,7)]
  
  cv5_YG = cv5_YG[,c(3,9)]
  cv5_GxE = cv5_GxE[,c(3,9)]
  cv5_VH = cv5_VH[,c(3,7)]
  cv5_YGVH = cv5_YGVH[,c(3,7)]
  
  cv6_YG = cv6_YG[,c(3,9)]
  cv6_GxE = cv6_GxE[,c(3,9)]
  cv6_VH = cv6_VH[,c(3,7)]
  cv6_YGVH = cv6_YGVH[,c(3,7)]}

{cv1_YG = cv1_YG[,c(3,5:7,9,10)]
  cv1_GxE = cv1_GxE[,c(3,5:7,9,10)]
  cv1_VH = cv1_VH[,c(2:5,7,8)]
  cv1_YGVH = cv1_YGVH[,c(2:5,7,8)]
  
  cv2_YG = cv2_YG[,c(3,5:7,9,10)]
  cv2_GxE = cv2_GxE[,c(3,5:7,9,10)]
  cv2_VH = cv2_VH[,c(2:5,7,8)]
  cv2_YGVH = cv2_YGVH[,c(2:5,7,8)]
  
  cv3_YG = cv3_YG[,c(3,5:7,9,10)]
  cv3_GxE = cv3_GxE[,c(3,5:7,9,10)]
  cv3_VH = cv3_VH[,c(2:5,7,8)]
  cv3_YGVH = cv3_YGVH[,c(2:5,7,8)]
  
  cv4_YG = cv4_YG[,c(3,5:7,9,10)]
  cv4_GxE = cv4_GxE[,c(3,5:7,9,10)]
  cv4_VH = cv4_VH[,c(2:5,7,8)]
  cv4_YGVH = cv4_YGVH[,c(2:5,7,8)]
  
  cv5_YG = cv5_YG[,c(3,5:7,9,10)]
  cv5_GxE = cv5_GxE[,c(3,5:7,9,10)]
  cv5_VH = cv5_VH[,c(2:5,7,8)]
  cv5_YGVH = cv5_YGVH[,c(2:5,7,8)]
  
  cv6_YG = cv6_YG[,c(3,5:7,9,10)]
  cv6_GxE = cv6_GxE[,c(3,5:7,9,10)]
  cv6_VH = cv6_VH[,c(2:5,7,8)]
  cv6_YGVH = cv6_YGVH[,c(2:5,7,8)]}



{write.table(cv1_YG, "scv1_YG.txt", row.names = F, quote = F)
write.table(cv1_GxE, "scv1_GxE.txt", row.names = F, quote = F)
write.table(cv1_VH, "scv1_VH.txt", row.names = F, quote = F)
write.table(cv1_YGVH, "scv1_YGVH.txt", row.names = F, quote = F)

write.table(cv2_YG, "scv2_YG.txt", row.names = F, quote = F)
write.table(cv2_GxE, "scv2_GxE.txt", row.names = F, quote = F)
write.table(cv2_VH, "scv2_VH.txt", row.names = F, quote = F)
write.table(cv2_YGVH, "scv2_YGVH.txt", row.names = F, quote = F)

write.table(cv3_YG, "scv3_YG.txt", row.names = F, quote = F)
write.table(cv3_GxE, "scv3_GxE.txt", row.names = F, quote = F)
write.table(cv3_VH, "scv3_VH.txt", row.names = F, quote = F)
write.table(cv3_YGVH, "scv3_YGVH.txt", row.names = F, quote = F)

write.table(cv4_YG, "scv4_YG.txt", row.names = F, quote = F)
write.table(cv4_GxE, "scv4_GxE.txt", row.names = F, quote = F)
write.table(cv4_VH, "scv4_VH.txt", row.names = F, quote = F)
write.table(cv4_YGVH, "scv4_YGVH.txt", row.names = F, quote = F)

write.table(cv5_YG, "scv5_YG.txt", row.names = F, quote = F)
write.table(cv5_GxE, "scv5_GxE.txt", row.names = F, quote = F)
write.table(cv5_VH, "scv5_VH.txt", row.names = F, quote = F)
write.table(cv5_YGVH, "scv5_YGVH.txt", row.names = F, quote = F)

write.table(cv6_YG, "scv6_YG.txt", row.names = F, quote = F)
write.table(cv6_GxE, "scv6_GxE.txt", row.names = F, quote = F)
write.table(cv6_VH, "scv6_VH.txt", row.names = F, quote = F)
write.table(cv6_YGVH, "scv6_YGVH.txt", row.names = F, quote = F)}



