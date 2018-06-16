# Create ped file
true$PAT = 0
true$MAT = 0
ped = true[,c(1,2,116,117,22,20,52:105)]
ped$Sex[ped$Sex == 1] = 2
ped$Sex[ped$Sex == 0] = 1
#ped[is.na(ped)] = -9

ped = ped[,rep(1:ncol(ped),each=2)]
ped = ped[,-c(2,4,6,8,10,12)]

ped2 = ped[,7:114]
ped2[ped2 == 1] = 2
ped2[ped2 == 0] = 1
ped = ped[,1:6]
ped = cbind(ped,ped2)

write.table(ped, "cnv.ped", quote = F, row.names = F, col.names = F)





a = as.data.frame(names(z[,7:60]))
map = cbind(map,a)
write.table(map[,c(1,4,2,3)], "cnv.map", row.names = F, col.names = F, quote = F)
