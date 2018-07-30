# Quantitative association

library(data.table)
chr = list.files(pattern = ".linear", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
d = d[which(d$TEST == "ADD"),]

write.table(d, file ="destination.txt", quote = F, row.names = F)
