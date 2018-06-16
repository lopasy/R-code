#### Effect size map ####

stats2 = stats1[order(stats1$chr),]
stats2$seq = 0
#stats3 = stats2[,c(1,9,11)]


# Extract SNPs for each CHR
chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
        "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
for (i in chr){
  assign(paste0(i), stats[which(stats$chr == i),])
}

{chr1$CHR = 1; chr2$CHR = 2; chr3$CHR = 3; chr4$CHR = 4; chr5$CHR = 5; chr6$CHR = 6; chr7$CHR = 7
chr8$CHR = 8; chr9$CHR = 9; chr10$CHR = 10; chr11$CHR = 11; chr12$CHR = 12; chr13$CHR = 13; chr14$CHR = 14
chr15$CHR = 15; chr16$CHR = 16; chr17$CHR = 17; chr18$CHR = 18; chr19$CHR = 19; chr20$CHR = 20; chr21$CHR = 21
chr22$CHR = 22}



# Make a contigious alignment of the genome
{chr1$seq=1:nrow(chr1)
chr2$seq=1:nrow(chr2); chr2$seq = chr2$seq + (range(chr1$seq)[2]+100000)
chr3$seq=1:nrow(chr3); chr3$seq = chr3$seq + (range(chr2$seq)[2]+100000)
chr4$seq=1:nrow(chr4); chr4$seq = chr4$seq + (range(chr3$seq)[2]+100000)
chr5$seq=1:nrow(chr5); chr5$seq = chr5$seq + (range(chr4$seq)[2]+100000)
chr6$seq=1:nrow(chr6); chr6$seq = chr6$seq + (range(chr5$seq)[2]+100000)
chr7$seq=1:nrow(chr7); chr7$seq = chr7$seq + (range(chr6$seq)[2]+100000)
chr8$seq=1:nrow(chr8); chr8$seq = chr8$seq + (range(chr7$seq)[2]+100000)
chr9$seq=1:nrow(chr9); chr9$seq = chr9$seq + (range(chr8$seq)[2]+100000)
chr10$seq=1:nrow(chr10); chr10$seq = chr10$seq + (range(chr9$seq)[2]+100000)
chr11$seq=1:nrow(chr11); chr11$seq = chr11$seq + (range(chr10$seq)[2]+100000)
chr12$seq=1:nrow(chr12); chr12$seq = chr12$seq + (range(chr11$seq)[2]+100000)
chr13$seq=1:nrow(chr13); chr13$seq = chr13$seq + (range(chr12$seq)[2]+100000)
chr14$seq=1:nrow(chr14); chr14$seq = chr14$seq + (range(chr13$seq)[2]+100000)
chr15$seq=1:nrow(chr15); chr15$seq = chr15$seq + (range(chr14$seq)[2]+100000)
chr16$seq=1:nrow(chr16); chr16$seq = chr16$seq + (range(chr15$seq)[2]+100000)
chr17$seq=1:nrow(chr17); chr17$seq = chr17$seq + (range(chr16$seq)[2]+100000)
chr18$seq=1:nrow(chr18); chr18$seq = chr18$seq + (range(chr17$seq)[2]+100000)
chr19$seq=1:nrow(chr19); chr19$seq = chr19$seq + (range(chr18$seq)[2]+100000)
chr20$seq=1:nrow(chr20); chr20$seq = chr20$seq + (range(chr19$seq)[2]+100000)
chr21$seq=1:nrow(chr21); chr21$seq = chr21$seq + (range(chr20$seq)[2]+100000)
chr22$seq=1:nrow(chr22); chr22$seq = chr22$seq + (range(chr21$seq)[2]+100000)}




# Combine
chr = list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,
        chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
stats1 = do.call(rbind,chr)




plot(1, type="n", xlab = "Alignment in the genome", ylab = "Effect size", xlim = c(0, 9527300), ylim = c(-4, 3))
marks <- c(0,2000000,4000000,6000000,8000000,10000000)
axis(1,at=marks,labels=marks)
points(stats1$seq,stats1$effalt, type = "l")






# There are duplicates. Remove
chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
        "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
for (i in chr){
  assign(paste0(i), stats2[which(stats2$chr == i),])
}


{chr1 = chr1[!duplicated(chr1$pos),]
chr2 = chr2[!duplicated(chr2$pos), ]
chr3 = chr3[!duplicated(chr3$pos), ]
chr4 = chr4[!duplicated(chr4$pos), ]
chr5 = chr5[!duplicated(chr5$pos), ]
chr6 = chr6[!duplicated(chr6$pos), ]
chr7 = chr7[!duplicated(chr7$pos), ]
chr8 = chr8[!duplicated(chr8$pos), ]
chr9 = chr9[!duplicated(chr9$pos), ]
chr10 = chr10[!duplicated(chr10$pos), ]
chr11 = chr11[!duplicated(chr11$pos), ]
chr12 = chr12[!duplicated(chr12$pos), ]
chr13 = chr13[!duplicated(chr13$pos), ]
chr14 = chr14[!duplicated(chr14$pos), ]
chr15 = chr15[!duplicated(chr15$pos), ]
chr16 = chr16[!duplicated(chr16$pos), ]
chr17 = chr17[!duplicated(chr17$pos), ]
chr18 = chr18[!duplicated(chr18$pos), ]
chr19 = chr19[!duplicated(chr19$pos), ]
chr20 = chr20[!duplicated(chr20$pos), ]
chr21 = chr21[!duplicated(chr21$pos), ]
chr22 = chr22[!duplicated(chr22$pos), ]}

chr = list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,
           chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
stats1 = do.call(rbind,chr)



