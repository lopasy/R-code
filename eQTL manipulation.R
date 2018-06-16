# eQTL manipulation

true_gtex = merge(true, true_GxE, "SNPID")
predicted_gtex = merge(predicted, predicted_GxE, "SNPID")
rm(true_GxE, predicted_GxE)


library(data.table)
chr = list.files(pattern = ".pairs.txt", full.names = TRUE)
chr = lapply(chr, fread)
pairs = do.call(rbind.data.frame, chr)
rm(chr)
corrected = merge(pairs, gtex, "variant_id")
#colnames(egenes)[19]="SNP"
# GTEx
#colnames(gtex)[2]="SNP"
stopwords = c("_b37")
pairs$variant_id = gsub(paste0(stopwords,collapse = "|"), "", pairs$variant_id)
pairs = as.data.frame(pairs)
# The above file contains all sig gene-variant pairs from GTEx with variants named as chr_pos_ref_alt

zed = merge(predicted, gtex, "SNP") # This is for adding ref allele to predicted output
zed = zed[,-c(2,6,8,10,12,14,15,16,17,19)] # This is to make computation faster
zed$variant_id = paste(zed$CHR, zed$BP, zed$V6, zed$A1, sep = "_") # Reference first
zed2 = zed
zed2$variant_id = paste(zed$CHR, zed$BP, zed$A1, zed$V6, sep = "_") # Reference second
zed_both = rbind(zed,zed2)
final = merge(pairs, zed_both, "variant_id") # This gets the number of eQTLs that are found in UKB

{pruned = read.table("~/Documents/UKB/gxescan/gtex/gtex_pruned.prune.in", quote="\"", comment.char="")
colnames(pruned)="SNP"
pruned = merge(test, pruned, "SNP")
length(unique(pruned$SNP))}


final = merge(pairs, gtex, "variant_id")
predicted=merge(predicted,YG, "SNPID")
predicted=predicted[,c(3,9)]
colnames(predicted)[1]="rs_id_dbSNP147_GRCh37p13"
predicted2=merge(final,predicted, "rs_id_dbSNP147_GRCh37p13")
length(unique(predicted2$gene_id))








































gtex$variant_id = paste(gtex$V1, gtex$V4, gtex$V6, gtex$V5, sep = "_")
colnames(gtex)[7] = "variant_id"
final2 = merge(pairs, gtex, "variant_id")
final_both = rbind(final, final2)
x = final_both[which(final_both$maf >0.049999999999999999999999),]
#d2 = d[unique(d$SNP),]
# To obtain the list of eGenes, select the rows with 'qval' < 0.05
qval = d[which(d$qval < 0.05),]

# How many SNPs without rsID
length(which(d$rs_id_dbSNP147_GRCh37p13 == "."))
length(which(qval$rs_id_dbSNP147_GRCh37p13 == "."))



qval$name = paste(qval$chr, qval$pos, qval$ref, qval$alt, sep = "_")
test = predicted.bim[,c(2,6)] # get reference allele from bim
colnames(test)=c("SNP","REF")
qtl = merge(predicted,test,"SNP")
qtl$name = paste(qtl$CHR, qtl$BP, qtl$REF, qtl$A1, sep = "_")

# Create a subset of SNPs for GxE testing
qvals = as.data.frame(qval$name);colnames(qvals) = "name"
subset = merge(qtl,qvals, "name")

# Match GTEx SNPs to those from UKB
colnames(qval)[19] = "SNP"
gtex_ukb = merge(predicted, qval, "SNP")
unique = gtex_ukb[!duplicated(gtex_ukb$SNP),]
final = unique[which(unique$A1 == unique$alt),]

subsets = subset[!duplicated(subset$SNP),]

comparison = merge(gtex_ukb, subset, "SNP")

# Compare alleles
alts = as.data.frame(unique$A1)
alts$gtex = unique$alt
length(which(alts[,1] == alts[,2]))


# Additional
qvald = qval[!duplicated(qval$SNP),]
length(which(!duplicated(qval$SNP)))
########################################################################
########################################################################
########################################################################

# True phenotype
test = true.bim[,c(2,6)] # get reference allele from bim
colnames(test)=c("SNP","REF")
qtl = merge(true,test,"SNP")
qtl$name = paste(qtl$CHR, qtl$BP, qtl$REF, qtl$A1, sep = "_")

# Create a subset of SNPs for GxE testing
qvals = as.data.frame(qval$name);colnames(qvals) = "name"
subset = merge(qtl,qvals, "name")

# Match GTEx SNPs to those from UKB
colnames(qval)[19] = "SNP"
gtex_ukb = merge(true, qval, "SNP")
unique = gtex_ukb[!duplicated(gtex_ukb$SNP),]
subsets = subset[!duplicated(subset$SNP),]

alts = as.data.frame(unique$A1)
alts$gtex = unique$alt
length(which(alts[,1] == alts[,2]))
final = unique[which(unique$A1 == unique$alt),]

