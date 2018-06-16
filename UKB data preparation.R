# Full UKB data manipulation


# Create variable that concatanate submitted and inffered sex
ukb$conc_sex <- paste(ukb$Submitted_Gender, ukb$Inferred_Gender, sep = "")                           # No of individuals 488377

# How many individuals with mismatched sex?
table(ukb$conc_sex)
'FF        FM     MF     MM 
264629    143    235    223370'   # 378


# Combine with family file to obtain FID and IID. Rows match
ukb <- cbind(fam, ukb)

# Remove heterozygosity and missing rates outliers
ukb2 <- ukb[which(ukb$het_missing_outliers == 0),]                                                    # No of outliers 968


# Select only those individuals with matching sex
toMatch <- c("MM", "FF")
matches <- ukb2[grepl(paste(toMatch, collapse="|"), ukb2$conc_sex), ]                                 # No of remaining individuals 487036



# Recode MM to 0, FF to 1
matches$conc_sex[matches$conc_sex == "MM"] <- 0; matches$conc_sex[matches$conc_sex == "FF"] <- 1


# Select columns that might be of interest
matches <- matches[,c(1,2,9,30,32:46,76)]
colnames(matches)[1] <- "FID"; colnames(matches)[2] <- "IID"; colnames(matches)[20] <- "Sex"           
colnames(matches)[3] <- "Array"; colnames(matches)[4] <- "WBA"

# Select only those who are white British
covar <- matches[which(matches$WBA == 1),]                                                             # No of individuals remaining    408664


# Bind with pheno file to get phenotypes and education status
covar <- merge(covar, edu, by = 1)                                                                     # No of individuals remaining    408661


# How many individuals are white (data from preliminary pheno file)
table(covar$Ethnicity)

'White 
408661'  # 100% agreement

covar <- covar[,c(1:20,24,27,30)]

# Do some data cleaning
colnames(covar)[2] <- "IID"; colnames(covar)[20] <- "Sex"
covar$Array[covar$Array == "UKBB"] <- 0; covar$Array[covar$Array == "UKBL"] <- 1


# leave only those with phenotypes
covar <- na.omit(covar)                                                                                # No of individuals remaining    86457


# How many people have uni degree? What is environmental prevalence?
table(covar$UniEdu)
'0     1 
55975 30311

30311/(55975 + 30311) 
[1] 0.3512853'


'# Remove individuals with other opthalmic conditions
final = merge(covar, out, by = 1)                                                                       # No of individuals remaining   86286
table(final$outlierMSE)
 0 
 86286 

final = final[,c(1:23)]
colnames(final)[2] <- "IID"; colnames(final)[20] <- "Sex"; colnames(final)[21] <- "Age"; colnames(final)[22] <- "avMSE"; colnames(final)[23] <- "UniEdu"
'
# Outliers have been removed by other screening steps

write.table(final, file = "covar.txt", quote = F, row.names = F)




# Create a list of individuals for creating bed,bim,fam files
id <- final[,1:2]
write.table(id, file = "ids_for_bed.txt", quote = F, row.names = F, col.names = F)


# Create a list of individuals for GRM constuction
GRM <- ukb[which(ukb$in_kinship_table == 0),] 
GRM <- GRM[,1:2]
colnames(GRM)[1] <- "FID"; colnames(GRM)[2] <- "IID"
forts <- merge(GRM, final, by = 1)                      # No of unrelated individuals 61448

GRM <- ukb[which(ukb$in_kinship_table == 1),] 
GRM <- GRM[,1:2]
colnames(GRM)[1] <- "FID"; colnames(GRM)[2] <- "IID"
fort <- merge(GRM, covar, by = 1)                       # No of individuals for relatedness testing 24838
fort <- fort[,1:2]
colnames(fort)[1] <- "FID"; colnames(fort)[2] <- "IID"

write.table(fort, file = "individuals_for_GRM.txt", quote = F, row.names = F)



# Total number of unrelated individuals
forts <- forts[,1:2]
colnames(forts)[1] <- "FID"; colnames(forts)[2] <- "IID"
total <- rbind(forts, unrelated)
write.table(total, file = "unrelated.txt", quote = F, row.names = F)

