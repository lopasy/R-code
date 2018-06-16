# Data pre-processing
"For this tutorial we use genotype data files formatted for use with PLINK software. We utilize the function, read.plink from snpStats, which allows the reading in of data formatted as .bed, .bim, and .fam files. The .bed file contains the genotype information, coded in binary. The .bim file contains information for each SNP with a respective column for each of the following information: chromosome number, SNP name (typically an rs #), genetic distance (not necessary for this tutorial), chromosomal position, identity of allele 1, and identity of allele 2. The assignment of allele 1 and allele 2, is related to the effect allele, or the allele that is being counted when we assign a numeric value to a genotype. This is typically assigned based on allele frequency, though not always. In this tutorial, allele 1 pertains to the minor, or less common allele. Lastly, the .fam file contains information for each samples with a respective column for each of the following information: family ID (this will be used to identify each sample when read into R), individual ID, paternal ID, maternal ID, sex (coded as 1 = male, 2 = female), and phenotype. In this tutorial we utilize a supplemental clinical file for outcome variables and additional covariates.

Alternatively, similar genotype information can also be formatted for PLINK software as .ped and .map files. The information of the .ped file can be thought of as a combination of the .bed and .fam files. It is a large table with the first six columns identical to a .fam file, followed by a columns containing the genotype data for each SNP. The .map file contains the first four columns of the .bim file, without the allele assignments. These files can be read in using the function, read.pedfile, from snpStats. More information about the formatting of these files can be found on the PLINK website.
"
install.packages("snpStats")
library(snpStats)

# Read in PLINK files
geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))
# The geno object contains a genotype member of type SnpMatrix where each column is a SNP and each row is a sample. For convenience, we assign that to the object, genotype. Within this object individual genotypes are assigned in the SnpMatrix specific RAW format. As a result, many functions dependent on the snpStats package specifically require SnpMatrix objects, while more basic functions will require conversion to other formats. geno also contains the SNP information from the .bim file within the map member that we assign to genoBim.

# Obtain the SnpMatrix object (genotypes) table from geno list
# Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
print(genotype)                  # 861473 SNPs read in for 1401 subjects
## A SnpMatrix with  1401 rows and  861473 columns
## Row names:  10002 ... 11596 
## Col names:  rs10458597 ... rs5970564
#Obtain the SNP information from geno list
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))
##            chr        SNP gen.dist position   A1 A2
## rs10458597   1 rs10458597        0   564621 <NA>  C
## rs12565286   1 rs12565286        0   721290    G  C
## rs12082473   1 rs12082473        0   740857    T  C
## rs3094315    1  rs3094315        0   752566    C  T
## rs2286139    1  rs2286139        0   761732    C  T
## rs11240776   1 rs11240776        0   765269    G  A
# Remove raw file to free up memory
rm(geno)
# Supplemental clinical data is found in a corresponding CSV file for each sample. It contains a column for the sample ID (Family ID in the .fam file) and a respective column for each of the following variables: coronary artery disease status (coded as 0 = control and 1 = affected), sex (coded as 1 = male, 2 = female), age (years), triglyceride level (mg/dL), high-density lipoprotein level (mg/dL), low-density lipoprotein level (mg/dL).

# Read in clinical file
clinical <- read.csv(clinical.fn,
colClasses=c("character", "factor", "factor", rep("numeric", 4)))
rownames(clinical) <- clinical$FamID
print(head(clinical))
##       FamID CAD sex age  tg hdl ldl
## 10002 10002   1   1  60  NA  NA  NA
## 10004 10004   1   2  50  55  23  75
## 10005 10005   1   1  55 105  37  69
## 10007 10007   1   1  52 314  54 108
## 10008 10008   1   1  58 161  40  94
## 10009 10009   1   1  59 171  46  92
# We filter the genotype data to only include samples with corresponding clinical data by indexing the genotype object using only row names that match the sample IDs.

# Subset genotype for subject data
genotype <- genotype[clinical$FamID, ]
print(genotype)  # Tutorial: All 1401 subjects contain both clinical and genotype data
## A SnpMatrix with  1401 rows and  861473 columns
## Row names:  10002 ... 11596 
## Col names:  rs10458597 ... rs5970564
# Write genotype, genoBim, clinical for future use
save(genotype, genoBim, clinical, file = working.data.fname(1))

# SNP level filtering

# Once the data is loaded, we are ready to remove SNPs that fail to meet minimum criteria due to missing data, low variability or genotyping errors. snpStats provides functions, col.summary and row.summary, that return statistics on SNPs and samples, respectively.

# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print(head(snpsum.col))
##            Calls Call.rate Certain.calls       RAF         MAF       P.AA
## rs10458597  1398 0.9978587             1 1.0000000 0.000000000 0.00000000
## rs12565286  1384 0.9878658             1 0.9483382 0.051661850 0.00433526
## rs12082473  1369 0.9771592             1 0.9985391 0.001460920 0.00000000
## rs3094315   1386 0.9892934             1 0.8217893 0.178210678 0.04761905
## rs2286139   1364 0.9735903             1 0.8621701 0.137829912 0.02199413
## rs11240776  1269 0.9057816             1 0.9988180 0.001182033 0.00000000
##                   P.AB      P.BB       z.HWE
## rs10458597 0.000000000 1.0000000          NA
## rs12565286 0.094653179 0.9010116 -1.26529432
## rs12082473 0.002921841 0.9970782  0.05413314
## rs3094315  0.261183261 0.6911977 -4.03172248
## rs2286139  0.231671554 0.7463343 -0.93146122
## rs11240776 0.002364066 0.9976359  0.04215743
# Using these summary statistics, we keep the subset of SNPs that meet our criteria for minimum call rate and minor allele frequency.

# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") #203287 SNPs will be removed
## 203287 SNPs will be removed due to low MAF or call rate.
# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print(genotype)                           # 658186 SNPs remain
## A SnpMatrix with  1401 rows and  658186 columns
## Row names:  10002 ... 11596 
## Col names:  rs12565286 ... rs5970564
# Write subsetted genotype data and derived results for future use
save(genotype, snpsum.col, genoBim, clinical, file=working.data.fname(2))

# Sample level filtering

# Sample level filtering

source("globals.R")

# load data created in previous snippets
load(working.data.fname(2))

library(snpStats)
"The second stage of data pre-processing involves filtering samples, i.e. removing individuals due to missing data, sample contamination, correlation (for population-based investigations), and racial/ethnic or gender ambiguity or discordance. In our study, we address these issues by filtering on call rate, heterozygosity, cryptic relatedness and duplicates using identity-by-descent, and we visually assess ancestry.

Basic sample filtering

Sample level quality control for missing data and heterozygosity is achieved using the row.summary function from snpStats. An additional heterozygosity F statistic calculation is carried out with the form, |F|=(1???O/E)|F|=(1???O/E), where OO is observed proportion of heterozygous genotypes for a given sample and EE is the expected proportion of heterozygous genotypes for a given sample based on the minor allele frequency across all non-missing SNPs for a given sample.
"

library(SNPRelate)                      # LD pruning, relatedness, PCA
library(plyr)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

head(snpsum.row)
##       Call.rate Certain.calls Heterozygosity          hetF
## 10002 0.9826554             1      0.3289825 -0.0247708291
## 10004 0.9891581             1      0.3242931 -0.0103236529
## 10005 0.9918427             1      0.3231825 -0.0062550972
## 10007 0.9861027             1      0.3241469 -0.0098475016
## 10008 0.9823333             1      0.3228218 -0.0075941985
## 10009 0.9913034             1      0.3213658 -0.0002633189
# We apply filtering on call rate and heterozygosity, selecting only those samples that meet our criteria.

# Setting thresholds
sampcall <- 0.95    # Sample call rate cut-off
hetcutoff <- 0.1    # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE    # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n") #0 subjects removed
## 0 subjects will be removed due to low sample call rate or inbreeding coefficient.
# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[ rownames(genotype), ]

# IBD analysis

# In addition to these summary statistics, we also want to filter on relatedness criteria. We use the SNPRelate package to perform identity-by-descent (IBD) analysis. This package requires that the data be transformed into a GDS format file. IBD analysis is performed on only a subset of SNPs that are in linkage equilibrium by iteratively removing adjacent SNPs that exceed an LD threshold in a sliding window using the snpgdsLDpruning function.

# Checking for Relatedness

ld.thresh <- 0.2    # LD cut-off
kin.thresh <- 0.1   # Kinship cut-off

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)
## Start snpgdsBED2GDS ...
##  BED file: "/Users/ericreed/Desktop/FoulkesLab/SIMFiles/GWAStutorial.bed" in the SNP-major mode (Sample X SNP)
##  FAM file: "/Users/ericreed/Desktop/FoulkesLab/SIMFiles/GWAStutorial.fam", DONE.
##  BIM file: "/Users/ericreed/Desktop/FoulkesLab/SIMFiles/GWAStutorial.bim", DONE.
## Wed Jun 24 16:51:03 2015     store sample id, snp id, position, and chromosome.
##  start writing: 1401 samples, 861473 SNPs ...
##      Wed Jun 24 16:51:03 2015    0%
##      Wed Jun 24 16:51:15 2015    100%
## Wed Jun 24 16:51:15 2015     Done.
## Optimize the access efficiency ...
## Clean up the fragments of GDS file:
##  open the file "/Users/ericreed/Desktop/FoulkesLab/SIMFiles/GWAStutorial.gds" (size: 308167905).
##  # of fragments in total: 39.
##  save it to "/Users/ericreed/Desktop/FoulkesLab/SIMFiles/GWAStutorial.gds.tmp".
##  rename "/Users/ericreed/Desktop/FoulkesLab/SIMFiles/GWAStutorial.gds.tmp" (size: 308167653).
##  # of fragments in total: 18.
genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)

# Automatically added "-1" sample suffixes are removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

#Prune SNPs for IBD analysis
set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
sample.id = geno.sample.ids, # Only analyze the filtered samples
snp.id = colnames(genotype)) # Only analyze the filtered SNPs
## Hint: it is suggested to call `snpgdsOpen' to open a SNP GDS file instead of `openfn.gds'.
## SNP pruning based on LD:
## Excluding 203287 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 658186 SNPs
##  Using 1 (CPU) core
##  Sliding window: 500000 basepairs, Inf SNPs
##  |LD| threshold: 0.2
## Chromosome 1: 8.25%, 5863/71038
## Chromosome 3: 8.10%, 4906/60565
## Chromosome 6: 8.06%, 4364/54176
## Chromosome 12: 8.59%, 3619/42124
## Chromosome 21: 9.40%, 1171/12463
## Chromosome 2: 7.67%, 5655/73717
## Chromosome 4: 8.23%, 4582/55675
## Chromosome 7: 8.51%, 3947/46391
## Chromosome 11: 7.90%, 3495/44213
## Chromosome 10: 8.01%, 3837/47930
## Chromosome 8: 7.68%, 3709/48299
## Chromosome 5: 8.08%, 4537/56178
## Chromosome 14: 8.79%, 2467/28054
## Chromosome 9: 8.25%, 3392/41110
## Chromosome 17: 11.17%, 2227/19939
## Chromosome 13: 8.36%, 2863/34262
## Chromosome 20: 9.40%, 2139/22753
## Chromosome 15: 9.25%, 2396/25900
## Chromosome 16: 9.30%, 2566/27591
## Chromosome 18: 8.90%, 2335/26231
## Chromosome 19: 13.01%, 1494/11482
## Chromosome 22: 10.96%, 1248/11382
## 72812 SNPs are selected in total.
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")  # Tutorial: expect 72812 SNPs
## 72812 will be used in IBD analysis
# The snpgdsIBDMoM function computes the IBD coefficients using method of moments. The result is a table indicating kinship among pairs of samples.

# Find IBD coefficients using Method of Moments procedure.  Include pairwise kinship.
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
sample.id = geno.sample.ids,
snp.id = snpset.ibd,
num.thread = 1)
## Hint: it is suggested to call `snpgdsOpen' to open a SNP GDS file instead of `openfn.gds'.
## IBD analysis (PLINK method of moment) on SNP genotypes:
## Excluding 788661 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 72812 SNPs
##  Using 1 (CPU) core
## PLINK IBD:   the sum of all working genotypes (0, 1 and 2) = 32757268
## PLINK IBD:   Wed Jun 24 16:51:53 2015    0%
## PLINK IBD:   Wed Jun 24 16:52:21 2015    100%
ibdcoeff <- snpgdsIBDSelection(ibd)     # Pairwise sample comparison
head(ibdcoeff)
##     ID1   ID2        k0         k1    kinship
## 1 10002 10004 0.9201072 0.07989281 0.01997320
## 2 10002 10005 0.9478000 0.05220002 0.01305001
## 3 10002 10007 0.9209875 0.07901253 0.01975313
## 4 10002 10008 0.9312527 0.06874726 0.01718682
## 5 10002 10009 0.9386937 0.06130626 0.01532656
## 6 10002 10010 0.9146065 0.08539354 0.02134839
# Using the IBD pairwise sample relatedness measure, we iteratively remove samples that are too similar using a greedy strategy in which the sample with the largest number of related samples is removed. The process is repeated until there are no more pairs of samples with kinship coefficients above our cut-off.

# Check if there are any candidates for relatedness
ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]

# iteratively remove samples with high kinship starting with the sample with the most pairings
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {

# count the number of occurrences of each and take the top one
sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
rm.sample <- sample.counts[1, 'x']
cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'],'other samples.\n')

# remove from ibdcoeff and add to list
ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n") 
## 0 similar samples removed due to correlation coefficient >= 0.1
print(genotype)                         # Tutorial: expect all 1401 subjects remain
## A SnpMatrix with  1401 rows and  658186 columns
## Row names:  10002 ... 11596 
## Col names:  rs12565286 ... rs5970564

# Ancestry

# To better understand ancestry, we plot the first two principal components of the genotype data. Principal component calculation is achieved via the snpgdsPCA function from SNPRelate. It is important to note that in this example we are reasonably confident that our samples are homogeneous, coming from european ancestry. Therefore, given that there are no clear outliers, we fail to remove any samples.

# Checking for ancestry

# Find PCA matrix
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=1)
## Hint: it is suggested to call `snpgdsOpen' to open a SNP GDS file instead of `openfn.gds'.
## Principal Component Analysis (PCA) on SNP genotypes:
## Excluding 788661 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 72812 SNPs
##  Using 1 (CPU) core
## PCA: the sum of all working genotypes (0, 1 and 2) = 32757268
## PCA: Wed Jun 24 16:54:27 2015    0%
## PCA: Wed Jun 24 16:54:57 2015    100%
## PCA: Wed Jun 24 16:54:57 2015    Begin (eigenvalues and eigenvectors)
## PCA: Wed Jun 24 16:54:58 2015    End (eigenvalues and eigenvectors)
# Create data frame of first two principal comonents
pctab <- data.frame(sample.id = pca$sample.id,
PC1 = pca$eigenvect[,1],    # the first eigenvector
PC2 = pca$eigenvect[,2],    # the second eigenvector
stringsAsFactors = FALSE)

# Plot the first two principal comonents
plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", main = "Ancestry Plot")
# plot of chunk code3-f

# Close GDS file
closefn.gds(genofile)

# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, file=working.data.fname(3))

# SNP Filtering - HWE filtering on control samples

# Finally, once samples are filtered, we return to SNP level filtering and apply a check of Hardy-Weinberg equilibrium. Rejection of Hardy-Weinberg equilibrium can be an indication of population substructure or genotyping errors. Given that we are performing a statistical test at every SNP, it is common to use a relatively lenient cut-off. In this example we only remove SNPs with p-values, corresponding to the HWE test statistic on CAD controls, of less than 1×10???61×10???6. We only test HWE on CAD controls due to possible violation of HWE caused by disease association.

# Hardy-Weinberg SNP filtering on CAD controls

hardy <- 10^-6      # HWE cut-off

CADcontrols <- clinical[ clinical$CAD==0, 'FamID' ]
snpsum.colCont <- col.summary( genotype[CADcontrols,] )
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)

HWEuse[is.na(HWEuse)] <- FALSE          # Remove NA's as well
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")  # 1296 SNPs removed
## 1296 SNPs will be removed due to high HWE.
# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]

print(genotype)                           # 656890 SNPs remain
## A SnpMatrix with  1401 rows and  656890 columns
## Row names:  10002 ... 11596 
## Col names:  rs12565286 ... rs28729663
# Overwrite old genotype with new filtered version
save(genotype, genoBim, clinical, file=working.data.fname(4))





"Genome-wide association analysis
Now that our data is loaded, filtered, and additional SNP genotypes imputed we are ready to perform genome-wide association analysis. This involves regressing each SNP separately on a given trait, adjusted for sample level clinical, environmental, and demographic factors. Due to the large number of SNPs and the generally uncharacterized relationships to the outcome, a simple single additive model will be employed.

The GWAA function requires two arguments. The genodata argument should specify the entire genotype data object in SnpMatrix format. The phenodata argument should be a data frame with a column of sample IDs, corresponding to the row names of genodata, and a columns for the continuous outcome variable. These columns must be named "id" and "phenotype", respectively. In order to fit the model, genotype data is converted to numeric format using the as function from snpStats. The genotypes of each SNP are then coded as continuous, thereby taking on the value of 0, 1, and 2. For this example, we wish for the value of the genotype to reflect the number of minor alleles. However, following conversion our values will reflect the opposite. To fix this a flip.matrix procedure is included in our GWAA function, which can be turned on or off using the flip argument.

Due to the large number of models that require fitting, the GWA analysis can be deployed in parallel across multiple processors or machines to reduce the running time. Here we demonstrate two basic methods for performing parallel processing using the doParallel package. This will be carried out differently depending on whether or not the analysis is run on a UNIX based system, though the arguments are the same. The user can specify the number of processes using the worker argument (set to 2 by default). Additional arguments include select.snps and nSplits. The former allows the user to subset the analysis via a vector of SNP IDs. The latter specifies a number of SNP-wise splits that are made to the genotype data. The function runs the GWA analysis on these smaller subsets of the genotype data one at a time. After each subset has finished running the function will print a progress update onto the R console. By default this is set to 10.

GWAA function"

# Genome-wide Association Analysis
# Parallel implementation of linear model fitting on each SNP

GWAA <- function(genodata=genotypes,  phenodata=phenotypes, family = gaussian, filename=NULL,
                 append=FALSE, workers=getOption("mc.cores",2L), flip=TRUE,
                 select.snps=NULL, hosts=NULL, nSplits=10)
{
  if (!require(doParallel)) { stop("Missing doParallel package") }
  
  #Check that a filename was specified
  if(is.null(filename)) stop("Must specify a filename for output.")
  
  #Check that the genotype data is of class 'SnpMatrix'
  if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")
  
  #Check that there is a variable named 'phenotype' in phenodata table
  if( !"phenotype" %in% colnames(phenodata))  stop("Phenotype data must have column named 'phenotype'")
  
  #Check that there is a variable named 'id' in phenodata table
  if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")
  
  #If a vector of SNPs is given, subset genotype data for these SNPs
  if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata)%in%select.snps)]
  
  #Check that there are still SNPs in 'SnpMatrix' object
  if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")
  
  #Print the number of SNPs to be checked
  cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
  
  #If append=FALSE than we will overwrite file with column names
  if(!isTRUE(append)) {
    columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
    write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  
  # Check sample counts
  if (nrow(phenodata) != nrow(genodata)) {
    warning("Number of samples mismatch.  Using subset found in phenodata.")
  }
  
  # Order genodata rows to be the same as phenodata
  genodata <- genodata[phenodata$id,]
  
  cat(nrow(genodata), "samples included in analysis.\n")
  
  # Change which allele is counted (major or minor)
  flip.matrix<-function(x) {
    zero2 <- which(x==0)
    two0 <- which(x==2)
    x[zero2] <- 2
    x[two0] <- 0
    return(x)
  }
  
  nSNPs <- ncol(genodata)
  genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
  
  snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
  snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group
  
  if (is.null(hosts)) {
    # On Unix this will use fork and mclapply.  On Windows it
    # will create multiple processes on localhost.
    cl <- makeCluster(workers)
  } else {
    # The listed hosts must be accessible by the current user using
    # password-less ssh with R installed on all hosts, all 
    # packages installed, and "rscript" is in the default PATH.
    # See docs for makeCluster() for more information.
    cl <- makeCluster(hosts, "PSOCK")
  }
  show(cl)                            # report number of workers and type of parallel implementation
  registerDoParallel(cl)
  
  foreach (part=1:nSplits) %do% {
    # Returns a standar matrix of the alleles encoded as 0, 1 or 2
    genoNum <- as(genodata[,snp.start[part]:snp.stop[part]], "numeric")
    
    # Flip the numeric values of genotypes to count minor allele
    if (isTRUE(flip)) genoNum <- flip.matrix(genoNum)
    
    # For each SNP, concatenate the genotype column to the
    # phenodata and fit a generalized linear model
    rsVec <- colnames(genoNum)
    res <- foreach(snp.name=rsVec, .combine='rbind') %dopar% {
      a <- summary(glm(phenotype~ . - id, family=family, data=cbind(phenodata, snp=genoNum[,snp.name])))
      a$coefficients['snp',]
    }
    
    # write results so far to a file
    write.table(cbind(rsVec,res), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 100*part/nSplits))
  }
  
  stopCluster(cl)
  
  return(print("Done."))
}

"Phenotype data preparation

First we create a data frame of phenotype features that is the concatenation of clinical features and the first ten principal components. The HDL feature is normalized using a rank-based inverse normal transform. We then remove variables that we are not including in the analysis, i.e. HDL(non-normalized), LDL, TG, and CAD. Finally, we remove samples with missing normalized HDL data.
"

library(GenABEL)
source("GWAA.R")

# Merge clincal data and principal components to create phenotype table
phenoSub <- merge(clinical,pcs)      # data.frame => [ FamID CAD sex age hdl pc1 pc2 ... pc10 ]

# We will do a rank-based inverse normal transformation of hdl
phenoSub$phenotype <- rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$phenotype, main="Histogram of Tranformed HDL", xlab="Transformed HDL")
# plot of chunk code7-a

# Remove unnecessary columns from table
phenoSub$hdl <- NULL
phenoSub$ldl <- NULL
phenoSub$tg <- NULL
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace=c(FamID="id"))

# Include only subjects with hdl data
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
# 1309 subjects included with phenotype data

print(head(phenoSub))
##      id sex age          pc1          pc2          pc3           pc4
## 2 10004   2  50 -0.012045108 -0.007231015 -0.003001290 -0.0107972693
## 3 10005   1  55 -0.016702930 -0.005347697  0.014449836 -0.0006151058
## 4 10007   1  52 -0.009537235  0.004556977  0.002683566  0.0166255657
## 5 10008   1  58 -0.015392106 -0.002446933  0.020508791 -0.0057241772
## 6 10009   1  59 -0.015123858 -0.002353917  0.021360452  0.0069156529
## 7 10010   1  54 -0.012816157  0.005126124  0.014654847 -0.0147082270
##             pc5           pc6           pc7          pc8          pc9
## 2 -0.0077705400 -0.0046457510  0.0018061075 -0.003087891 -0.001833242
## 3  0.0345170160  0.0387085513  0.0205790788 -0.012265508  0.003592690
## 4 -0.0002363142  0.0055146271  0.0159588869  0.027975455  0.029777180
## 5 -0.0039696226  0.0053542437 -0.0007269312  0.027014714  0.010672162
## 6  0.0400677558  0.0232224781  0.0152485234  0.013296852  0.022746352
## 7 -0.0008190769 -0.0003831342 -0.0131606658 -0.013647709 -0.008912913
##           pc10  phenotype
## 2 -0.004538746 -2.2877117
## 3 -0.002287043 -0.4749316
## 4 -0.007461255  0.8855512
## 5 -0.003352997 -0.1644639
## 6  0.013143889  0.3938940
## 7 -0.056187339  1.7109552

"Parallel model fitting

Using this phenotype data, we perform model fitting on each of the typed SNPs in thegenotype object and write the results to a .txt file."

# Run GWAS analysis
# Note: This function writes a file, but does not produce an R object
start <- Sys.time()
GWAA(genodata=genotype, phenodata=phenoSub, filename=gwaa.fname)
## Loading required package: doParallel
## Loading required package: foreach
## Loading required package: iterators
## Loading required package: parallel
## 656890  SNPs included in analysis.
## 1309 samples included in analysis.
## socket cluster with 2 nodes on host 'localhost'
## GWAS SNPs 1-65689 (10% finished)
## GWAS SNPs 65690-131378 (20% finished)
## GWAS SNPs 131379-197067 (30% finished)
## GWAS SNPs 197068-262756 (40% finished)
## GWAS SNPs 262757-328445 (50% finished)
## GWAS SNPs 328446-394134 (60% finished)
## GWAS SNPs 394135-459823 (70% finished)
## GWAS SNPs 459824-525512 (80% finished)
## GWAS SNPs 525513-591201 (90% finished)
## GWAS SNPs 591202-656890 (100% finished)
## [1] "Done."
end <- Sys.time()
print(end-start)
## Time difference of 2.259378 hours
# Add phenosub to saved results
save(genotype, genoBim, clinical, pcs, imputed, target, rules, phenoSub, support, file=working.data.fname(7))


"Model fitting of non-typed SNPs

We also perform association testing on additional SNPs from genotype imputation. Here we use thesnp.rhs.tests function from snpStats to perform the analysis based on the imputation "rules" we calculated previously. We need to specify the variables from the phenoSub data frame that we are including in the model with row names corresponding to the sample IDs.

The resulting SNPs are combined with the chromosome position information to create a table of SNPs, location and p-value. Finally, we take ???log10???log10 of the p-value for plotting.
"
# Carry out association testing for imputed SNPs using snp.rhs.tests()
rownames(phenoSub) <- phenoSub$id

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                     family = "Gaussian", data = phenoSub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs by calling methods on the returned GlmTests object.
results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results <- results[!is.na(results$p.value),]

#Write a file containing the results
write.csv(results, impute.out.fname, row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
imputeOut<-merge(results, support[, c("SNP", "position")])
imputeOut$chr <- 16

imputeOut$type <- "imputed"

# Find the -log_10 of the p-values
imputeOut$Neg_logP <- -log10(imputeOut$p.value)

# Order by p-value
imputeOut <- arrange(imputeOut, p.value)
print(head(imputeOut))
##          SNP      p.value position chr    type Neg_logP
## 1  rs1532624 9.805683e-08 57005479  16 imputed 7.008522
## 2  rs7205804 9.805683e-08 57004889  16 imputed 7.008522
## 3 rs12446515 1.430239e-07 56987015  16 imputed 6.844591
## 4 rs17231506 1.430239e-07 56994528  16 imputed 6.844591
## 5   rs173539 1.430239e-07 56988044  16 imputed 6.844591
## 6   rs183130 1.430239e-07 56991363  16 imputed 6.844591

"Mapping associated SNPs to genes

Using a separate data file containing the chromosome and coordinate locations of each protein coding gene, we can locate coincident genes and SNPs.

We use the following function to extract the subset of SNPs that are near a gene of interest."

# Returns the subset of SNPs that are within extend.boundary of gene
# using the coords table of gene locations.
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000) {
  coordsSub <- coords[coords$gene == gene,] #Subset coordinate file for spcified gene
  
  coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries
  coordsSub$stop <- coordsSub$stop + extend.boundary
  
  SNPsub <- SNPs[SNPs$position >= coordsSub$start & SNPs$position <= coordsSub$stop &
                   SNPs$chr == coordsSub$chr,] #Subset for SNPs in gene
  
  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}
"The SNP with the lowest p-value in both the typed and imputed SNP analysis lies within the boundaries of the cholesteryl ester transfer protein gene, CETP. We can call the map2gene function for "CETP" to filter the imputed genotypes and extract only those SNPs that are near CETP. This will be used for post-analytic interrogation to follow.
"
# Read in file containing protein coding genes coords
genes <- read.csv(protein.coding.coords.fname, stringsAsFactors = FALSE)

# Subset for CETP SNPs
impCETP <- map2gene("CETP", coords = genes, SNPs = imputeOut)

# Filter only the imputed CETP SNP genotypes 
impCETPgeno <- imputed[, impCETP$SNP ]
save(genotype, genoBim, clinical, pcs, imputed, target, rules,
     phenoSub, support, genes, impCETP, impCETPgeno, imputeOut, file = working.data.fname(8))






"Post-analytic visualization and genomic interrogation
We now have generated and fit both typed and imputed genotypes. The next step is to combine the results, and isolate just those SNPs in our region of interest. Following similar steps as for imputed SNPs, the typed SNPs are loaded from a file generated by the GWAA function. We follow similar steps to attach chromosome and position to each SNP, order by significance, and take ???log10???log10 of the p-value.
"
# Read in GWAS output that was produced by GWAA function
GWASout <- read.table(gwaa.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])
rm(genoBim)

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))
##          SNP   Estimate  Std.Error   t.value      p.value Neg_logP chr
## 1  rs1532625  0.2024060 0.03756207  5.388575 8.452365e-08 7.073022  16
## 2   rs247617  0.2119357 0.03985979  5.317030 1.243480e-07 6.905361  16
## 3 rs10945761  0.1856564 0.04093602  4.535282 6.285358e-06 5.201670   6
## 4  rs3803768 -0.3060086 0.06755628 -4.529685 6.451945e-06 5.190309  17
## 5  rs4821708 -0.1816673 0.04020915 -4.518058 6.825085e-06 5.165892  22
## 6  rs9647610  0.1830434 0.04072772  4.494320 7.607161e-06 5.118777   6
##    position
## 1  57005301
## 2  56990716
## 3 162065367
## 4  80872028
## 5  38164106
## 6 162066421

"Isolate CETP-specific SNPs

The two tables of typed and imputed genotypes are combined into a single table. In addition, we also concatenate just the SNPs near CETP and display them all here.
"
# Combine typed and imputed
GWASout$type <- "typed"

GWAScomb<-rbind.fill(GWASout, imputeOut)
head(GWAScomb)
##          SNP   Estimate  Std.Error   t.value      p.value Neg_logP chr
## 1  rs1532625  0.2024060 0.03756207  5.388575 8.452365e-08 7.073022  16
## 2   rs247617  0.2119357 0.03985979  5.317030 1.243480e-07 6.905361  16
## 3 rs10945761  0.1856564 0.04093602  4.535282 6.285358e-06 5.201670   6
## 4  rs3803768 -0.3060086 0.06755628 -4.529685 6.451945e-06 5.190309  17
## 5  rs4821708 -0.1816673 0.04020915 -4.518058 6.825085e-06 5.165892  22
## 6  rs9647610  0.1830434 0.04072772  4.494320 7.607161e-06 5.118777   6
##    position  type
## 1  57005301 typed
## 2  56990716 typed
## 3 162065367 typed
## 4  80872028 typed
## 5  38164106 typed
## 6 162066421 typed
tail(GWAScomb)
##               SNP Estimate Std.Error t.value   p.value     Neg_logP chr
## 818521 rs62048372       NA        NA      NA 0.9999838 7.048600e-06  16
## 818522  rs8056666       NA        NA      NA 0.9999838 7.048600e-06  16
## 818523  rs8057313       NA        NA      NA 0.9999838 7.048600e-06  16
## 818524  rs8061812       NA        NA      NA 0.9999838 7.048600e-06  16
## 818525  rs9940700       NA        NA      NA 0.9999838 7.048600e-06  16
## 818526 rs13334556       NA        NA      NA 0.9999843 6.825503e-06  16
##        position    type
## 818521 53775940 imputed
## 818522 53794830 imputed
## 818523 53794855 imputed
## 818524 53794856 imputed
## 818525 53795409 imputed
## 818526  5463800 imputed
# Subset for CETP SNPs
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP SNPs from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP)[,c("SNP","p.value","Neg_logP","chr","position","type","gene")]
print(CETP)
##            SNP      p.value   Neg_logP chr position    type gene
## 1    rs1532625 8.452365e-08 7.07302173  16 57005301   typed CETP
## 2     rs289742 3.788738e-04 3.42150548  16 57017762   typed CETP
## 3     rs289715 4.299934e-03 2.36653823  16 57008508   typed CETP
## 4    rs6499863 1.382602e-02 1.85930275  16 56992017   typed CETP
## 5    rs1800777 8.833782e-02 1.05385333  16 57017319   typed CETP
## 6    rs4783962 1.039467e-01 0.98318933  16 56995038   typed CETP
## 7   rs12708980 6.375740e-01 0.19546941  16 57012379   typed CETP
## 8    rs1532624 9.805683e-08 7.00852215  16 57005479 imputed CETP
## 9    rs7205804 9.805683e-08 7.00852215  16 57004889 imputed CETP
## 10  rs17231506 1.430239e-07 6.84459142  16 56994528 imputed CETP
## 11    rs183130 1.430239e-07 6.84459142  16 56991363 imputed CETP
## 12   rs3764261 1.430239e-07 6.84459142  16 56993324 imputed CETP
## 13    rs821840 1.430239e-07 6.84459142  16 56993886 imputed CETP
## 14  rs11508026 1.151771e-06 5.93863373  16 56999328 imputed CETP
## 15  rs12444012 1.151771e-06 5.93863373  16 57001438 imputed CETP

write.csv(CETP, CETP.fname, row.names=FALSE) # save for future use

save(genotype, clinical, pcs, imputed, target, rules, phenoSub, support, genes,
     impCETP, impCETPgeno, imputeOut, GWASout, GWAScomb, CETP, file=working.data.fname(9))

"Visualization and QC

Several plots allow us both to visualize the GWA analysis findings while performing quality control checks. Specifically, we are interested in identifying data inconsistencies and potential systemic biases.

Manhattan plot

Manhattan plots are used to visual GWA significant results by chromosome location. We will call the GWAS_Manhattan function to plot ???log10???log10 of the p-value against SNP position across the entire set of typed and imputed SNPs. The plot will show two horizontal lines. The higher of the two is the commonly used "Bonferroni" adjusted significance cut-off of ???log10(5×10???8)???log10(5×10???8), while the lower is less stringent ("Candidate") cut-off of ???log10(5×10???6)???log10(5×10???6). Typed and imputed SNPs will be represented by black and blue, respectively. We label the typed SNPs with signals that have surpassed the less stringent cutoff.
"
# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
                           col.detected=c("black"), col.imputed=c("blue"), col.text="black",
                           title="GWAS Tutorial Manhattan Plot", display.text=TRUE,
                           bonferroni.alpha=0.05, bonferroni.adjustment=1000000,
                           Lstringent.adjustment=10000) {
  
  bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
  Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment)
  xscale <- 1000000
  
  manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
  
  #sort the data by chromosome and then location
  manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
  manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
  
  ##Finding the maximum position for each chromosome
  max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
  max.pos2 <- c(0, cumsum(max.pos))                  
  
  #Add spacing between chromosomes
  max.pos2 <- max.pos2 + c(0:21) * xscale * 10
  
  #defining the positions of each snp in the plot
  manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
  
  # alternate coloring of chromosomes
  manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
  
  # draw the chromosome label roughly in the middle of each chromosome band
  text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
  
  # Plot the data
  plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
       pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab=NA,
       ylab="Negative Log P-value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
  #Add x-label so that it is close to axis
  mtext(side = 1, "Chromosome", line = 1.25)
  
  points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="imputed"],
         pch=20, cex=.4, col = col.imputed)
  
  points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
         pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
  
  axis(2)
  abline(h=0)
  
  SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
  
  #Add legend
  legend("topright",c("Bonferroni corrected threshold (p = 5E-8)", "Candidate threshold (p = 5E-6)"),
         border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1),
         lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.6)
  
  #Add chromosome number
  text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=.8)
  
  #Add bonferroni line
  abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
  
  #Add "less stringent" line
  abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
  
  #Plotting detected genes
  #Were any genes detected?
  if (length(SigNifSNPs)>0){
    
    sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
    
    points(manhat.ord$pos[sig.snps]/xscale,
           manhat.ord$Neg_logP[sig.snps],
           pch=20,col=col.detected, bg=col.detected,cex=0.5)
    
    text(manhat.ord$pos[sig.snps]/xscale,
         manhat.ord$Neg_logP[sig.snps],
         as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.5)
  }
}
# Create Manhattan Plot
GWAS_Manhattan(GWAScomb)
"plot of chunk code10-a

Quantile-quantile plots and the ????-statistic

Q-Q plots are used to visualize the relationship between the expected and observed distributions of SNP level test statistics. Here we compare these statistics for the unadjusted model (left) compared with the model adjusted for confounders by incorporating the first ten principal components along with clinical covariates.

A new set of models is generated with only the phenotype (HDL) and no additional factors. The results are plotted using the GenABEL package's estlambda function.
"
# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")] # remove all extra factors, leave only phenotype

GWAA(genodata=genotype, phenodata=phenoSub2, filename=gwaa.unadj.fname)
## 656890  SNPs included in analysis.
## 1309 samples included in analysis.
## socket cluster with 2 nodes on host 'localhost'
## GWAS SNPs 1-65689 (10% finished)
## GWAS SNPs 65690-131378 (20% finished)
## GWAS SNPs 131379-197067 (30% finished)
## GWAS SNPs 197068-262756 (40% finished)
## GWAS SNPs 262757-328445 (50% finished)
## GWAS SNPs 328446-394134 (60% finished)
## GWAS SNPs 394135-459823 (70% finished)
## GWAS SNPs 459824-525512 (80% finished)
## GWAS SNPs 525513-591201 (90% finished)
## GWAS SNPs 591202-656890 (100% finished)
## [1] "Done."
GWASoutUnadj <- read.table(gwaa.unadj.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
par(mfrow=c(1,2))
lambdaAdj <- estlambda(GWASout$t.value^2,plot=TRUE,method="median")
lambdaUnadj <- estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")
# plot of chunk code10-b

cat(sprintf("Unadjusted lambda: %s\nAdjusted lambda: %s\n", lambdaUnadj$estimate, lambdaAdj$estimate))
## Unadjusted lambda: 1.01417377078806
## Adjusted lambda: 1.00214021515846
# Calculate standardized lambda
lambdaAdj_1000<-1+(lambdaAdj$estimate-1)/nrow(phenoSub)*1000
lambdaUnadj_1000<-1+(lambdaUnadj$estimate-1)/nrow(phenoSub)*1000
cat(sprintf("Standardized unadjusted lambda: %s\nStandardized adjusted lambda: %s\n", lambdaUnadj_1000, lambdaAdj_1000))
## Standardized unadjusted lambda: 1.0108279379588
## Standardized adjusted lambda: 1.00163500012105
# We see here that the tail of the distribution is brought closer to the y=x line after accounting for confounding by race/ethnicity in the modeling framework. If the data in this figure were shifted up or down from the y=xy=x line, then we would want to investigate some form of systemic bias. The degree of deviation from this line is measured formally by the ????-statistic, where a value close to 1 suggests appropriate adjustment for the potential admixture. A slight deviation in the upper right tail from the y=xy=x line suggests crudely that some form of association is present in the data. There is only a slight improvement in ???? between the unadjusted model and the model with PCs indicating that the population is relatively homogenous.

#Heatmap

# Heatmaps are typically used in the context of GWA analysis to visualize the linkage disequilibrium pattern between significant SNPs other SNPs in nearby regions. Here we include our most significant SNP from our analysis and other SNPs near CETP. The darker shading indicates higher LD. The plot also includes ???log10(p)???log10(p) values to illustrate their connection with physical location and LD.

library(LDheatmap)
library(rtracklayer)

# Add "rs247617" to CETP
CETP <- rbind.fill(GWASout[GWASout$SNP == "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
subgen <- cbind(genotype[,colnames(genotype) %in% CETP$SNP], impCETPgeno)     # CETP subsets from typed and imputed SNPs

# Subset SNPs for only certain genotypes
certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c(0,1,2,NA)) })
subgen <- subgen[,certain]

# Subset and order CETP SNPs by position
CETP <- CETP[CETP$SNP %in% colnames(subgen),]
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames(subgen),CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen, subgen, stats="R.squared") # Find LD map of CETP SNPs

ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add plot of -log(p)
library(ggplot2)

plot.new()
llQplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .34)
pushViewport(viewport(x = 0.483, y= 0.76, width = .91 ,height = .4))

grid.draw(ggplotGrob({
qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", xlim = range(CETP$position),
asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) + 
theme(axis.text.x = element_blank(),
axis.title.y = element_text(size = rel(0.75)), legend.position = "none", 
panel.background = element_blank(), 
axis.line = element_line(colour = "black")) +
scale_color_manual(values = c("red", "black"))
}))
"plot of chunk code10-c

Regional Association

Similar to the LD heatmap above, a regional association plot allows for visualization of SNP-wise signal accross a segment of a particular chromsome with the pairwise correlation between SNPs. However regional assoication plots typically show a larger window of the genome. Therefore, for plot legibility, LD calculations to be displayed can be selected based on pairwise SNP proximity and minimum LD. In this example we demonstrate a regional plot create by the regionplot function from postgwas. This function can use HapMap data downloaded from Ensembl, for LD calculations. By default it will use the most recent Genome Reference Consortium human genome build. Therefore, if we wish to use build GRCh37 (hg19) we will have to create a custom biomartConfigs object to retrieve the appropriate data.
"
# Create regional association plot
# Create data.frame of most significant SNP only
library(postgwas)

snps<-data.frame(SNP=c("rs1532625"))

# Change column names necessary to run regionalplot function
GWAScomb <- rename(GWAScomb, c(p.value="P", chr="CHR", position="BP"))


# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hg19

myconfig <- biomartConfigs$hsapiens
myconfig$hsapiens$gene$host <- "grch37.ensembl.org"
myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL"
myconfig$hsapiens$snp$host <- "grch37.ensembl.org"
myconfig$hsapiens$snp$mart <- "ENSEMBL_MART_SNP"


# Run regionalplot using HAPMAP data (pop = CEU)
regionalplot(snps, GWAScomb, biomart.config = myconfig, window.size = 400000, draw.snpname = data.frame(
snps = c("rs1532625", "rs247617"), 
text = c("rs1532625", "rs247617"),
angle = c(20, 160),
length = c(1, 1), 
cex = c(0.8)
),
ld.options = list(
gts.source = 2, 
max.snps.per.window = 2000, 
rsquare.min = 0.8, 
show.rsquare.text = FALSE
),
out.format = list(file = NULL, panels.per.page = 2))