# Required libraries
library(ggplot2)


# Defaults
N = 100
maf = 0.05
SNPs = 10
num.phen = 2
prevalence = 0.3
pheno.name = "avMSE"
mean = 0
variance = 1
#########################################
############# Phenotypes ################
#########################################

### Univariate distributions
disease = function(N, prevalence, num.phen){
  empty = matrix(nrow = N, ncol = num.phen)
  samples = replicate(num.phen, rbinom(N, 1, prevalence))
  samples = as.data.frame(samples) # For if only trait selected
}

continious = function(N, mean, variance, num.phen){
  empty = matrix(nrow = N, ncol = num.phen)
  samples = replicate(num.phen, rnorm(N, mean, variance))
  samples = as.data.frame(samples) # For if only trait selected
}

### Multivariate distributions ### So far only mixture of normals normal
continious.mvn = function(N, n.dist, prob.dist1, mean1, sd1, num.phen) {
  components = matrix(nrow = N, ncol = num.phen)
  components <- sample(1:n.dist,prob=c(prob.dist1),size=N,replace=TRUE)
  mus <- c(mean1)
  sds <- sqrt(c(sd1))
  samples <- replicate(num.phen,rnorm(n=N,mean=mus[components],sd=sds[components]))
  samples = as.data.frame(samples)
}
#########################################
############# Genotypes #################
#########################################

# Under a random mating assumption, the number of alternative alleles for a SNP and individual follows a binomial distribution with 2 draws (one from mum and one from dad) and probability equal to the minor allele frequency of that SNP.


# Let's assume that the minor allele frequency of our MM SNPs come from a uniform distribution between 0 and 0.5
maf = function(N, n1, n2){
  runif(N, n1, n2)
}

### The following generates unscaled matrix
genotype = function(N, MAF, SNPs, MONOMORPHIC){
  if(MONOMORPHIC == FALSE){
  genotypes = t(replicate(N, rbinom(10*SNPs, 2, c(MAF, MAF))))
  genotypes.polymorphic = apply(genotypes, 2, var) > 0
  genotypes = genotypes[,genotypes.polymorphic][,1:SNPs]
  # To check if polymorphic
  # A1 = c(MAF, MAF)[genotypes.polymorphic][1:SNPs]
  # round(A1[1:10], 2)
  } else {
    genotypes = t(replicate(N, rbinom(SNPs, 2, c(MAF, MAF))))
  # To check if polymorphic
    # round(amaf[1:10], 2) # Ignore amaf. It's a placeholder for MAF obtained with maf function
  }
}
# To check if MONOMORPHIC or not
round(A1[1:10], 2)

### The following generates unscaled matrix
multiallelic = function(N, MAF, SNPs, MONOMORPHIC,n.alleles){
  if(MONOMORPHIC == FALSE){
    multialleles = t(replicate(N, rbinom(10*SNPs, n.alleles, c(MAF, MAF))))
    multialleles.polymorphic = apply(multialleles, 2, var) > 0
    multialleles = multialleles[,multialleles.polymorphic][,1:SNPs]
    # To check if polymorphic
    # A1 = c(MAF, MAF)[multialleles.polymorphic][1:SNPs]
    # round(A1[1:10], 2)
  } else {
    multialleles = t(replicate(N, rbinom(SNPs, n.alleles, c(MAF, MAF))))
    # To check if polymorphic
    # round(amaf[1:10], 2)
  }
}
# To check if MONOMORPHIC or not
round(A1[1:10], 2)

amaf = round(amaf[1:10], 2)
### Estimate MAF for uknown cases
maf_est = colMeans(geno)/2
qplot(A1, maf_est, col=amaf) +
  geom_abline() + xlab('p') + ylab(expression(hat(p)))



#########################################
####### Environment or Covariates #######
#########################################
cov = function(N, prevalence, num.phen){
  empty = matrix(nrow = N, ncol = num.phen)
  samples = replicate(num.phen, rbinom(N, 1, prevalence))
  samples = as.data.frame(samples) # For if only trait selected
}

covs = function(N, mean, variance, num.phen){
  empty = matrix(nrow = N, ncol = num.phen)
  samples = replicate(num.phen, rnorm(N, mean, variance))
  samples = as.data.frame(samples) # For if only trait selected
}





############ Data arrangment ############

### This is only for univariate
create.frame = function(pheno, geno, pheno.name, SNPs, cov, covs){
  if(cov == T | covs == T){
    colnames(pheno)[1] = pheno.name
    colnames(geno)[1:SNPs] <- paste("rs", 1:SNPs, sep = "")
    dataset = as.data.frame(cbind(pheno, geno))
    colnames(cov)[1:ncol(cov)] <- paste("Cov", 1:ncol(cov), sep = "")
    dataset = as.data.frame(cbind(dataset,cov))
    colnames(covs)[1:ncol(covs)] <- paste("Covariate", 1:ncol(covs), sep = "")
    dataset = as.data.frame(cbind(dataset,covs))
  }else if(cov == T){
    colnames(pheno)[1] = pheno.name
    colnames(geno)[1:SNPs] <- paste("rs", 1:SNPs, sep = "")
    dataset = as.data.frame(cbind(pheno, geno))
    colnames(cov)[1:ncol(cov)] <- paste("Cov", 1:ncol(cov), sep = "")
    dataset = as.data.frame(cbind(dataset,cov))
  }else if(covs == T){
    colnames(pheno)[1] = pheno.name
    colnames(geno)[1:SNPs] <- paste("rs", 1:SNPs, sep = "")
    dataset = as.data.frame(cbind(pheno, geno))
    colnames(covs)[1:ncol(covs)] <- paste("Covariate", 1:ncol(covs), sep = "")
    dataset = as.data.frame(cbind(dataset,covs))
  }else{
    colnames(pheno)[1] = pheno.name
    colnames(geno)[1:SNPs] <- paste("rs", 1:SNPs, sep = "")
    dataset = as.data.frame(cbind(pheno, geno))
  }
}


dataset

