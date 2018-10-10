num_snps = 100
N = 1000
MAF = 0.05


snp_effect_ols = -0.05
snp_effect_cqr = rnorm(9, -0.05, 1/num_snps)

x012 = t(replicate(N, rbinom(num_snps, 2, MAF)))
names(x012) = 