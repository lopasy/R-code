reps = 1000
B = 5

N = 6000
SNPs = 10000
qtl = 0.1635
freq = runif(6000, 0.1, 0.4)
env1 = env2 = 0.3
marginalG = marginalE1 = 0.0017
marginalGE1 = seq(0,0.004, 0.002)
marginalE2 = marginalGE2 = 0.002

marginalG = effectG^2*(qtl(2-qtl)(1-qtl)^2)
marginalE1 = effectE1^2*(env1(1-env1))
marginalGE1 = effectGE1^2*(2*env1(1-env1)(qtl(2-qtl)(1-qtl)^2))
marginalE2 = effectE2^2*(env2(1-env2))
marginalGE2 = effectGE2^2*(2*env2(1-env2)(qtl(2-qtl)(1-qtl)^2))
variance_scenario1 = 1-(marginalG+marginalE1+marginalGE1)
variance_scenario2 = 1-(marginalG+marginalE1+marginalGE1+marginalE2+marginalGE2)
error1 = rnorm(6000, 0, variance_scenario1)
error2 = rnorm(6000, 0, variance_scenario2)

meanG = qtl(2-qtl)


geno = matrix(rbinom(6000*10000, 1, 0.235), ncol = 6000, nrow = 10000)

