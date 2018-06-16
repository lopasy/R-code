library(msm)
cream$cc_0[cream$mse_at_visit <= 0]= 1; cream$cc_0[cream$mse_at_visit > 0]= 0

m =glm(cc_0 ~ rs12898755_A*ReadingBinary + sex + age_at_visit, data = cream,  family=binomial(link="logit"))
summary(m)
additive_interactions(m, cream, recode = T,  monotone = 2)
