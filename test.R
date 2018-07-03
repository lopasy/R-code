{test1 = lme(mse_at_visit ~ rs75770582_CAT + poly(I(age_at_visit - 7.5),2) +
            rs75770582_CAT:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))

test2 = lme(mse_at_visit ~ rs75770582_CAT + poly(I(age_at_visit - 7.5),3) +
            rs75770582_CAT:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))

test3 = lme(mse_at_visit ~ rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
           rs75770582_CAT:I(age_at_visit - 7.5),
           random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
           na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))}

anova(test1, test2, test3)
summary(test3)
 
test4 = lme(mse_at_visit ~ sex + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))   

anova(test3, test4)

{test5 = lme(mse_at_visit ~ ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))   

test6 = lme(mse_at_visit ~ ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:ReadingBinary,
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))  

test7 = lme(mse_at_visit ~ ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 

test8 = lme(mse_at_visit ~ ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5) +
            rs75770582_CAT:ReadingBinary:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))}

anova(test5, test6, test7, test8)
summary(test7)

test9 = lme(mse_at_visit ~ sex + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 

anova(test7, test9)




test10 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
            rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 

summary(test10)


test11 = lme(mse_at_visit ~ rs511217_T + Outdoor1171 + poly(I(age_at_visit - 7.5),4) +
               rs511217_T:Outdoor1171 + Outdoor1171:I(age_at_visit - 7.5) +rs511217_T:I(age_at_visit - 7.5) +
             rs511217_T:Outdoor1171:I(age_at_visit - 7.5),
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 
anova(test10, test11)

test12 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
             rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5) +
             rs1150687_C:ReadingBinary + rs1150687_C:I(age_at_visit - 7.5) +
             rs75770582_CAT:ReadingBinary:I(age_at_visit - 7.5),
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 

test13 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
             rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5) +
             rs1150687_C:ReadingBinary + rs1150687_C:I(age_at_visit - 7.5) +
             rs75770582_CAT:ReadingBinary:I(age_at_visit - 7.5) +
             rs1150687_C:ReadingBinary:I(age_at_visit - 7.5),
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 
anova(test10,test11,test12,test13)
anova(test10,test11)
anova(test10,test12)
anova(test10,test13)

test10 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
               rs75770582_CAT:ReadingBinary + rs75770582_CAT:I(age_at_visit - 7.5) + ReadingBinary:I(age_at_visit - 7.5),
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 

summary(test10)
#############
# SNP - SNP #
#############
test1 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
               rs1150687_C:rs75770582_CAT,
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 

test2 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
              rs1150687_C:rs75770582_CAT + rs1150687_C:ReadingBinary + rs75770582_CAT:ReadingBinary,
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))
anova(test1, test2)

test3 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
              rs1150687_C:rs75770582_CAT + rs1150687_C:I(age_at_visit - 7.5) + rs75770582_CAT:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1))
anova(test1, test3)

test4 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
               rs1150687_C:rs75770582_CAT + rs1150687_C:ReadingBinary + rs75770582_CAT:ReadingBinary +
               rs1150687_C:rs75770582_CAT:ReadingBinary,
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 
anova(test1, test4)

test5 = summary(lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
               rs1150687_C:rs75770582_CAT + rs1150687_C:I(age_at_visit - 7.5) + rs75770582_CAT:I(age_at_visit - 7.5) +
               rs1150687_C:rs75770582_CAT:I(age_at_visit - 7.5),
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) )
anova(test1, test5)
anova(test1, test2, test3, test4, test5)



test6 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
             rs1150687_C:rs75770582_CAT + rs1150687_C:I(age_at_visit - 7.5) + rs75770582_CAT:I(age_at_visit - 7.5) + rs1150687_C:rs75770582_CAT + rs1150687_C:ReadingBinary + ReadingBinary:I(age_at_visit - 7.5),
             random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
             na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 
anova(test5, test6)

test7 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
              rs1150687_C:rs75770582_CAT + rs1150687_C:I(age_at_visit - 7.5) + rs75770582_CAT:I(age_at_visit - 7.5) + rs1150687_C:rs75770582_CAT + rs1150687_C:ReadingBinary + ReadingBinary:I(age_at_visit - 7.5) +
              rs1150687_C:rs75770582_CAT:ReadingBinary + rs1150687_C:rs75770582_CAT:I(age_at_visit - 7.5) + rs1150687_C:ReadingBinary:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 
anova(test6, test7)

test8 = lme(mse_at_visit ~ rs1150687_C + ReadingBinary + rs75770582_CAT + poly(I(age_at_visit - 7.5),4) +
              rs1150687_C:rs75770582_CAT + rs1150687_C:I(age_at_visit - 7.5) + rs75770582_CAT:I(age_at_visit - 7.5) + rs1150687_C:rs75770582_CAT + rs1150687_C:ReadingBinary + ReadingBinary:I(age_at_visit - 7.5) +
              rs1150687_C:rs75770582_CAT:ReadingBinary + rs1150687_C:rs75770582_CAT:I(age_at_visit - 7.5) + rs1150687_C:ReadingBinary:I(age_at_visit - 7.5) +
              rs1150687_C:rs75770582_CAT:ReadingBinary:I(age_at_visit - 7.5),
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, method = "ML", cor=corAR1(form=~ visit | alfred_ID1)) 
anova(test7, test8)














