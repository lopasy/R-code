p = ggplot(z, aes(factor(CNV_22q11.2dup), AgeSpexWear, fill=factor(CNV_22q11.2dup))) +
  geom_violin(trim = FALSE, width = 1.1) + 
  geom_boxplot(width = 0.08,outlier.shape=NA, fill = "white") +
  scale_fill_brewer(palette="Blues") +
  guides(fill=FALSE) +
  labs(x = "CNV_22q11.2dup", y = "Years") +
  theme_classic() +
  theme(axis.text= element_text(size=10,colour = "black"), axis.title=element_text(size=11,face="bold"))
