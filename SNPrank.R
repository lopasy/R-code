# SNPrank
library(stringr)
YG$SNP = str_split_fixed(YG$SNP, "_",3)

YG2 = YG[1:10,]

colnames(merged) = c("CHR","SNP","S","POS","REF", "EFF")
a = merge(YG2,merged, by="SNP")
barplot(a1$SNPrank, a1$CHR, xlab = "10")
barplot(a2$SNPrank, a2$CHR, add = T, xlab = "13")

library(ggplot2)
ggplot(data=a, aes(x=interaction(CHR,POS, lex.order = T), y=SNPrank, colour = CHR)) +
  coord_cartesian(ylim = c(0.000598, 0.00125), expand = F) +
#  scale_x_discrete(a$CHR)+
#  annotate(geom = "text", x = a$CHR, y = 0, label = a$CHR, size = 4) +
  geom_bar(stat="identity", position = "identity", fill = NA) +
  theme(legend.position = "none")

g2 <- ggplot_gtable(ggplot_build(a))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid::grid.draw(g2)