plot(-log10(final$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
     ylab = substitute(-log[10](pval)), cex.lab = 1.2, col = "red", ylim = range(3:10), type = "p")
points(-log10(epistasis4$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
     ylab = substitute(-log[10](pval)), cex.lab = 1.2, col = "blue")
points(-log10(epistasis5$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
     ylab = substitute(-log[10](pval)), cex.lab = 1.2, col = "green")
points(-log10(epistasis6$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
       ylab = substitute(-log[10](pval)), cex.lab = 1.2, col = "black")
points(-log10(epistasis7$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
       ylab = substitute(-log[10](pval)), cex.lab = 1.2, col = "violet")
points(-log10(epistasis8$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
       ylab = substitute(-log[10](pval)), cex.lab = 1.2, col = "black")

ggplot(epistasis4,aes(BETA_INT,-log10(P)))+geom_point(aes(color=""))+ 
  geom_point(data=epistasis5,aes(color="Second line"))+
  geom_point(data=epistasis6,aes(color="Third line"))+
  geom_point(data=epistasis7,aes(color="Fourth line"))+
  geom_point(data=epistasis8,aes(color="Fifth line"))+
  labs(color="Legend")

