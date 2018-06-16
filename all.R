library(readr)
SNP8YG <- read_delim("G:/Work/Epistasis/SNP_5e-08YG.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP7YG <- read_delim("G:/Work/Epistasis/SNP_5e-07YG.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP6YG <- read_delim("G:/Work/Epistasis/SNP_5e-06YG.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP5YG <- read_delim("G:/Work/Epistasis/SNP_5e-05YG.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP4YG <- read_delim("G:/Work/Epistasis/SNP_5e-04YG.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP3YG <- read_delim("G:/Work/Epistasis/SNP_5e-03YG.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP2YG <- read_delim("G:/Work/Epistasis/SNP_5e-02YG.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

SNP8VH <- read_delim("G:/Work/Epistasis/SNP_5e-08VH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP7VH <- read_delim("G:/Work/Epistasis/SNP_5e-07VH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP6VH <- read_delim("G:/Work/Epistasis/SNP_5e-06VH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP5VH <- read_delim("G:/Work/Epistasis/SNP_5e-05VH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP4VH <- read_delim("G:/Work/Epistasis/SNP_5e-04VH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP3VH <- read_delim("G:/Work/Epistasis/SNP_5e-03VH.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP2VH <- read_delim("G:/Work/Epistasis/SNP_5e-02VH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)

SNP8YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-08YGVH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP7YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-07YGVH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP6YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-06YGVH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP5YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-05YGVH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP4YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-04YGVH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP3YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-03YGVH.epi.qt", "\t", escape_double = FALSE, trim_ws = TRUE)
SNP2YGVH <- read_delim("G:/Work/Epistasis/SNP_5e-02YGVH.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


{par(mfrow=c(3,2))
  plot(SNP2YG$BETA_INT, -log10(SNP2YG$P), main='YG', ylab = substitute(-log[10](pval)),
       xlab = expression(hat(beta)[GxG]), cex.lab = 1.3, col = "yellowgreen", pch = 16)
  points(SNP3YG$BETA_INT, -log10(SNP3YG$P), main='YG', ylab = substitute(-log[10](pval)),
    xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "slategrey", pch = 16)
  points(SNP4YG$BETA_INT, -log10(SNP4YG$P), main='YG', ylab = substitute(-log[10](pval)),
    xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "dodgerblue", pch = 16)
  points(SNP5YG$BETA_INT, -log10(SNP5YG$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
       xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "forestgreen", pch = 16)
  points(SNP6YG$BETA_INT, -log10(SNP6YG$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "firebrick3", pch = 16)
  points(SNP7YG$BETA_INT, -log10(SNP7YG$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "gold1", pch = 16)
  points(SNP8YG$BETA_INT, -log10(SNP8YG$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "darkviolet", pch = 16)
  
  
  abline(h = -log10(0.000089), col="darkviolet", lwd=3, lty=1)
  text(-0.0025,4.3,expression("8.9x10"^"-5"))
  abline(h = -log10(0.0000344), col="gold1", lwd=4, lty=2)
  text(-0.0025,4.55,expression("3.44x10"^"-5"))
  abline(h = -log10(0.00000949), col="firebrick3", lwd=3, lty=1)
  text(-0.0025,5.25,expression("9.49x10"^"-6"))
  abline(h = -log10(0.000001539), col="forestgreen", lwd=4, lty=2)
  text(-0.0025,6.1,expression("1.539x10"^"-6"))
  abline(h = -log10(0.0000001169), col="dodgerblue", lwd=3, lty=1)
  text(-0.0025,7.2,expression("1.169x10"^"-7"))
  abline(h = -log10(0.000000004294), col="slategrey", lwd=4, lty=2)
  text(-0.0025,8.6,expression("4.294x10"^"-9"))
  abline(h = -log10(0.0000000005), col="yellowgreen", lwd=3, lty=1)
  text(-0.0025,9.35,expression("5x10"^"-10"))}


{par(mfrow=c(3,2))
  plot(SNP2YGVH$BETA_INT, -log10(SNP2YGVH$P), main='YGVH', ylab = substitute(-log[10](pval)),
       xlab = expression(hat(beta)[GxG]), cex.lab = 1.3, col = "yellowgreen", pch = 16)
  points(SNP3YGVH$BETA_INT, -log10(SNP3YGVH$P), main='YGVH', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "slategrey", pch = 16)
  points(SNP4YGVH$BETA_INT, -log10(SNP4YGVH$P), main='YGVH', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "dodgerblue", pch = 16)
  points(SNP5YGVH$BETA_INT, -log10(SNP5YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "forestgreen", pch = 16)
  points(SNP6YGVH$BETA_INT, -log10(SNP6YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "firebrick3", pch = 16)
  points(SNP7YGVH$BETA_INT, -log10(SNP7YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "gold1", pch = 16)
  points(SNP8YGVH$BETA_INT, -log10(SNP8YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "darkviolet", pch = 16)
  
  
  abline(h = -log10(0.000037), col="darkviolet", lwd=3, lty=1)
  #text(-0.0025,4.3,expression("3.7x10"^"-5"))
  abline(h = -log10(0.0000111), col="gold1", lwd=4, lty=2)
  #text(-0.0025,4.55,expression("1.11x10"^"-5"))
  abline(h = -log10(0.00000224), col="firebrick3", lwd=3, lty=1)
  #text(-0.0025,5.25,expression("2.24x10"^"-6"))
  abline(h = -log10(0.000000234), col="forestgreen", lwd=4, lty=2)
  #text(-0.0025,6.1,expression("2.34x10"^"-7"))
  abline(h = -log10(0.0000000157), col="dodgerblue", lwd=3, lty=1)
  #text(-0.0025,7.2,expression("1.57x10"^"-8"))
  abline(h = -log10(0.000000000748), col="slategrey", lwd=4, lty=2)
  #text(-0.0025,8.6,expression("7.48x10"^"-10"))
  abline(h = -log10(0.00000000005), col="yellowgreen", lwd=3, lty=1)
  #text(-0.0025,9.35,expression("5x10"^"-11"))}
  
  
  plot(SNP2VH$BETA_INT, -log10(SNP2VH$P), main='VH', ylab = substitute(-log[10](pval)),
       xlab = expression(hat(beta)[GxG]), cex.lab = 1.3, col = "yellowgreen", pch = 16)
  points(SNP3VH$BETA_INT, -log10(SNP3VH$P), main='YGVH', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "slategrey", pch = 16)
  points(SNP4VH$BETA_INT, -log10(SNP4VH$P), main='YGVH', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "dodgerblue", pch = 16)
  points(SNP5VH$BETA_INT, -log10(SNP5VH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "forestgreen", pch = 16)

  
  
  abline(h = -log10(0.000037), col="darkviolet", lwd=3, lty=1)
  #text(-0.0025,4.3,expression("3.7x10"^"-5"))
  abline(h = -log10(0.0000111), col="gold1", lwd=4, lty=2)
  #text(-0.0025,4.55,expression("1.11x10"^"-5"))
  abline(h = -log10(0.00000224), col="firebrick3", lwd=3, lty=1)
  #text(-0.0025,5.25,expression("2.24x10"^"-6"))
  abline(h = -log10(0.000000234), col="forestgreen", lwd=4, lty=2)
  #text(-0.0025,6.1,expression("2.34x10"^"-7"))
  abline(h = -log10(0.0000000157), col="dodgerblue", lwd=3, lty=1)
  #text(-0.0025,7.2,expression("1.57x10"^"-8"))
  abline(h = -log10(0.000000000748), col="slategrey", lwd=4, lty=2)
  #text(-0.0025,8.6,expression("7.48x10"^"-10"))
  abline(h = -log10(0.00000000005), col="yellowgreen", lwd=3, lty=1)
  #text(-0.0025,9.35,expression("5x10"^"-11"))
  
  
#  hist(SNP4YG$P, 200, col='black', xlab = "P-value", main = "Histogram of P values", cex.lab = 1.2)
#  qqPlot(SNP4YG$P, main='QQ-plot', ylab = substitute(-log[10]("observed P")),
#         xlab = substitute(-log[10]("expected P")))
  plot(-log10(SNP4YG$P), main='YG', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "dodgerblue", pch = 16)
  points(-log10(SNP5YG$P), main='Manhattan plot', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "forestgreen", pch = 16)
  points(-log10(SNP6YG$P), main='Manhattan plot', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "firebrick3", pch = 16)
  points(-log10(SNP7YG$P), main='Manhattan plot', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "gold1", pch = 16)
  points(-log10(SNP8YG$P), main='Manhattan plot', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "darkviolet", pch = 16)
  
  abline(h = -log10(0.0005), col="blue", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-5"))
  abline(h = -log10(0.00005), col="dodgerblue", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-6"))
  abline(h = -log10(0.000005), col="forestgreen", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-7"))
  abline(h = -log10(0.00000005), col="firebrick3", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-8"))
  abline(h = -log10(0.000000005), col="gold1", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-9"))
  abline(h = -log10(0.0000000005), col="darkviolet", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-10"))
  arrows(720,9.19,800,9.33, length=0.1, lwd = 2)
  text(680,9.10,"chr4 rs73244033\nchr6 rs17196517", cex = 0.9)
  arrows(450,9.16,510,9.06, length=0.1, lwd = 2)
  text(400,9.24,"Intron in DSE gene", cex = 0.9)
#  hist(table(SNP4YG$SNP1), 20, main = "Frequency of SNP 1 ", cex.lab = 1.2, xlab = "Number of interactions")
#  hist(table(SNP4YG$SNP2), 20, main = "Frequency of SNP 2 ", cex.lab = 1.2, xlab = "Number of interactions")
  
  
  
  plot(SNP4VH$BETA_INT, -log10(SNP4VH$P), main='VH', ylab = substitute(-log[10](pval)),
       xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "dodgerblue", pch = 16)
  points(SNP5VH$BETA_INT, -log10(SNP5VH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "forestgreen", pch = 16)
  
  
  abline(h = -log10(0.0005), col="blue", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-5"))
  abline(h = -log10(0.00005), col="dodgerblue", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-6"))
  abline(h = -log10(0.000005), col="forestgreen", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-7"))
  abline(h = -log10(0.00000005), col="firebrick3", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-8"))
  abline(h = -log10(0.000000005), col="gold1", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-9"))
  abline(h = -log10(0.0000000005), col="darkviolet", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-10"))
  
  
  #  hist(SNP4VH$P, 200, col='black', xlab = "P-value", main = "Histogram of P values", cex.lab = 1.2)
  #  qqPlot(SNP4VH$P, main='QQ-plot', ylab = substitute(-log[10]("observed P")),
  #         xlab = substitute(-log[10]("expected P")))
  plot(-log10(SNP4VH$P), main='VH', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "dodgerblue", pch = 16)
  points(-log10(SNP5VH$P), main='Manhattan plot', xlab = "GxG Comparison Number",
         ylab = '', cex.lab = 1.2, col = "forestgreen", pch = 16)
  
  
  abline(h = -log10(0.0005), col="blue", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-5"))
  abline(h = -log10(0.00005), col="dodgerblue", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-6"))
  abline(h = -log10(0.000005), col="forestgreen", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-7"))
  abline(h = -log10(0.00000005), col="firebrick3", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-8"))
  abline(h = -log10(0.000000005), col="gold1", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-9"))
  abline(h = -log10(0.0000000005), col="darkviolet", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-10"))
  arrows(720,9.19,800,9.33, length=0.1, lwd = 2)
  text(680,9.10,"chr4 rs73244033\nchr6 rs17196517", cex = 0.9)
  arrows(450,9.16,510,9.06, length=0.1, lwd = 2)
  text(400,9.24,"Intron in DSE gene", cex = 0.9)
  
  
  plot(SNP4YGVH$BETA_INT, -log10(SNP4YGVH$P), main='YGVH', ylab = substitute(-log[10](pval)),
       xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "dodgerblue", pch = 16)
  points(SNP5YGVH$BETA_INT, -log10(SNP5YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "forestgreen", pch = 16)
  points(SNP6YGVH$BETA_INT, -log10(SNP6YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "firebrick3", pch = 16)
  points(SNP7YGVH$BETA_INT, -log10(SNP7YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "gold1", pch = 16)
  points(SNP8YGVH$BETA_INT, -log10(SNP8YGVH$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
         xlab = expression(hat(beta)[GxE]), cex.lab = 1.3, col = "darkviolet", pch = 16)
  
  
  abline(h = -log10(0.0005), col="blue", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-5"))
  abline(h = -log10(0.00005), col="dodgerblue", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-6"))
  abline(h = -log10(0.000005), col="forestgreen", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-7"))
  abline(h = -log10(0.00000005), col="firebrick3", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-8"))
  abline(h = -log10(0.000000005), col="gold1", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-9"))
  abline(h = -log10(0.0000000005), col="darkviolet", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-10"))
  
  
  #  hist(SNP4YGVH$P, 200, col='black', xlab = "P-value", main = "Histogram of P values", cex.lab = 1.2)
  #  qqPlot(SNP4YGVH$P, main='QQ-plot', ylab = substitute(-log[10]("observed P")),
  #         xlab = substitute(-log[10]("expected P")))
  plot(-log10(SNP4YGVH$P), main='YGVH', xlab = "GxG Comparison Number",
       ylab = '', cex.lab = 1.2, col = "dodgerblue", pch = 16)
  points(-log10(SNP5YGVH$P), main='Manhattan plot', xlab = "GxG Comparison Number",
         ylab = '', cex.lab = 1.2, col = "forestgreen", pch = 16)
  points(-log10(SNP6YGVH$P), main='Manhattan plot', xlab = "GxG Comparison Number",
         ylab = '', cex.lab = 1.2, col = "firebrick3", pch = 16)
  points(-log10(SNP7YGVH$P), main='Manhattan plot', xlab = "GxG Comparison Number",
         ylab = '', cex.lab = 1.2, col = "gold1", pch = 16)
  points(-log10(SNP8YGVH$P), main='Manhattan plot', xlab = "GxG Comparison Number",
         ylab = '', cex.lab = 1.2, col = "darkviolet", pch = 16)
  
  abline(h = -log10(0.0005), col="blue", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-5"))
  abline(h = -log10(0.00005), col="dodgerblue", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-6"))
  abline(h = -log10(0.000005), col="forestgreen", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-7"))
  abline(h = -log10(0.00000005), col="firebrick3", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-8"))
  abline(h = -log10(0.000000005), col="gold1", lwd=4, lty=2)
  text(20,8.4,expression("5x10"^"-9"))
  abline(h = -log10(0.0000000005), col="darkviolet", lwd=3, lty=1)
  text(24,9.39,expression("5x10"^"-10"))
  arrows(720,9.19,800,9.33, length=0.1, lwd = 2)
  text(680,9.10,"chr4 rs73244033\nchr6 rs17196517", cex = 0.9)
  arrows(450,9.16,510,9.06, length=0.1, lwd = 2)
  text(400,9.24,"Intron in DSE gene", cex = 0.9)}

par(mfrow=c(1,1))



