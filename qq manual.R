library(plyr)
library(reshape)
library(ggplot2)
library(gridExtra)

My.Qq.log <- function(pvalues) {
  # pvalues has to be a vector (e.g. a column in a data frame)
  observed <- -log10(sort(pvalues, decreasing = FALSE))
  expected <- -log10(1:length(observed) / length(observed) )
  out <- data.frame(expected = expected,
                    observed = observed)
}
plot(My.Qq.log(cassi$P), pch = 19,xlab = expression(Expected ~ ~ -log[10](italic(p))),
     ylab = expression(Observed ~ ~ -log[10](italic(p))))
abline(0,1, col = "red")
points(My.Qq.log(YG_3_0$LIN_COVAR_P), pch = 19, col = "blue")



qq(cassi$P)

theme_my <- theme_bw() + theme(
  axis.line        = element_line(colour = "black"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.border     = element_blank(),
  strip.background = element_blank(),
  legend.key       = element_blank()
)


MyQQPlots <- function(data, name = "Name") {
  # default for the column that holds the data set name is "Name"
  require(ggplot2)
  current.Name <- unique(data[, names(data) == name]) # get the Name of the current subset
  
  # axes are hard-coded to "expected" and "observed"
  p <- ggplot(data, aes(x = expected, y = observed))
  p <- p + geom_abline(intercept  = 0, 
                       slope  = 1, 
                       colour = "black", linetype = "dashed")
  p <- p + geom_point(aes(colour  = variable))
  p <- p + scale_colour_manual(
    name = "", 
    values = c("red", "green", "blue"),
    labels = c("LM", "LM (X)", "GLM (X+K)"))
  p <- p + labs(x = expression(Expected~~-log[10](italic(p))),
                y = expression(Observed~~-log[10](italic(p))),
                title = current.Name)
  p <- p + theme_my 
  p <- p + theme(legend.position = "none")
  return(p)  
}


cassi2 = cassi[,12]
cassi2$new = "cassi"
YG_5_0_2_250_epi$new = "plink"
colnames(YG_5_0_2_250_epi) = colnames(cassi2)
cassi3 = rbind(cassi2,YG_5_0_2_250_epi)
cassi3 = as.data.frame(cassi3)
my.melt <- melt(cassi3, id = "new")
my.melt = my.melt[which(my.melt$variable == "LIN_COVAR_P"),]


my.qq.log <- ddply(my.melt,
                   .(new,variable),
                   function(x) My.Qq.log(cassi3$LIN_COVAR_P))
my.plots <- dlply(my.qq.log,
                  .(new,variable),        
                  function(x) MyQQPlots(x))
do.call("grid.arrange", c(my.plots, ncol = 2))


