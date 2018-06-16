library(data.table); library(qqman); library(ggplot2);library(gmodels)
{chr = list.files(pattern = ".log", full.names = TRUE)
chr = lapply(chr, fread)
d = do.call(rbind.data.frame, chr)
}


s = t[grep("min lambda", t$PLINK.v1.90b3c.64.bit..2.Feb.2015.),]
write.csv(s, "t.csv", quote = F, row.names = F)


a = as.data.frame(noquote(t$b.2015.))
a = gsub("([0-9].+)*.", "\\1", a$`noquote(t$b.2015.)`)
a = as.numeric(a)

{h = hist(a, breaks = 60)
xfit = seq(min(a), max(a), length = 40) 
yfit = dnorm(xfit, mean=mean(a), sd=sd(a)) 
yfit = yfit*diff(h$mids[1:2])*length(a) 
lines(xfit, yfit, col = "blue", lwd = 2)
summary(a)}
ci(a)
