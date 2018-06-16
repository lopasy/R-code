# Split dataset by row
predicted = read.csv("~/Documents/Useless/master_predicted.txt", sep="")
l = as.data.frame(sample(predicted[,1], 160228, F))
l$du = l
n_cv = 5
split = nrow(pred)/n_cv
l = split(pred, rep(1:n_cv, each = split))

# Save
lapply(names(l), function(x, l) write.table(l[[x]], paste("cv", x, ".txt", sep = ""), col.names=F, 
                                                row.names=F, sep="\t", quote=F), l)






