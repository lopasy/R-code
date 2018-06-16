na <-sapply(testing, function(y) sum(length(which(is.na(y)))))
na <- data.frame(na)



m <- "intialize" # initialize, sometimes better as just `list()`
for(i in 2:ncol(testing)-1){
  tryCatch(if(table(testing[i])[1] > table(testing[i])[2] & table(testing[i])[1] > table(testing[i])[3]){
    # paste into position i of vector m
    testing[i][is.na(testing[i])] <- 0
  } 
    else if(table(testing[i])[2] > table(testing[i])[1] & table(testing[i])[2] > table(testing[i])[3]){
    # paste into position i of vector m
    testing[i][is.na(testing[i])] <- 1
  } 
    else {
    testing[i][is.na(testing[i])] <- 2
  },error=function(error_message) {
    message(error_message)
    return(NULL)})
}
m

summary(lm(V1~.,testing))
summary(testing)


mean2 <-sapply(testin2, function(y) mean(y, na.rm = T))
mean <-sapply(testing, function(y) mean(y, na.rm = T))
plot(density(mean2))
lines(density(mean))

write.table(testing, "112.dat", quote = F, col.names = F, row.names = F)
