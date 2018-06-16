{df = matrix(nrow = x, ncol = y)
df = as.data.frame(df)
i = 1
for (i in 1:y){
  n = i*2
  test2[,(n-1):n] = replicate(2,df[,i])
  print(i)
}}
rownames(df) = rownames(df)

write.table(df, file = "", col.names = F, row.names = F, quote = F)
