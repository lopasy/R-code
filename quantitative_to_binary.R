# For MSE threshold of 0

a = ifelse(gcta$X3 > 0, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_1.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.1

a = ifelse(gcta$X3 > -0.1, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.1.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.2

a = ifelse(gcta$X3 > -0.2, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.2.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.3

a = ifelse(gcta$X3 > -0.3, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.3.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.4

a = ifelse(gcta$X3 > -0.4, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.4.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.5

a = ifelse(gcta$X3 > -0.5, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.5.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.6

a = ifelse(gcta$X3 > -0.6, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.6.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.7

a = ifelse(gcta$X3 > -0.7, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.7.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.8

a = ifelse(gcta$X3 > -0.8, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.8.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of controls to cases
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -0.9

a = ifelse(gcta$X3 > -0.9, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_0.9.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]


# For MSE threshold of -1

a = ifelse(gcta$X3 > -1, 0, 1)
a = as.data.frame(a)
alf =gcta[,c(1:2)]
alf$new = a[,1]
write.table(alf, file = "pheno_minus_1.txt", col.names = F, row.names = F, quote = F)

# Find the odds/ratio of cases to controls
z = table(a)
z = as.data.frame(z)
z = z[1,]/z[2,]
z[,2]
