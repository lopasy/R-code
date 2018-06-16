# Simulating genetic drift

plot(1, type = "n", xlab = "generations", ylab = "frequency of allele A", xlim = range(1:50), ylim = 0:1,
     main = "Allele frequency change due to genetic drift")
set.seed(1)
{for (xxx in 1:10){
  print(xxx)
  
  N = 50
  pA = c()
  pA[1] = 0.5
  i = 1
  
  while((pA[i] < 1) & (pA[i] > 0)){
    nA = 0
    
    for (j in 1:N){
      random = runif(1)
      if(random < pA[i]){nA = nA + 1}
    }
    
    pA[i+1] = nA/N
    i = i+1
    
  }
  lines(1:j, pA[1:j])}
}


plot(1, type = "n", xlab = "generations", ylab = "frequency of allele A", xlim = range(1:50), ylim = 0:1,
     main = "Allele frequency change due to genetic drift")
set.seed(1)
{for (xxx in 1:10){
  print(xxx)
  
  N = 50
  pA = c()
  pA[1] = 0.5
  i = 1
  
  while((pA[i] < 1) & (pA[i] > 0)){
    nA = 0
    
    for (j in 1:N){
      random = runif(1)
      if(random < pA[i]){nA = nA + 1}
    }
    
    pA[i+1] = nA/N
    i = i+1
    
  }
  lines(1:j, pA[1:j])}
}



# Estimating time to fixation.
N = c(20, 50, 100, 200, 400, 800)
reps = 100
matrix = matrix(0, nrow=reps, ncol = length(N))
  for(t in 1:100){
    Fixation = c()
      for(r in N){
        pA = c()
        pA[1] = 0.5
        i = 1
      
           while((pA[i] < 1) & (pA[i] > 0)){
             nA = rbinom(n = 1, size = r, prob = pA[i])
             pA[i+1] = nA/r
             i = i+1}
      
        Fixation = c(Fixation, i)
    }
   matrix[t,] = Fixation
}

Means = colMeans(matrix)
plot(N, Means, xlab = "Population Size", ylab = "Average time to fixation")  
abline(lm(Means~N))
  
  
  
  
  