# Advanced genetic drift simulator

GenDriftSim = function(Pop = Pop, Gen = Gen, NM, NF, P, graph = "y", histo = "y"){
  P = (2*(NM+NF))*P
  NE = round((4*NM*NF)/(NM+NF),0)
  SR = round(NM/NF,2)
  Na = NM+NF
  if(graph=="y"){
    plot(c(0,0),type = "n", main = bquote('N'[M]~'/ N'[F]~'='~.(SR)*', N'[A]~'='~.(Na)*', N'[E]~'='~.(NE)), cex.main = 1, 
         xlim = c(1,Gen), ylim=c(0,1), xlab = "Generations", ylab = "Fequency of focal allele")
  }else{}
  for (i in 1:Pop){
    N = NM+NF
    startA = as.vector(c(rep(1, times = P),rep(0, times = (2*N)-P)))
    Population = matrix(c(
      c(sample(startA, size = 2*N, replace = FALSE)),
      c(rep("M", times = NM), rep("F", times = NF))),
      ncol = 3)
    SimResults[(Gen*i)+1-Gen, 3] <<- sum(as.numeric(Population[,1:2]))/(N*2)
    for(j in 1:(Gen-1)){
      
      Population = matrix(c(
        c(sample(sample(Population[(1:NM),1:2], replace = TRUE),N, replace = TRUE)),
        c(sample(sample(Population[(1+NM):N,1:2], replace = TRUE),N, replace = TRUE)),
        c(rep("M", times = NM), rep("F", times = NF))), ncol = 3)
      SimResults[(Gen*i)+1+j-Gen, 3] <<- sum(as.numeric(Population[,1:2]))/(N*2)
    }
    s = (i*Gen)-Gen+1; e = i*Gen
    r = as.vector(SimResults[s:e, 3])
    if(graph=="y"){
      points(r~c(1:Gen), type = "l")
    }else{}
  }
  if(histo == "y"){SimResults[,1] = rep(1:Pop, each = Gen)
  SimResults[,2] = rep(1:Gen, times = Pop)
  hist(SimResults[,3][SimResults[,2]==Gen], breaks = 100, cex.lab = 0.7, cex.axis = 0.7, xlim = c(0,1), cex.main = 1, main = bquote('N'[M]~'/ N'[F]~'='~.(SR)*', N'[A]~'='~.(Na)*', N'[E]~'='~.(NE)), xlab = paste0("Frequency of focal allele after ",Gen," Generations"))
  }else{}
}

Pop = 10
#Gen = c(10,20,50,100,200,400, 800, 1000, 2000)
#for (z in Gen){
#  SimResults = matrix(data = NA, ncol = 3, nrow = Gen*Pop)
#  GenDriftSim(Pop = Pop, Gen = Gen, NM = 100, NF = 100, P = P, graph = "y",  histo = "n")
#}
Gen = 100
P = 0.5

SimResults = matrix(data = NA, ncol = 3, nrow = Gen*Pop)
par(mfrow=c(3,3))
GenDriftSim(Pop = Pop, Gen = Gen, NM = 100, NF = 100, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 150, NF = 150, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 200, NF = 900, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 300, NF = 180, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 350, NF = 900, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 400, NF = 380, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 450, NF = 1900, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 500, NF = 2180, P = P, graph = "y",  histo = "n")
GenDriftSim(Pop = Pop, Gen = Gen, NM = 550, NF = 180, P = P, graph = "y",  histo = "n")
#dev.off()