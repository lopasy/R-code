fractal = function(start, reps, expo){
n = matrix(nrow = 1, ncol = reps)
i = 1
n[,i] = ((start*2/4 + 5) / sin(39))/exp(expo)
  for (i in seq_along(2:reps)){
      if(i < reps) n[,i+1] = ((n[,i]*2/4 + 5) / sin(39))/exp(expo)
      else {
        break
      }
      
  }
return(n) 
}
{n = fractal(143,20, 5)
n1 = fractal(20,20, 1)
n2 = fractal(143,2000, 0)
#plot(density(n))

plot(1, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,5))
points(n,n1, type = "l")
points(n1, type = "l")
points(n2,n2, type = "l")}

hist(n)



