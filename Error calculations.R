# Simulating Type 1 error rate

n=10000 # testing 10,000 times
t1err=0
for (i in 1:n){
  x=rnorm(1000, 0, 1)
  if (((t.test(x, mu=0))$p.value)<=0.05) (t1err=t1err+1) 
}
cat("Type I error rate in percentage is", (t1err/n)*100,"%")


# Simulating Type 2 error rate
n=10000 # testing 10,000 times
t2err=0
for (i in 1:n){
  x=rnorm(100, 2, 1)
  if (((t.test(x, mu=0))$p.value)>0.05) (t2err=t2err+1) 
}
cat("Type II error rate in percentage is", (t2err/n)*100,"%")





# Anecdotal tests
# type I error
x=rnorm(10^3)
sum(x>1.68)/10^3
# type II error
y=rnorm(10^3,2)
sum(y<1.68)/10^3





