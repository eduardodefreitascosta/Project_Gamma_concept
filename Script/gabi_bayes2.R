
if(!require(here)){
  install.packages("here")
}
library(here)


if(!require(BRugs)){
  install.packages("BRugs")
}
library(BRugs)


if(!require(coda)){
  install.packages("coda")
}
library(coda)


if(!require(rjags)){
  install.packages("rjags")
}
library(rjags)





cat(
  'model{
  
  for(i in 1:20){
  for(j in 1:2){
  
  Y[i,j] ~ dpois(D[i])
  }
  D[i]<-N[i]*d[i]*0.1 
  
  N[i]~dgamma(A[i],G[i])
  G[i]<-0.05*d[i]
  A[i]<-5
  }
  
  }' , file={f1<-tempfile()} )




# Defining burn-in, number of simulations, thin and parameters to be shown
burn = 1000;nsim = 1000000;nthin = 250

# Defining the data as a list
dadosjags1 = list(d=(sal1$d), 
                  Y = array(cbind(sal1$placa1, sal1$placa2), dim=c(20,2)))

parms = c("N")

# initializing for adoptation
m1 <- jags.model(f1, dadosjags1,  n.chains=3,n.adapt=1000)


# updating burn-in
update(m1,burn)

# Final sampling
mcmc2 <- coda.samples(m1, parms, n.iter=nsim,thin=nthin)

#plot(mcmc2)
#effectiveSize(mcmc2)
#autocorr.diag(mcmc2)
summary(mcmc2)
#gelman.plot(mcmc2)
