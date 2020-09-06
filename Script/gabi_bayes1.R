

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
  
  for(i in 1:23){
  for(j in 1:2){
  for (k in 1:2){
  Y[i,j,k] ~ dpois(D[i,j])
  }
  D[i,j]<-N[i]*d[i,j]*0.1 
  
  }
  N[i]~dgamma(A[i],G[i])
  G[i]<-0.05*d[i,1]
  A[i]<-5
  }
  
  }' , file={f1<-tempfile()} )



# Defining burn-in, number of simulations, thin and parameters to be shown
burn = 1000;nsim = 1000000;nthin = 250

# Defining the data as a list
dadosjags = list(d=array(cbind(sal$d[1:23],sal$d[24:46]),dim=c(23,2)), 
                 Y = array(cbind(sal$placa1, sal$placa2), dim=c(23,2,2)))

parms = c("N")

# initializing for adoptation
m1 <- jags.model(f1, dadosjags,  n.chains=3,n.adapt=1000)


# updating burn-in
update(m1,burn)

# Final sampling
mcmc1 <- coda.samples(m1, parms, n.iter=nsim,thin=nthin)


#plot(mcmc1)
#effectiveSize(mcmc1)
#autocorr.diag(mcmc1)
#summary(mcmc1)
#gelman.plot(mcmc1)


