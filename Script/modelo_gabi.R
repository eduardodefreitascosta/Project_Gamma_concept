#########################
#######Modelo Gabi#######
#########################


###Vetores de tempo
tempo_hora<-c(0,4,8,16,21,29,42,54,66,78,88,102,140,188,245,307,379,454,523,643,715,837,941)
delta<-c(0,4,4,8,5,8,13,12,12,12,10,14,38,48,57,62,72,75,69,120,72,122,104)
length(tempo_hora)
length(delta)
###Vetor de Log10 UFC observado

UFC<-c(7.095810988,7.659378483,7.660246366,8.313832357,8.333019113,8.240877032,8.256059775,8.200465326,8.166984212,
       8.143855305,8.112755705,8.287416353,8.256511737,6.967132647,6.877658539,6.460341676,5.886811163,5.173619004,
       4.871754036,4.562146084,4.223370944,3.658713137,3.375805848)
length(UFC)
###Importando dados de ambiente

ph<-c(6.693333,6.700000,6.486667,6.126667,5.955000,5.791667,5.688333,5.551667,5.400000,5.630000,5.058333,5.176667,5.646667,5.646667,
6.005000,5.386667,5.706667,5.510000,6.298333,5.610000,5.780000,5.800000,5.831667)

aw<-(0.957*exp(-tempo_hora*0.000243903660924917))

tem<-c(30,30,30,30,30,30,30,30,30,20,20,22,21,21,20,20,20,20,20,20,21,19,21)


length(ph)
length(aw)
length(tem)

##Parametros do ambiente para crescimento
###Set parametros de crescimento
ph_max<-10.42
ph_min<-3.58
ph_opt<-7
aw_max<-1
aw_min<-0.951
aw_opt<-0.997
t_max<-46.48
t_min<-5.06
t_opt<-39.3


##Gamma e fi pH

gama_ph<-as.numeric()
fi_ph<-as.numeric()
for (i in 1:23){
    if(ph[i]<=ph_min){
       gama_ph[i]<-0
    }else if(ph[i]>ph_max){
        gama_ph[i]<-0
     }else{
        gama_ph[i]<-((ph[i]-ph_max)*(ph[i]-ph_min))/
        ((ph_opt-ph_min)*(ph[i]-ph_opt)-(ph_opt-ph_max)*(ph_min-ph[i]))  
      }
        fi_ph[i]<-(1-gama_ph[i])^3
}

##Gama e fi atividade de agua  
  
gama_aw<-as.numeric()
fi_aw<-as.numeric()
for (i in 1:23){
  
  if(aw[i]<=aw_min){
    gama_aw[i]<-0
  }else if(aw[i]>aw_max){
    gama_aw[i]<-0
  }else{
    gama_aw[i]<-(((aw[i]-aw_max)*(aw[i]-aw_min))/
      ((aw_opt-aw_min)*(aw[i]-aw_opt)-(aw_opt-aw_max)*(aw_min-aw[i]))  )
    
  }
  
  fi_aw[i]<-(1-gama_aw[i])^3
}

##Gamma e fi temperatura

gama_t<-as.numeric()
fi_t<-as.numeric()
for (i in 1:23){
  
  if(tem[i]<=t_min){
    gama_t[i]<-0
  }else if(tem[i]>t_max){
    gama_t[i]<-0
  }else{
    gama_t[i]<-((tem[i]-t_max)*(tem[i]-t_min)^2)/
      ((t_opt-t_min)*((t_opt-t_min)*(tem[i]-t_opt)-(t_opt-t_max)*(t_opt+t_min-2*tem[i])))  
    
  }
  
fi_t[i]<-(1-sqrt(gama_t[i]))^3
}



##Xi zeta e interação

psi<-as.numeric()

for(i in 1:23){
  
  if (fi_aw[i]==1){psi[i]<-1
                  }else  {psi[i]<-(fi_ph[i]/(2*(1-fi_aw[i])*(1-fi_t[i])))+
                            (fi_aw[i]/(2*(1-fi_ph[i])*(1-fi_t[i])))+
                            (fi_t[i]/(2*(1-fi_aw[i])*(1-fi_ph[i])))
                          }
  
}

zeta<-as.numeric()

for (i in 1:23){
    if(psi[i]<=0.5){
      zeta[i]<-1
    }else if(psi[i]>1){
      zeta[i]<-0
    }else {zeta[i]<-2*(1-psi[i])}      
        
    
  
  
}

gama_zeta<-as.numeric()
for(i in 1:23){
gama_zeta[i]<-gama_ph[i]*gama_aw[i]*gama_t[i]*zeta[i]
}



####Parametros do ambiente para inativação
###Set parametros de inativação

Z_ph1<-(-2.3)
Z_ph2<-(-2.38)
ph_ast<-6
ph_opt_reduc<-5.91
Z_aw1<-0.034
Z_aw2<-0.039
aw_opt_reduc<-0.997
Z_t<-24.1
t_ast<-12
t_c<-22
t_opt_reduc<-10.5

##Parametro da weibull fixo
p<-2.27

##Lambdas

lamb_aw_1<-as.numeric()
lamb_aw_2<-as.numeric()
lamb_ph_1<-as.numeric()
lamb_ph_2<-as.numeric()
lamb_t<-as.numeric()

for(i in 1:23){
lamb_ph_1[i]<-((ph[i]-ph_opt_reduc)/Z_ph1)^3-
    ((ph_ast-ph_opt_reduc)/Z_ph1)^3
    
lamb_ph_2[i]<-((ph[i]-ph_opt_reduc)/Z_ph2)^3-
    ((ph_ast-ph_opt_reduc)/Z_ph2)^3



lamb_aw_1[i]<-((aw[i]-aw_opt_reduc)/Z_aw1)^2

lamb_aw_2[i]<-((aw[i]-aw_opt_reduc)/Z_aw2)^2


if (tem[i]<=t_c){
    lamb_t[i]<-2*((t_c-t_opt_reduc)/Z_t^2)*(tem[i]-t_ast)

     }else{lamb_t[i]<-((tem[i]-t_opt_reduc)/Z_t)^2-((t_c-t_opt_reduc)*
                        (2*t_ast-t_c-t_opt_reduc)/Z_t^2)

}

}

#Modelagem

CFU<-as.numeric()
CFU[1]<-UFC[1]
beta1<-numeric()
beta2<-numeric()
modelo<-function(param){
  
for (i in 2:23){
  
  beta1[i]<-10^(param[3]-(lamb_aw_1[i]+lamb_ph_1[i]+lamb_t[i]) )
  beta2[i]<-10^(param[4]-(lamb_aw_2[i]+lamb_ph_2[i]+lamb_t[i]) ) 
  
  
          if(gama_zeta[i]>0){CFU[i]<-CFU[i-1]+log10(exp(param[1]*gama_zeta[i]*delta[i]) )
                            }else 
                              {CFU[i]<- CFU[i-1]-log10(1+10^param[2])  +  
                      
                      log10(   10^-( ((delta[i]/beta1[i])^p )+param[2] )    + 
                               
                               10^-( (delta[i]/beta2[i])^p )   )    
    
                              }

               }   

 return( sum((UFC-CFU)^2))
}

#Otimizacao
inici=c(3.6,0.2,6,39)
fit=try(nlm(modelo,inici,hessian=TRUE))

RMSE1=sqrt(sum((UFC-CFU)^2)/23)
##Plot

for (i in 2:23){
  
  beta1[i]<-(fit$estimate[3]-(lamb_aw_1[i]+lamb_ph_1[i]+lamb_t[i]) )
  beta2[i]<-(fit$estimate[4]-(lamb_aw_2[i]+lamb_ph_2[i]+lamb_t[i]) ) 
  
  
  if(gama_zeta[i]>0){CFU[i]<-CFU[i-1]+log10(exp(fit$estimate[1]*gama_zeta[i]*delta[i]) )
  }else 
  {CFU[i]<- CFU[i-1]-log10(1+10^fit$estimate[2])  +  
    
    log10(   10^-( ((delta[i]/10^beta1[i])^p) +fit$estimate[2] )    + 
               
               10^-( (delta[i]/10^beta2[i])^p )   )    
  
  }
  
}
par(mfrow=c(1,2))
par(xpd=FALSE)
par(mar = c(5, 5, 1, 1))
plot(tempo_hora,UFC,xlim=c(0,1000),ylim=c(3,9), ylab='log(CFU/g)',xlab='Time (hours)',col = "black", 
     pch = 20,lwd=2,cex.lab=2,cex.axis=1.5,cex=2)
par(xpd=FALSE)
lines(tempo_hora,CFU,type='l',col='black',pch = 20,lwd=2,cex=2)
#text(800,8,"A",cex=3,lwd=2)



reg<-lm(UFC~CFU)
par(xpd=F)
par(mar = c(5, 5, 1, 1))
plot(UFC,CFU,ylab='Observed log(CFU/g)',xlab='Predicted log(CFU/g)',col='black',pch=20,lwd=2,cex=2,cex.axis=1.5,cex.lab=2)
abline(reg, col="black",lwd=2)
mylabel = bquote(italic(R)^2 == .(format('0.967', digits = 3)))
text(x = 4.5, y = 8, labels = mylabel)
text(4.5,7.5,c('y= -0.34+1.05'))



##Intervalos de confian?a

error<-sqrt(diag(solve(fit$hessian)))


Intervalos<-rbind(
c(fit$estimate[1],fit$estimate[1]-1.96*error[1],fit$estimate[1]+1.96*error[1]),
c(fit$estimate[2],fit$estimate[2]-1.96*error[2],fit$estimate[2]+1.96*error[2]),
c(fit$estimate[3],fit$estimate[3]-1.96*error[3],fit$estimate[3]+1.96*error[3]),
c(fit$estimate[4],fit$estimate[4]-1.96*error[4],fit$estimate[4]+1.96*error[4])
)

rownames(Intervalos)<-c( expression((mu[opt])),
                        expression(alpha),
                        expression((beta[1])),
                        expression((beta[2])))

colnames(Intervalos)<-c('Estimated','Min', 'Max')

library(knitr)
kable(Intervalos,digits = 3,align = 'c')





############################
# Analise de sensibilidade #
############################
tempo_hora2<-c(0,4,8,14,20,26,33,40.5,64.5,89.5,161.5,233.5,305.5,353.5,425.5,521.5,689.5,809.5,1001.5,1121.5)
delta2<-numeric()
delta2[1]<-0
for(i in 2: length(tempo_hora2))
delta2[i]<-tempo_hora2[i]-tempo_hora2[i-1]

UFC2<-c(6.836765909,6.417383185,6.415412014,6.604785063,7.056801123,7.379198572,7.265907878,7.208778583,
7.019579637,6.992571499,6.682268214,6.199567189,5.282010673,5.393627883,5.139824895,4.787727592,3.441734324,
3.363568304,2.264423854,1.821826812)

ph2<-c(6.50600, 6.41000,6.42400,6.31400, 6.10000, 5.75400, 5.65400,5.29400,6.14400,
       5.15200, 5.09200,5.03600,5.32200,5.68500,5.38167,5.59000,5.38667,5.80667,
       5.77667,6.06333)


#aw2<-c(0.955000,0.954208,0.953416,0.952230,0.951046,0.949863,0.948485,0.947010,
 #     0.942307,0.937433,0.923535,0.909844,0.896356, 0.887475,0.874318,
#      0.857078, 0.827723,  0.807373,0.761835, 0.796732)
      
aw2<-0.955*exp(-tempo_hora2*0.000203620182165236)


tem2<-c(25,25,25,25,25,25,25,25,17,17,17,17,17,17,17,17,17,17,17,17)

length(tempo_hora2)
length(UFC2)
length(ph2)
length(aw2)
length(tem2)
length(delta2)

##Parametros do ambiente para crescimento
###Set parametros de crescimento
ph_max<-10.42
ph_min<-3.58
ph_opt<-7
aw_max<-1
aw_min<-0.951
aw_opt<-0.997
t_max<-46.48
t_min<-5.06
t_opt<-39.3


##Gamma e fi pH

gama_ph2<-as.numeric()
fi_ph2<-as.numeric()
for (i in 1:20){
  if(ph2[i]<=ph_min){
    gama_ph2[i]<-0
  }else if(ph2[i]>ph_max){
    gama_ph2[i]<-0
  }else{
    gama_ph2[i]<-((ph2[i]-ph_max)*(ph2[i]-ph_min))/
      ((ph_opt-ph_min)*(ph2[i]-ph_opt)-(ph_opt-ph_max)*(ph_min-ph2[i]))  
  }
  fi_ph2[i]<-(1-gama_ph2[i])^3
}

##Gama e fi atividade de agua  

gama_aw2<-as.numeric()
fi_aw2<-as.numeric()
for (i in 1:20){
  
  if(aw2[i]<=aw_min){
    gama_aw2[i]<-0
  }else if(aw2[i]>aw_max){
    gama_aw2[i]<-0
  }else{
    gama_aw2[i]<-(((aw2[i]-aw_max)*(aw2[i]-aw_min))/
                   ((aw_opt-aw_min)*(aw2[i]-aw_opt)-(aw_opt-aw_max)*(aw_min-aw2[i]))  )
    
  }
  
  fi_aw2[i]<-(1-gama_aw2[i])^3
}

##Gamma e fi temperatura

gama_t2<-as.numeric()
fi_t2<-as.numeric()
for (i in 1:20){
  
  if(tem2[i]<=t_min){
    gama_t2[i]<-0
  }else if(tem2[i]>t_max){
    gama_t2[i]<-0
  }else{
    gama_t2[i]<-((tem2[i]-t_max)*(tem2[i]-t_min)^2)/
      ((t_opt-t_min)*((t_opt-t_min)*(tem2[i]-t_opt)-(t_opt-t_max)*(t_opt+t_min-2*tem2[i])))  
    
  }
  
  fi_t2[i]<-(1-sqrt(gama_t2[i]))^3
}



##Xi zeta e interação

psi2<-as.numeric()

for(i in 1:20){
  
  if (fi_aw2[i]==1){psi2[i]<-1
  }else  {psi2[i]<-(fi_ph2[i]/(2*(1-fi_aw2[i])*(1-fi_t2[i])))+
    (fi_aw2[i]/(2*(1-fi_ph2[i])*(1-fi_t2[i])))+
    (fi_t2[i]/(2*(1-fi_aw2[i])*(1-fi_ph2[i])))
  }
  
}

zeta2<-as.numeric()

for (i in 1:20){
  if(psi2[i]<=0.5){
    zeta2[i]<-1
  }else if(psi2[i]>1){
    zeta2[i]<-0
  }else {zeta2[i]<-2*(1-psi2[i])}      
  
  
  
  
}

gama_zeta2<-as.numeric()
for(i in 1:20){
  gama_zeta2[i]<-gama_ph2[i]*gama_aw2[i]*gama_t2[i]*zeta2[i]
}



####Parametros do ambiente para inativação
###Set parametros de inativação

Z_ph1<-(-2.3)
Z_ph2<-(-2.38)
ph_ast<-6
ph_opt_reduc<-5.91
Z_aw1<-0.034
Z_aw2<-0.039
aw_opt_reduc<-0.997
Z_t<-24.1
t_ast<-12
t_c<-22
t_opt_reduc<-10.5

##Parametro da weibull fixo
p<-2.27

##Lambdas

lamb_aw_12<-as.numeric()
lamb_aw_22<-as.numeric()
lamb_ph_12<-as.numeric()
lamb_ph_22<-as.numeric()
lamb_t2<-as.numeric()

for(i in 1:20){
  lamb_ph_12[i]<-((ph2[i]-ph_opt_reduc)/Z_ph1)^3-
    ((ph_ast-ph_opt_reduc)/Z_ph1)^3
  
  lamb_ph_22[i]<-((ph2[i]-ph_opt_reduc)/Z_ph2)^3-
    ((ph_ast-ph_opt_reduc)/Z_ph2)^3
  
  
  
  lamb_aw_12[i]<-((aw2[i]-aw_opt_reduc)/Z_aw1)^2
  
  lamb_aw_22[i]<-((aw2[i]-aw_opt_reduc)/Z_aw2)^2
  
  
  if (tem2[i]<=t_c){
    lamb_t2[i]<-2*((t_c-t_opt_reduc)/Z_t^2)*(tem2[i]-t_ast)
    
  }else{lamb_t2[i]<-((tem2[i]-t_opt_reduc)/Z_t)^2-((t_c-t_opt_reduc)*
                                                   (2*t_ast-t_c-t_opt_reduc)/Z_t^2)
  
  }
  
}

#Modelagem

##Plot

CFU2<-as.numeric()
CFU2[1]<-UFC2[1]
beta12<-numeric()
beta22<-numeric()
for (i in 2:20){
  
  beta12[i]<-10^(fit$estimate[3]-(lamb_aw_12[i]+lamb_ph_12[i]+lamb_t2[i]) )
  beta22[i]<-10^(fit$estimate[4]-(lamb_aw_22[i]+lamb_ph_22[i]+lamb_t2[i]) ) 
  
  
  if(gama_zeta2[i]>0){CFU2[i]<-CFU2[i-1]+log10(exp(fit$estimate[1]*gama_zeta2[i]*delta2[i]) )
  }else 
  {CFU2[i]<- CFU2[i-1]-log10(1+10^fit$estimate[2])  +  
    
    log10(   10^-( ((delta2[i]/beta12[i])^p )+fit$estimate[2] )    + 
               
               10^-( (delta2[i]/beta22[i])^p )   )    
  
  }
  
}

RMSE2=sqrt(sum((UFC2-CFU2)^2)/20)
par(mfrow=c(1,2))
par(mar=c(5,5,1,1))
plot(tempo_hora2,UFC2,xlim=c(0,1200),ylim=c(1,7.5), ylab='log(CFU/g)',xlab='Time (hours)',col = "black", 
     pch = 20,cex=2,cex.axis=1.5,cex.lab=2)
lines(tempo_hora2,CFU2,type='l',col='black',pch = 20,lwd=2)
#text(1000,6,"B",lwd=2,cex=3)


reg2<-lm(UFC2~CFU2)
par(mar=c(5,5,1,1))
plot(UFC2,CFU2,ylab='Observed log(CFU/g)',xlab='Predicted log(CFU/g)',col='black',
     pch = 20,cex=2,cex.axis=1.5,cex.lab=2)
par(xpd=F)
abline(reg2, col="black",lwd=2)
mylabel = bquote(italic(R)^2 == .(format('0.9176', digits = 3)))
text(x = 3.5, y = 7, labels = mylabel)
text(3.5,6.5,c('y= -0.4+1.073'))




##Extra plots
par(xpd=TRUE)
corners = par("usr")
par(mar = c(4, 5, 3, 7))
plot(tempo_hora,zeta,type='l',xlim=c(0,1000),ylim=c(0,1),ylab=expression(paste(gamma[x],(x[i]))),xlab='Time (hours)',lwd=2,cex.axis=1.5,cex.lab=2)
lines(tempo_hora,gama_aw,lty=3,lwd=2)
lines(tempo_hora,gama_ph,lty=4,lwd=2)
lines(tempo_hora,gama_t,lty=5,lwd=2)
legend(1070,1,inset=c(-0.2,0), c(expression(xi),'Aw','pH','T'), cex=1.5, lty=c(1,3,4,5), title=" ",bty = 'n',lwd=2)
#text(x= corners[3]-200, y = (corners[4])+0.085,'A',cex = 2.5)

par(xpd=TRUE)
corners = par("usr")
par(mar = c(4, 5, 3, 7))
plot(tempo_hora2,zeta2,type='l',xlim=c(0,1200),ylim=c(0,1),ylab=expression(paste(gamma[x],(x[i]))),xlab='Time (hours)',lwd=2,cex.axis=1.5,cex.lab=2)
lines(tempo_hora2,gama_aw2,lty=3,lwd=2)
lines(tempo_hora2,gama_ph2,lty=4,lwd=2)
lines(tempo_hora2,gama_t2,lty=5,lwd=2)
legend(1270,1,inset=c(-0.2,0), c(expression(xi),'Aw','pH','T'), cex=1.5, lty=c(1,3,4,5), title=" ",bty = 'n',lwd=2)
#text(x= corners[3]-200, y = (corners[4])+0.085,'B',cex=2.5)





######################
##Modelo estocastico##
######################



niter=10000
CFU3<-matrix(nrow=23,ncol=niter)
for (i in 2:23){
  for (j in 1:niter){
  
  CFU3[1,j]<-sample(c(-1.187086643,-0.259637311,0.740362689,1.740362689,2.447158031),1,
                    prob = c(0.3607,0.5082,0.0328,0.0656,0.0328))
  beta1[i]<-(fit$estimate[3]-(lamb_aw_1[i]+lamb_ph_1[i]+lamb_t[i]) )
  beta2[i]<-(fit$estimate[4]-(lamb_aw_2[i]+lamb_ph_2[i]+lamb_t[i]) ) 
  
  
  if(gama_zeta[i]>0){CFU3[i,j]<-CFU3[i-1,j]+log10(exp(fit$estimate[1]*gama_zeta[i]*delta[i]) )
  }else 
  {CFU3[i,j]<- CFU3[i-1,j]-log10(1+10^fit$estimate[2])  +  
    
    log10(   10^-( ((delta[i]/10^beta1[i])^p) +fit$estimate[2] )    + 
               
               10^-( (delta[i]/10^beta2[i])^p )   )    
  
    }
  
  }
}


media<-apply(CFU3, 1, mean)

perc2_5<-numeric()
perc97_5<-numeric()
for(i in 1:23){
perc2_5[i]<-quantile(CFU3[i,],0.025)
perc97_5[i]<-quantile(CFU3[i,],0.975)
}

par(xpd=NA)
par(mar = c(5, 5, 1, 1))
plot(tempo_hora,media,xlim=c(0,1000),ylim=c(-5,4), ylab='log(CFU/g)',xlab='Time (hours)',col = "black", 
     pch = 20,lwd=2,cex.lab=2,cex.axis=1.5,type='l')
lines(tempo_hora,perc2_5,pch = 20,lwd=2,lty=2)
lines(tempo_hora,perc97_5,pch = 20,lwd=2,lty=2)
#legend(1070,4,inset=c(-0.2,0), c('2.5%','Mean','97.5%'), cex=0.8, lty=c(2,1,2), title=" ",bty = 'n',lwd=2,xpd=TRUE)

abline(h=0,xpd=FALSE)


