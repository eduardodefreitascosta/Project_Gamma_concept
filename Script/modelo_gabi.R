#########################
#######General model#####
#########################

###Time vector
tempo_hora<-c(0,4,8,16,21,29,42,54,66,78,88,102,140,188,245,307,379,454,523,643,715,837,941)
delta<-c(0,4,4,8,5,8,13,12,12,12,10,14,38,48,57,62,72,75,69,120,72,122,104)
length(tempo_hora)
length(delta)
###Vetor of Log10 UFC observed


UFC<-as.numeric(log10(summary(mcmc1)$quantiles)[,3])


length(UFC)
###Environmental data

ph<-c(6.693333,6.700000,6.486667,6.126667,5.955000,5.791667,5.688333,5.551667,5.400000,5.630000,5.058333,5.176667,5.646667,5.646667,
6.005000,5.386667,5.706667,5.510000,6.298333,5.610000,5.780000,5.800000,5.831667)

aw<-(0.957*exp(-tempo_hora*0.000243903660924917))


#time_T<-c(0,24,48,72,96,120,144,168,192,216,240,
#264,288,312,336,360,384,408,432,456,480,504,528,
#552,576,600,624,648,672,696,720,744,768,792,816,
#840,864,888,912,936)

tem<-c(30,30,30,30,30,30,30,30,30,20.12,21,21.5,20.5,20.8,20,20,19.8,19.5,19.2,20.2,20.2,20,20.5)
         
         
#temp<-approxfun(time_T,tem_1)

#tem<-temp(tempo_hora)

#tem[23]<-20.5

#tem<-c(30,30,30,30,30,30,30,30,30,20,20,22,21,21,20,20,20,20,20,20,21,19,21)

length(ph)
length(aw)
length(tem)

##Environmental parameters for growth 
###Set parameters for growth
ph_max<-10.42
ph_min<-3.58
ph_opt<-7
aw_max<-1
aw_min<-0.951
aw_opt<-0.997
t_max<-46.48
t_min<-5.06
t_opt<-39.3


##Gamma and fi pH

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

##Gama and fi Aw 
  
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

##Gamma and fi temperature

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



##Xi zeta and interaction

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



####Environment parameters for inactivation
###Set parameters for inactivation

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

#Modeling

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

#Optimization
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



theme_set(theme_minimal())
p1<-ggplot()+
  geom_point(aes(x=tempo_hora,y=UFC),size=2)+
  geom_line(aes(x=tempo_hora,y=CFU),size=1)+
  labs(title= "A",
       y='log(CFU/g)', x = 'Time(hours)')+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  theme(text = element_text(size = 20))+
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.3, "cm"))


summary(lm(UFC~CFU)->reg)

p2<-ggplot()+
  geom_point(aes(x=UFC,y=CFU),size=2)+
  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2],size=1)+
  labs(title= "B",
       y='Observed log(CFU/g)', x = 'Predicted log(CFU/g)')+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  theme(text = element_text(size = 20))+
  annotate(geom="text", x=5, y=8, label=paste("y=", round(reg$coefficients[1],2),"+",round(reg$coefficients[2],2),"x",sep=""),
           color="black",size=5)+
  annotate(geom="text",x=5,y=7.5,label=bquote(italic(R)^2 == .(format(summary(reg)$r.squared, digits = 3))),size=5)+
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))


tiff(file=here("Figura","modelo.tiff"), height = 4, width = 7, units = 'in', res=300)

grid.arrange(p1, p2, widths = c(0.5, 0.5),layout_matrix = rbind(c(1, 2),
                                                                c(1, 2)))

dev.off()

##Confidence intervals
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
# Validation model         #
############################
tempo_hora2<-c(0,4,8,14,20,26,33,40.5,64.5,89.5,161.5,233.5,305.5,353.5,425.5,521.5,689.5,809.5,1001.5,1121.5)
delta2<-numeric()
delta2[1]<-0
for(i in 2: length(tempo_hora2))
delta2[i]<-tempo_hora2[i]-tempo_hora2[i-1]


UFC2<-as.numeric(log10(summary(mcmc2)$quantiles)[,3])
  


ph2<-c(6.50600, 6.41000,6.42400,6.31400, 6.10000, 5.75400, 5.65400,5.29400,6.14400,
       5.15200, 5.09200,5.03600,5.32200,5.68500,5.38167,5.59000,5.38667,5.80667,
       5.77667,6.06333)

aw2<-0.955*exp(-tempo_hora2*0.000203620182165236)



#time_T2<-c(0,24,144,168,192,216,240,264,288,312,336,360,384,408,432,456,480,504,528,
#          552,576,600,624,648,672,720,744,768,792,816,888,912,936,960,984,1008,1056,1080,1104)

tem2<-c(25,25,25,25,25,25,25,25,25,17.23,17,16,17,17.5,16.17,17,17,17.5,17.8,17)
        

#temp22<-approxfun(time_T2,tem_12)

#tem2<-temp22(tempo_hora2)

#tem2[c(1:5)]<-25
#tem2[20]<-17


#tem1.2<-c(25,25,25,25,25,25,25,25,17,17,17,17,17,17,17,17,17,17,17,17)


length(tempo_hora2)
length(UFC2)
length(ph2)
length(aw2)
length(tem2)
length(delta2)

##Environmental parameters for growth
###Set parameters for growth
ph_max<-10.42
ph_min<-3.58
ph_opt<-7
aw_max<-1
aw_min<-0.951
aw_opt<-0.997
t_max<-46.48
t_min<-5.06
t_opt<-39.3


##Gamma and fi pH

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

##Gama and fi Aw

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

##Gamma and fi temperature

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



##Xi zeta and interação

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



####Environment parameters for inactivation
###Set parameters for inactivation

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

##Fix Weibull parameters
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

#Modeling

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



theme_set(theme_minimal())
p3<-ggplot()+
  geom_point(aes(x=tempo_hora2,y=UFC2),size=2)+
  geom_line(aes(x=tempo_hora2,y=CFU2),size=1)+
  labs(title= "A",
       y='log(CFU/g)', x = 'Time(hours)')+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  theme(text = element_text(size = 20))+
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.3, "cm"))



summary(lm(UFC2~CFU2)->reg2)

p4<-ggplot()+
  geom_point(aes(x=UFC2,y=CFU2),size=2)+
  geom_abline(intercept = reg2$coefficients[1], slope = reg2$coefficients[2],size=1)+
  labs(title= "B",
       y='Observed log(CFU/g)', x = 'Predicted log(CFU/g)')+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  theme(text = element_text(size = 20))+
  annotate(geom="text", x=4, y=7, label=paste("y=", round(reg2$coefficients[1],2),"+",round(reg2$coefficients[2],2),"x",sep=""),
           color="black",size=5)+
  annotate(geom="text",x=4,y=6.5,label=bquote(italic(R)^2 == .(format(summary(reg2)$r.squared, digits = 3))),size=5)+
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))


tiff(file=here("Figura","validacao.tiff"), height = 4, width = 7, units = 'in', res=300)

grid.arrange(p3, p4, widths = c(0.5, 0.5),layout_matrix = rbind(c(1, 2),
                                                                c(1, 2)))

dev.off()



##Extra plots

env1<-cbind.data.frame( rep(tempo_hora,4),c(zeta,gama_aw,gama_ph,gama_t),
                        c(rep("z",23),rep("aw",23),rep("ph",23),rep("t",23)))

colnames(env1)<-c("time","value","factor")

tiff(file=here("Figura","environment1.tiff"), height = 4, width = 6, units = 'in', res=300)

ggplot(env1,aes(x=time,y=value,group=factor))+
  geom_line(aes(linetype=factor),size=1.5)+
    labs( y=expression(paste(gamma[x],(x[i]))), x = 'Time (hours)')+
    
  scale_linetype_manual(values = c("solid","dashed", "dotted","dotdash"),
                        labels = c('Aw','pH','T',expression(xi)),
                        name="")+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5),guide=F)+
  theme(legend.key.width = unit(2,"cm"))+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  
  theme(text = element_text(size = 20))+
  theme(plot.margin = margin(0.2, 0.3, 0.3, 0.3, "cm"))

dev.off()


env2<-cbind.data.frame( rep(tempo_hora2,4),c(zeta2,gama_aw2,gama_ph2,gama_t2),
                        c(rep("z",20),rep("aw",20),rep("ph",20),rep("t",20)))

colnames(env2)<-c("time","value","factor")


tiff(file=here("Figura","environment2.tiff"), height = 4, width = 6, units = 'in', res=300)

ggplot(env2,aes(x=time,y=value,group=factor))+
  geom_line(aes(linetype=factor),size=1.5)+
  labs( y=expression(paste(gamma[x],(x[i]))), x = 'Time (hours)')+
  
  scale_linetype_manual(values = c("solid","dashed", "dotted","dotdash"),
                        labels = c('Aw','pH','T',expression(xi)),
                        name="")+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5),guide=F)+
  theme(legend.key.width = unit(2,"cm"))+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  
  theme(text = element_text(size = 20))+
  theme(plot.margin = margin(0.2, 0.3, 0.3, 0.3, "cm"))

dev.off()




######################
##Statistical model ##
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

stoc<-cbind.data.frame(rep(tempo_hora,3),c(media,perc2_5,perc97_5),
                       c(rep("mean",23),rep("lower",23),rep("upper",23)))


colnames(stoc)<-c("time","value","stat")


tiff(file=here("Figura","estocastico.tiff"), height = 4, width = 6, units = 'in', res=300)

ggplot(stoc,aes(x=time,y=value,group=stat))+
  geom_line(aes(linetype=stat),size=1.5)+
  labs( y="log(CFU/g)", x = 'Time (hours)')+
  
  scale_linetype_manual(values = c("dashed","solid", "dashed"),
                        name="")+
  scale_size_manual(values=c(1.5,1.5,1.5),guide=F)+
  theme(legend.position="none")+
  
  theme(axis.text.x = element_text(colour = 'black', size = 20),
        axis.title.x=element_text(size = 20, 
                                  hjust = 0.5, vjust = 0.2),axis.title = element_text(size = 20)) +
  
  theme(axis.text.y = element_text(colour = 'black', size = 20), 
        axis.title.y = element_text(size = 20, 
                                    hjust = 0.5, vjust = 3),axis.title = element_text(size = 20)) +
  
  theme(text = element_text(size = 20))+
  theme(plot.margin = margin(0.2, 0.3, 0.3, 0.3, "cm"))+
  
  annotate(geom="text", x=185, y=0.2, label="132 hours",
           color="black")+
  
  annotate(geom="text", x=715, y=0.2, label="670 hours",
           color="black")+
  
  geom_segment(aes(x = 132, y = -5, xend = 132, yend = 0))+
  
  geom_segment(aes(x = 670, y = -5, xend = 670, yend = 0))+
  
  geom_segment(aes(x = 0, y = 0, xend = 670, yend = 0))
  

dev.off()


#Save the final model
save.image(file=here("Resultado" ,"modelo.RData"))


###Observed values in the Bayesian analysis

write.csv(
cbind(
log10(summary(mcmc1)$statistics[,1]),
log10(summary(mcmc1)$quantiles[,c(1,3,5)])
),here("Resultado","bayes_out1.csv")

)

write.csv(
cbind(
  log10(summary(mcmc2)$statistics[,1]),
  log10(summary(mcmc2)$quantiles[,c(1,3,5)])
),here("Resultado","bayes_out2.csv")

)
