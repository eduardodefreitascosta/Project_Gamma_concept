
library(multcompView)
library(lsmeans)
library(lme4)
library(car)

head(gabi)
gabi$time<-as.factor(gabi$time)

summary(lmer(mean~tecnica+(1|time/par),data=gabi)->a)
lsmeans(a, pairwise ~ tecnica,adjust="tukey")
Anova(a)
