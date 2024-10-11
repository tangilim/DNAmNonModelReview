#Epigenetic clocks in non-model species meta-analysis
#https://onlinelibrary.wiley.com/doi/10.1111/mec.17065
#Final model and graphs
#Marianthi Tangili

setwd("C:/Users/MWP-/Desktop/PhD/Review_DNAmWild/Figures/FINALfigures")


#load required packages
library(dbplyr)
library(Matrix)
library(tidyverse)
library(brms)
library(readxl)
library(tidybayes)
library(ggdist)
library(forestplot)
library(bayesplot)
library(sjstats)
library(ggplot2)
library(emmeans)
library(ggpubr)
library(GGMnonreg)
library(cowplot)
library(gtable)
library(gridExtra)

#simple graph theme
graph_theme<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour =NA ,fill=NA,size=0.5),
  axis.title.x=element_text(size=15,hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=15,hjust=0.5,vjust=0.5,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=12),
  axis.text.y=element_text(colour="black",angle=0,size=12),
  axis.ticks=element_line(colour="black",size=0.8),
  axis.line.x=element_line(size=0.8),
  axis.line.y=element_line(size=0.8))

## Import the data
d <- read_xlsx("C:/Users/MWP-/Desktop/SumRevOct2022.xlsx",na="NA")

## Turn some variables into factors and change order of levels for some
d <- d %>% mutate(WildCap = factor(WildCap),
                  WildCap = fct_relevel(WildCap, "Captive", "Wild", "Both"),
                  Method = factor(Method),
                  Method = fct_relevel(Method,"Microarray","Targeted BS", "RRBS"),
                  Tissue = factor(Tissue),
                  Tissue = fct_relevel(Tissue,"Blood","Multiple", "Other"),
                  TissueP = as.factor(Prol),
                  TissueP= fct_relevel(Prol, "Proliferative", "Non-proliferative", "Multiple"),
                  Class = factor(Class),
                  Class = fct_relevel(Class,"Mammal","Fish","Bird"),
                  Species = factor(Species),
                  Study = factor(1:nrow(d)))

## We can't use data with missing Rsq or sample size.
## Sample size is required to compute standard errors for correlations
d <- d %>% filter(!is.na(Rsq),!is.na(SSInd))


## Compute Fisher's z-transformation (z=(1/2)log((1+r)/(1-r))) for correlations 
## and its standard error 1/sqrt(n-3). 
## These will be the response variables in the models
d <- d %>% mutate(z_R = 0.5*log((1+R)/(1-R)), 
                  se_zR = 1/sqrt(SSInd-3))


#final big model with all variables
my_priors1 <- prior(normal(0,10), class = Intercept) +
  prior(normal(0,1), class = b)



bf10 <- bf(z_R | se(se_zR) ~ +1 + Method +  WildCap + Class +  TissueP + (1|Study)) + gaussian()

m10 <- brm(bf10, 
           data = d,
           prior = my_priors1,
           warmup = 1000,
           iter = 3500,
           chains = 4, 
           cores = 4,
           backend = "cmdstanr",
           threads = threading(2),
           control = list(adapt_delta = 0.8,
                          max_treedepth = 15),
           seed = 666,
           file = "m10")

pp_check(m10, ndraws=50)

#find probability of direction (pd)
pos<- bayestestR::describe_posterior(m10)

library(bayestestR)

bayestestR::p_direction(pos, method = "direct")
p_direction(m10)

###EMMEANS###

# # estimate estimated marginal means
# # estimated marginal means: show what the marginal means would be if the data set was balanced 
# # "pairs" estimates estimated mean difference
#make reference grid (data type to insert in emmeans)
rg<- ref_grid(m10)
m10d<-as.data.frame(m10)
#calculate emmeans for Methods
em<- emmeans(rg, "Method")
#summary
summary(em, point.est = mean)
#default plot
plot(em)
#plot pairwise comparisons
plot(pairs(em))
#turn into data frames
md<- as.data.frame(em)
mp<- as.data.frame(pairs(em))
#ggplots
m1<- ggplot()+geom_pointrange(data=md, aes(x=emmean, y=Method, xmin=lower.HPD, xmax=upper.HPD), size=0.8)+ylab("")+xlab("Estimated marginal mean")+theme(aspect.ratio = 1)
m1

##the above is still in the z-scale
##we need to transform everything to R and then R^2

#backtransform z to Rsq
##Methods
#z to R
#function to transform z to R
ztor<- function(x){
  r=(exp(2*z)-1)/(exp(2*z)+1)
  return(r)
}

#function to transform CIR to CIR2
CI <- function(CIR){
  CIR2 <- c(0.0,0.0)
  if (sign(CIR[1])*sign(CIR[2]) < 0){ # CIR contains 0
    CIR2[2] <- max(abs(CIR))^2
  }
  else{ # CIR does not contain 0
    CIR2[1] <- min(abs(CIR))^2
    CIR2[2] <- max(abs(CIR))^2
  }
  return(CIR2)
}

#turn the emmeans from z to R
z<- md[,2:4]
mdt<- ztor(x)
#R to R^2
#ONLY DQUARE R, NOT CI
mdt$emmean<-mdt$emmean^2
#turn into df
mdt<- as.data.frame(mdt)
#CI for R to CI for R^2
c1<-CI(c(0.9555524, 0.9926982))
c2<- CI(c(0.8681837,  0.9846428))
c3<- CI(c(0.8153254,0.9722469))
c4<- CI(c(0.7822901, 0.9728621))
c1<- as.data.frame(c1)
c2<- as.data.frame(c2)
c3<- as.data.frame(c3)
c4<- as.data.frame(c4)
c1<-t(c1)
c2<-t(c2)
c3<-t(c3)
c4<-t(c4)

#add method column
#put them together
mdt<- cbind(mdt, md$Method)
mci<- rbind(c1,c2,c3,c4)
mdt2<- cbind (mdt, mci)
colnames(mdt2)<- c("emmeans", "lower.HPD", "upper.HPD", "Method", "lower", "upper")
#graph
mt1<- ggplot()+geom_pointrange(data=mdt2, aes(x=emmeans, y=Method, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated marginal means")+theme(aspect.ratio = 1)+ theme(legend.position="none",
                                                                                                                                                                              axis.text.x=element_blank(),
                                                                                                                                                                              axis.ticks.x=element_blank(),
                                                                                                                                                                              plot.margin=unit(c(1,1,-0.5,1), "cm"))+xlim(0,1)
mt1

#METHOD
emm10 <- emmeans(m10, ~ Method)
d10 <- gather_emmeans_draws(emm10)
d10 <- mutate(d10,R2 = tanh(.value)^2) # inverse Fisher, squared, i.e. z->R^2

## Separate columns for methods:
d10w <- pivot_wider(d10, id_cols = .draw, names_from = Method, values_from = R2)

## Compute posteriors for contrasts:
d10wm <- mutate(d10w,
               TR = `Targeted BS` - RRBS,
               TO = `Targeted BS` - Other,
               RO = RRBS - Other,
               MT = Microarray - `Targeted BS`,
               MR = Microarray - RRBS,
               MO = Microarray - Other)

## Compute the CI limits and means
quantile(d10wm$TR,c(0.025,0.5,0.975))
mean(d10wm$TR)

quantile(d10wm$TO,c(0.025,0.5,0.975))
mean(d10wm$TO) 

quantile(d10wm$RO,c(0.025,0.5,0.975))
mean(d10wm$RO) 

quantile(d10wm$MT,c(0.025,0.5,0.975))
mean(d10wm$MT)

quantile(d10wm$MR,c(0.025,0.5,0.975))
mean(d10wm$MR) 

quantile(d10wm$MO,c(0.025,0.5,0.975))
mean(d10wm$MO) 


mp<- data.frame(matrix(ncol = 4, nrow = 6))
colnames(mp)<-c("contrast", "mean", "lower", "upper")
mp$contrast<- c("Targeted BS - RRBS",
                "Targeted BS - Other",
                "RRBS - Other",
                "Microarray - Targeted BS",
                "Microarray - RRBS",
                "Microarray - Other")
mp$mean<- c(mean(d10wm$TR),mean(d10wm$TO), mean(d10wm$RO), mean(d10wm$MT), mean(d10wm$MR), mean(d10wm$MO))
mp$lower<- c(-0.09166744 ,-0.12245596,-0.19251284, 0.003120076, 0.0246005, 0.01515179) 
mp$upper<- c(0.21191198, 0.29901777 , 0.26018594 ,0.198533842 ,0.2800393 , 0.35325423  )

mt2<- ggplot()+geom_pointrange(data=mp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)
mt2

###METHODS FINAL FIGURE###
ggarrange(mt1, mt2, nrow=2,ncol=1,align = "v")

###CLASS###
c<- emmeans(rg, "Class")
plot(c)
cd<- as.data.frame(c)

z<- cd[,2:4]
cdt<- ztor(x)
cdt$emmean<-cdt$emmean^2
c1<-CI(c(0.7966327, 0.9770496))
c2<- CI(c(0.7139015, 0.9953444))
c3<- CI(c(0.4839445, 0.9899830))
c4<- CI(c(0.8786488, 0.9962127))
c1<- as.data.frame(c1)
c2<- as.data.frame(c2)
c3<- as.data.frame(c3)
c4<- as.data.frame(c4)
c1<-t(c1)
c2<-t(c2)
c3<-t(c3)
c4<-t(c4)

#add method column
#put them together
cdt<- cbind(cdt, cd$Class)
cci<- rbind(c1,c2,c3,c4)
cdt2<- cbind (cdt, cci)
colnames(cdt2)<- c("emmeans", "lower.HPD", "upper.HPD", "Class", "lower", "upper")


ct1<- ggplot()+geom_pointrange(data=cdt2, aes(x=emmeans, y=Class, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated marginal mean")+theme(aspect.ratio = 1)+xlim(0,1)
ct1


#contrast
emc10 <- emmeans(m10, ~ Class)
c10 <- gather_emmeans_draws(emc10)
cc10 <- mutate(c10,R2 = tanh(.value)^2) # inverse Fisher, squared, i.e. z->R^2

## Separate columns for methods:
c10w <- pivot_wider(cc10, id_cols = .draw, names_from = Class, values_from = R2)

## Compute posteriors for contrasts:
c10wm <- mutate(c10w,
                MO= Mammal-Other,
                MF = Mammal-Fish,
                MB = Mammal-Bird,
                FO = Fish - Other,
                FB = Fish-Bird,
                BO = Bird-Other)

## Compute the CI limits and means
quantile(c10wm$MO,c(0.025,0.5,0.975))

quantile(c10wm$MF,c(0.025,0.5,0.975))

quantile(c10wm$MB,c(0.025,0.5,0.975))

quantile(c10wm$FO,c(0.025,0.5,0.975))

quantile(c10wm$FB,c(0.025,0.5,0.975))

quantile(c10wm$BO,c(0.025,0.5,0.975))

cp<- data.frame(matrix(ncol = 4, nrow = 6))
colnames(mp)<-c("contrast", "mean", "lower", "upper")
cp$contrast<- c("Mammal-Other",
                "Mammal-Fish",
                "Mammal-Bird",
                "Fish - Other",
                "Fish-Bird",
                "Bird-Other")
cp$mean<- c(mean(c10wm$MO),mean(c10wm$MF), mean(c10wm$MB), mean(c10wm$FO), mean(c10wm$FB), mean(c10wm$BO))
cp$lower<- c(-0.3077245 ,-0.33388491,-0.1783943,-0.48221740, -0.39461272, -0.7247790) 
cp$upper<- c(0.0798490, 0.41769365 , 0.5695092 ,0.19005979 ,0.72688448 , 0.1038128)

ct2<- ggplot()+geom_pointrange(data=cp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)
ct2

###CLASS FINAL FIGURE###
ggarrange(ct1,ct2, nrow=2, ncol=1, align="v")

###TISSUE###
t<- emmeans(rg, "TissueP")
plot(t)
td<- as.data.frame(t)

z<- td[,2:4]
tdt<- ztor(x)
tdt$emmean<-tdt$emmean^2
c1<-CI(c(0.8453860, 0.9712804))
c2<- CI(c(0.6677186,  0.9962868))
c3<- CI(c(0.8650932, 0.9849398))

c1<- as.data.frame(c1)
c2<- as.data.frame(c2)
c3<- as.data.frame(c3)

c1<-t(c1)
c2<-t(c2)
c3<-t(c3)

tdt<- cbind(tdt, td$TissueP)
tci<- rbind(c1,c2,c3)
tdt2<- cbind (tdt, tci)
colnames(tdt2)<- c("emmeans", "lower.HPD", "upper.HPD", "Tissue", "lower", "upper")


tt1<- ggplot()+geom_pointrange(data=tdt2, aes(x=emmeans, y=Tissue, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated marginal mean")+theme(aspect.ratio = 1)+xlim(0,1)
tt1

##contrast
emt10 <- emmeans(m10, ~ TissueP)
t10 <- gather_emmeans_draws(emt10)
tc10 <- mutate(t10,R2 = tanh(.value)^2) # inverse Fisher, squared, i.e. z->R^2

## Separate columns for methods:
t10w <- pivot_wider(tc10, id_cols = .draw, names_from = TissueP, values_from = R2)

## Compute posteriors for contrasts:
t10wm <- mutate(t10w,
                PM= Proliferative-Multiple,
                PN= Proliferative-`Non-proliferative`,
                NM= `Non-proliferative`- Multiple)

## Compute the CI limits and means
quantile(t10wm$PM,c(0.025,0.5,0.975))

quantile(t10wm$PN,c(0.025,0.5,0.975))

quantile(t10wm$NM,c(0.025,0.5,0.975))



tp<- data.frame(matrix(ncol = 4, nrow = 3))

colnames(mp)<-c("contrast", "mean", "lower", "upper")
tp$contrast<- c( "Proliferative-Multiple",
                "Proliferative- Non-proliferative",
               " Non-proliferative- Multiple")

tp$mean<- c(mean(t10wm$PM),mean(t10wm$PN), mean(t10wm$NM))

tp$lower<- c(-0.14052198,-0.26378752,-0.47679450) 
tp$upper<- c(0.06234251,0.44786158,0.22576250)

tt2<- ggplot()+geom_pointrange(data=tp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)
tt2

###FINAL TISSUE GRAPH###
ggarrange(tt1,tt2, nrow=2, ncol=1, align="v")

##WILD AND CAPTIVE##
w<- emmeans(rg, "WildCap")
plot(w)
wd<- as.data.frame(w)

z<- wd[,2:4]
wdt<- ztor(x)
wdt$emmean<-wdt$emmean^2
c1<-CI(c(0.9104035,0.9818999))
c2<- CI(c(0.8733297, 0.9768131))
c3<- CI(c(0.8263826, 0.9856564))

c1<- as.data.frame(c1)
c2<- as.data.frame(c2)
c3<- as.data.frame(c3)

c1<-t(c1)
c2<-t(c2)
c3<-t(c3)

wdt<- cbind(wdt, wd$WildCap)
wci<- rbind(c1,c2,c3)
wdt2<- cbind (wdt, wci)
colnames(wdt2)<- c("emmeans", "lower.HPD", "upper.HPD", "WildCap", "lower", "upper")


wt1<- ggplot()+geom_pointrange(data=wdt2, aes(x=emmeans, y=WildCap, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated marginal mean")+theme(aspect.ratio = 1)+xlim(0,1)
wt1

##contrast
emw10 <- emmeans(m10, ~ WildCap)
w10 <- gather_emmeans_draws(emw10)
wc10 <- mutate(w10,R2 = tanh(.value)^2) # inverse Fisher, squared, i.e. z->R^2

## Separate columns for methods:
w10w <- pivot_wider(wc10, id_cols = .draw, names_from = WildCap, values_from = R2)

## Compute posteriors for contrasts:
w10wm <- mutate(w10w,
                CW= Captive-Wild,
                BW= Both-Wild,
                BC= Both- Captive)

## Compute the CI limits and means
quantile(w10wm$CW,c(0.025,0.5,0.975))

quantile(w10wm$BW,c(0.025,0.5,0.975))

quantile(w10wm$BC,c(0.025,0.5,0.975))



wp<- data.frame(matrix(ncol = 4, nrow = 3))

colnames(wp)<-c("contrast", "mean", "lower", "upper")
wp$contrast<- c( "Captive-Wild",
                 "Both-Wild",
                 "Both- Captive")

wp$mean<- c(mean(w10wm$CW),mean(w10wm$BW), mean(w10wm$BC))

wp$lower<- c(-0.02823365,-0.205834985,-0.22848214) 
wp$upper<- c(0.11848058,0.153149008 ,0.08671781)

wt2<- ggplot()+geom_pointrange(data=wp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)
wt2

###FINAL GRID OF ALL GRAPHS###
mt1<- ggplot()+geom_pointrange(data=mdt2, aes(x=emmeans, y=Method, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("")+theme(aspect.ratio = 1)+ coord_cartesian(xlim=c(0, 1))+scale_x_continuous(breaks=seq(0, 1, 0.3))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
mt1

mt2<- ggplot()+geom_pointrange(data=mp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)+xlim(-0.8,0.8)+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
mt2


grid.arrange(mt1, ct1, heights=c(1.2,1.2), widths=c(1.4))


ct1<- ggplot()+geom_pointrange(data=cdt2, aes(x=emmeans, y=Class, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab(NULL)+theme(aspect.ratio = 1)+ coord_cartesian(xlim=c(0, 1))+scale_x_continuous(breaks=seq(0, 1, 0.3))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
ct1

ct2<- ggplot()+geom_pointrange(data=cp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)+xlim(-0.8,0.8)+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
ct2

tt1<- ggplot()+geom_pointrange(data=tdt2, aes(x=emmeans, y=Tissue, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab(NULL)+theme(aspect.ratio = 1)+ coord_cartesian(xlim=c(0, 1))+scale_x_continuous(breaks=seq(0, 1, 0.3))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
tt1

tt2<- ggplot()+geom_pointrange(data=tp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)+xlim(-0.8,0.8)+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
tt2

wt1<- ggplot()+geom_pointrange(data=wdt2, aes(x=emmeans, y=WildCap, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated marginal mean R-sq")+theme(aspect.ratio = 1)+ coord_cartesian(xlim=c(0, 1))+scale_x_continuous(breaks=seq(0, 1, 0.3))
wt1

wt2<- ggplot()+geom_pointrange(data=wp, aes(x=mean, y=contrast, xmin=lower, xmax=upper), size=0.8)+ylab("")+xlab("Estimated mean difference R-sq")+theme(aspect.ratio = 1)+geom_vline(xintercept=0, linetype="dashed",  color = "black", size=0.4)+xlim(-0.8,0.8)
wt2

tiff("emmeansAll.tiff", units="in", width=6, height=10, res=300)
plot_grid(mt1, tt1, wt1, ncol=1, align="v",label_size=18,labels=c("a"), label_x=0.25, rel_heights=c(0.25,0.25,0.25,0.32), rel_widths =c(0.25,0.25,0.25,0.32))
dev.off()

tiff("meandifAll.tiff", units="in", width=8, height=10, res=300)
plot_grid(mt2,ct2,tt2,wt2, align="v" ,ncol = 1,label_size=18,labels=c("b"), label_x=0.3, rel_heights=c(0.25,0.25,0.25,0.32), rel_widths =c(0.25,0.25,0.25,0.32))
dev.off()

##Supplementaty information
#conditional affects +-se
ce<- conditional_effects(m10)
plot(ce)
Method<- as.data.frame(ce$Method)

method<- ggplot(data=Method, aes(x=Method, y=estimate__))+geom_errorbar(aes(ymin=estimate__-se__, ymax=estimate__+se__), width=0.2)+geom_point(size=3)+ylab("Conditional effects of \n z-transformed coefficients of determination")+xlab("")+ylim(0.8,3.5)

WildCap<- as.data.frame(ce$WildCap)
wildcap<- ggplot(data=WildCap, aes(x=WildCap, y=estimate__))+geom_errorbar(aes(ymin=estimate__-se__, ymax=estimate__+se__), width=0.2)+geom_point(size=3)+ylab("Conditional effects of \n z-transformed coefficients of determination")+xlab("")+ylim(0.8,3.5)
wildcap

Class<- as.data.frame(ce$Class)
class<- ggplot(data=Class, aes(x=Class, y=estimate__))+geom_errorbar(aes(ymin=estimate__-se__, ymax=estimate__+se__), width=0.2)+geom_point(size=3)+ylab("Conditional effects of \n z-transformed coefficients of determination")+xlab("")+ylim(0.8,3.5)
class

Tissue<- as.data.frame(ce$TissueP)
tissue<- ggplot(data=Tissue, aes(x=TissueP, y=estimate__))+geom_errorbar(aes(ymin=estimate__-se__, ymax=estimate__+se__), width=0.2)+geom_point(size=3)+ylab("Conditional effects of \n z-transformed coefficients of determination")+xlab("")+ylim(0.8,3.5)
tissue

tiff("ConditionalEffects.tiff", units="in", width=12, height=12, res=300)
ggarrange(method, wildcap, class, tissue, nrow=2, ncol=2,labels="AUTO")
dev.off()
