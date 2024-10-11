#This is an R script of the graphs included in the Review paper on epigenetic clocks in non-model animals
#https://onlinelibrary.wiley.com/doi/10.1111/mec.17065
#Author: Marianthi Tangili
#Last updated 18/04/2023

setwd("C:/Users/MWP-/Desktop/PhD/Review_DNAmWild/Figures/FINAFiguresJune2023")

#load necessary packages
library(readxl)
library(ggplot2)
library(ggthemes)
library(janitor)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(dplyr)
library(smplot)
library(devtools)

#load graph theme
graph_theme<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour =NA ,fill=NA,size=0.5),
  axis.title.x=element_text(size=24,hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=24,hjust=0.5,vjust=0.5,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=24),
  axis.text.y=element_text(colour="black",angle=0,size=24),
  axis.ticks=element_line(colour="black",size=0.8),
  axis.line.x=element_line(colour="black",size=0.8),
  axis.line.y=element_line(colour="black",size=0.8))

#read in data file all studies
#"SumRev____.xlsx" is a summary I made because the original Excel file was hard to work with in R
#This is the file I keep updated (adding new studies etc.)
SumRev <- read_excel("C:/Users/MWP-/Desktop/PhD/Review_DNAmWild/Review_Analysis/DNAm_NonModel_Database.xlsx",na="NA", sheet="AllClocksS")
View(SumRev)
names(SumRev)
str(SumRev)


#read in data file for only papers
SinglePap <- read_excel("C:/Users/MWP-/Desktop/SumRevOct2022.xlsx", 
                        sheet = "SinglePaper")
View(SinglePap)

#turn "NA" text in the file to NA=missing values
SumRev[SumRev=="NA"] <- NA
SinglePap[SinglePap=="NA"] <- NA

#transform variables
SumRev$Rsq<- as.numeric(SumRev$Rsq)
SumRev$ClockCPGs<- as.numeric(SumRev$ClockCPGs)
SumRev$SSInd<- as.numeric(SumRev$SSInd)


### FIGURE 3 ###
# Year, Class, Goal 

#PAPERS PER YEAR
#Make "Year" numeric
SinglePap$Year<- as.numeric(SinglePap$Year)
#find how many studies have available Year for N=
sum(!is.na(SinglePap$Year))
#find sample size per year to annotate
table(SinglePap$Year)
#histogram of papers on epigenetic clocks per year (ggplot)
y<- ggplot(SinglePap, aes(x=Year))+ geom_bar(fill="#006666", color="black",size=2) + ylab("Number of studies")+xlab(NULL) + annotate(geom="text", x=2014, y=4, label = paste("2"), size=12)+annotate(geom="text", x=2016, y=3, label = paste("1"), size=12)+annotate(geom="text", x=2017, y=5, label = paste("3"), size=12)+annotate(geom="text", x=2018, y=4, label = paste("2"), size=12)+annotate(geom="text", x=2019, y=7, label = paste("5"), size=12)+annotate(geom="text", x=2020, y=10, label = paste("8"), size=12)+annotate(geom="text", x=2021, y=20, label = paste("18"), size=12)+annotate(geom="text", x=2022, y=14, label = paste("12"), size=12)+scale_y_continuous(expand=c(0,0), limits=c(0,22),breaks=seq(0,22, by=5)) +scale_x_continuous(labels=c("2014","2015", "2016","2017", "2018", "2019", "2020", "2021", "2022 \n (6 months)"), breaks = 2014:2022)+ theme(axis.text.x = element_text(size = 24),axis.text.y = element_text(size = 24), axis.title.x=element_text(size=24), axis.title.y=element_text(size=24))
y

##CLASS
#sample size per Class
Class<- tabyl(SinglePap, Class)
#barplot of No epigenetic clocks per class
c<- ggplot(data=Class, aes(x=Class, y=n))+geom_bar(stat="identity",fill="#006666", color="black",size=2)+ylim(0,40)+ylab("Number of Studies")+xlab(NULL)+scale_y_continuous(expand=c(0,0), limits=c(0,46),breaks=seq(0,46, by=5)) + theme(axis.text.x = element_text(size = 22, angle = 45, hjust=1),axis.text.y = element_text(size = 22), axis.title.x=element_text(size=24), axis.title.y=element_text( margin = margin(t = 0, r = 15, b = 0, l = 0), size=24), axis.text = element_text(size = 20))+annotate(geom="text", x=1, y=9, label = paste("6"), size=12)+annotate(geom="text", x=2, y=7, label = paste("4"), size=12)+annotate(geom="text", x=3, y=41, label = paste("38"), size=12)+annotate(geom="text", x=4, y=6, label = paste("3"), size=12)
c
##GOALS
Goals<- tabyl(SinglePap, Goal)
#barplot No epigenetic clocks per Goal
g<- ggplot(data=Goals, aes(x=Goal, y=n))+geom_bar(stat="identity",fill="#006666", color="black",size=2)+ylim(0,40)+ylab(NULL)+xlab(NULL)+scale_y_continuous(expand=c(0,0), limits=c(0,46),breaks=seq(0,46, by=5))+ theme(axis.text.x = element_text(size = 22, angle = 45, hjust=1),axis.text.y = element_text(size = 22), axis.title.x=element_text(size=24), axis.title.y=element_text( margin = margin(t = 0, r = 15, b = 0, l = 0), size=24), axis.text = element_text(size = 20))+annotate(geom="text", x=1, y=9, label = paste("6"), size=12)+annotate(geom="text", x=2, y=8, label = paste("5"), size=12)+annotate(geom="text", x=3, y=39 ,label = paste("36"), size=12)+annotate(geom="text", x=4, y=7, label = paste("4"), size=12)
g

#arrange all 3 figures in 1

gt<- ggarrange(y,
               ggarrange(c, g, ncol = 2, align = "h",widths = c(1.5,2), labels=c("(b)", "(c)"), label.x = c(0.22,0.08),font.label = list(size = 22)),
               nrow = 2, 
               heights = c(1, 1.5, 1.5), labels=c("(a)"), font.label = list(size = 22), label.x = 0.085)

#save Fig.3
CairoPDF("Figure_3.pdf", width=10, height=12)
gt
dev.off()



### FIGURE 4 ###
##Method, tissue, WildCap##
##METHODS
Methods<- tabyl(SumRev, Method)
m<- ggplot(data=Methods, aes(x=Method, y=n))+geom_bar(stat="identity",fill="#006666", color="black", size=2)+ylim(0,30)+ylab("Frequency")+xlab(NULL)+theme(axis.text = element_text(size =22))+scale_y_continuous(expand=c(0,0), limits=c(0,45),breaks=seq(0,45, by=5)) + theme(axis.text.x = element_text(size = 22,angle = 45, hjust=1),axis.text.y = element_text(size = 22), axis.title.x=element_text(size=24), axis.title.y=element_text(size=24))+annotate(geom="text", x=1, y=31, label = paste("28"), size=12)+annotate(geom="text", x=2, y=7, label = paste("4"), size=12)+annotate(geom="text", x=3, y=11, label = paste("8"), size=12)+annotate(geom="text", x=4, y=14, label = paste("11"), size=12)+annotate(geom="text", x=5, y=7, label = paste("4"), size=12)
m

TissueP<- tabyl(SumRev, Prol)
TissueP$Prol<- factor(TissueP$Prol, levels=c("Proliferative", "Non-proliferative", "Multiple"))
t<- ggplot(data=TissueP, aes(x=Prol, y=n))+geom_bar(stat="identity",fill="#006666", color="black", size=2)+ylab("Number of Epigenetic Clocks")+xlab(NULL)+scale_y_continuous(expand=c(0,0), limits=c(0,45),breaks=seq(0,45, by=5)) + theme(axis.text.x = element_text(size = 22, angle = 45, hjust=1),axis.text.y = element_text(size = 24), axis.title.x=element_text(size=24), axis.title.y=element_text(size=24))+ annotate(geom="text", x=3, y=14, label = paste("11"), size=12)+ annotate(geom="text", x=2, y=8, label = paste("5"), size=12)+ annotate(geom="text", x=1, y=42, label = paste("39"), size=12)
t

WC<- tabyl(SumRev,WildCap)
WC$WildCap<- factor(WC$WildCap, levels=c("Captive", "Wild", "Both"))
w<- ggplot(data=WC, aes(x=WildCap, y=n))+geom_bar(stat="identity",fill="#006666", color="black", size=2)+ylab(NULL)+xlab(NULL)+scale_y_continuous(expand=c(0,0), limits=c(0,45),breaks=seq(0,45, by=5)) + theme(axis.text.x = element_text(size = 22, angle = 45, hjust=1),axis.text.y = element_text(size = 24), axis.title.x=element_text(size=24), axis.title.y=element_text(size=24))+ annotate(geom="text", x=3, y=7, label = paste("4"), size=12)+ annotate(geom="text", x=1, y=35, label = paste("32"), size=12)+ annotate(geom="text", x=2, y=22, label = paste("19"), size=12)
w

#arrange 6 figures in 1
ycgmtw<- ggarrange(y, 
          ggarrange(c, g, ncol = 2, align="h", labels=c("(b)", "(c)"), label.x = c(0.15,0.07),font.label = list(size=24)), 
          ggarrange(m, t, w, ncol = 3, align="h", labels=c("(d)", "(e)", "(f)"), label.x=c(0.16,0.15,0.11), font.label = list(size=24)), ncol = 1, labels="(a)", label.x=0.07, font.label = list(size=24))

ggsave(file="Figure_4.pdf",ycgmtw, width=5000, height=5500, units="px")

###FIGURE 5###
##Rsq histogram##
#find how many R-sq are actually reported
sum(!is.na(SumRev$Rsq))
#histogram of R-sq distribution
rsq<- ggplot(SumRev, aes(x=Rsq))+ geom_histogram(fill="#006666", color="black",size=1.5, bins = 14) +xlab("Coefficient of determination" ~ (R^2))+ylab("Number of epigenetic clocks")+theme(axis.text = element_text(size =20))+scale_y_continuous(expand=c(0,0),breaks=seq(0,24, by=2),limits=c(0,26))+scale_x_continuous(limits=c(0,1.02),breaks=seq(0,1.02, by=0.2)) +annotate(geom="text", x=0.5, y=14, label = paste("N=42"), size=10)+ theme(aspect.ratio=1/1, panel.border = element_rect(colour = "black", fill=NA, size=2))
rsq
sum(!is.na(SumRev$PropMAE))

mad<- ggplot(SumRev, aes(x=PropMAE))+ geom_histogram(fill="#006666", color="black",size=1.5, bins = 14) +xlab("Scaled Mean/Median Absolute Deviation (MAD)")+ylab(NULL)+theme(axis.text = element_text(size =24))+scale_y_continuous(expand=c(0,0),breaks=seq(0,24, by=2),limits=c(0,26))+scale_x_continuous(limits=c(0,1.02),breaks=seq(0,1.02, by=0.2)) +annotate(geom="text", x=0.5, y=14, label = paste("N=36"), size=10)+ theme(aspect.ratio=1/1, panel.border = element_rect(colour = "black", fill=NA, size=2))
mad

sum(!is.na(SumRev$ScaleMAD))
mads<- ggplot(SumRev, aes(x=ScaleMAD))+ geom_histogram(fill="#006666", color="black",size=1.5, bins = 14) +xlab("Scaled Mean/Median Absolute Deviation (MAD)")+ylab(NULL)+theme(axis.text = element_text(size =24))+scale_y_continuous(expand=c(0,0),breaks=seq(0,10, by=2),limits=c(0,10))+scale_x_continuous(limits=c(-1.8,-0.25),breaks=seq(-1.8,-0.25, by=0.25)) +annotate(geom="text", x=0.5, y=14, label = paste("N=36"), size=10)+ theme(aspect.ratio=1/1, panel.border = element_rect(colour = "black", fill=NA, size=2))
mads

sum(!is.na(SumRev$MADar))
madR<- ggplot(SumRev, aes(x=MADar))+ geom_histogram(fill="#006666", color="black",size=1.5, bins = 14) +xlab("Scaled Mean/Median Absolute Deviation (MAD)")+ylab(NULL)+theme(axis.text = element_text(size =24))+scale_y_continuous(expand=c(0,0),breaks=seq(0,24, by=2),limits=c(0,26))+scale_x_continuous(limits=c(0,1),breaks=seq(0,1, by=0.25)) +annotate(geom="text", x=0.5, y=14, label = paste("N=36"), size=10)+ theme(aspect.ratio=1/1, panel.border = element_rect(colour = "black", fill=NA, size=2))
madR

tiff("RsqMADPanel.tiff", units="in", width=14, height=8, res=300)
rm<- ggarrange(rsq, madR, labels=c("(a)", "(b)"), align="h",font.label = list(size=24))+  theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
dev.off()


ggsave(file="Figure_4.pdf",rm,  width=4200, height=1800, units="px")

cor.test(noW$Rsq, noW$ScaleMAD, method = c("spearman"))


### FIGURE 5 ###

##ClockCpGs Histogram##
sum(!is.na(SumRev$ClockCPGs))
ggplot(SumRev, aes(x=ClockCPGs))+ geom_histogram(fill="white", color="black",size=1, binwidth = 14)+annotate(geom="text", x=300, y=3, label = paste("N=31"), size=6)+scale_x_continuous(expand=c(0.02,0), breaks=seq(0,640, by=100))+ylab("Number of Epigenetic Clocks")+xlab("Number of CpGs used for clock construction") +theme(axis.text = element_text(size =20))+scale_y_continuous(expand=c(0,0), limits=c(0,4.5), breaks=seq(0,4.5, by=1))

### FIGURE 6 ###
##clockCpGs histogram##

h<- ggplot(SumRev, aes(x=ClockCPGs))+ geom_histogram(fill="#006666", color="black",size=1, binwidth = 14)+annotate(geom="text", x=300, y=3, label = paste("N=31"), size=8)+scale_x_continuous(expand=c(0.02,0), breaks=seq(0,640, by=100))+ylab("Number of Epigenetic Clocks")+xlab("Number of CpGs used for clock construction") +scale_y_continuous(expand=c(0,0), limits=c(0,6.5), breaks=seq(0,6.5, by=1))+ theme(aspect.ratio=1/1, panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(axis.title.y = element_text(vjust = 3), axis.title.x = element_text(vjust = 0))
h

##clockCpGs boxplot##
b<-ggplot(SumRev, aes(y=ClockCPGs))+geom_boxplot(size=1, outlier.shape = NA, color="black")+geom_jitter(aes(x=0), size=4, alpha=0.8, col="#006666")+ 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+ylab("CpG sites used for clock construction")+xlab("")+ theme(text = element_text(size = 20))

b

tiff("ClockCpGs2.tiff", units="in", width=14, height=10, res=300)
ggarrange(h,b, ncol=2, nrow=1, widths = c(2, 1), labels = c("(a)", "(b)"),font.label = list(size=26), label.x = c(-0.01,-0.1))+  theme(plot.margin = unit(c(1, 0, 0, 0), "cm"))
dev.off()

ggsave(file="Figure_5.pdf",h,  width=2200, height=1800, units="px")



### Figure 7 ###
#Rsq + MAE vs. clockCpGs, SSInd and AgeRange
noW<- subset(SumRev, Method!="WGBS")


library(viridis)

colors<- viridis(4)
shapes<- c(16,17,18,15)


##Sample size graphs
rss<- ggplot(data=noW, aes(x=as.numeric(SSInd), y=as.numeric(Rsq)))+geom_point(size=8 , aes(shape=Method, col=Method), alpha=0.8)+ scale_shape_manual(values = shapes)+ theme(axis.line = element_line(colour = "gray", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+ xlab(NULL) + scale_x_continuous(trans='log10')+ylab("Coefficient of determination" ~ (R^2))+ ylim(0.5,1)+geom_smooth(method=lm,se=FALSE,fullrange=TRUE,linetype="dashed",color="black")+ theme(plot.margin = unit(c(3,0,0,0), "lines"))+  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 28)) + scale_color_manual(values = colors)+ scale_shape_manual(values = shapes)
rss

cor.test(noW$Rsq, log(noW$SSInd), method = c("pearson"))

mss<- ggplot(noW, aes(x=SSInd, y=ScaleMAD))+geom_point(size=8, aes(shape=Method, col=Method), alpha=0.8)+ theme(axis.line = element_line(colour = "gray", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+ xlab("Individuals sampled") + scale_x_continuous(trans='log10')+ ylab("Scaled MAD")+ ylim(-2,-0.4)+geom_smooth(method=lm,se=FALSE,fullrange=TRUE,linetype="solid",color="black")+ theme(plot.margin = unit(c(3,0,0,0), "lines"))+  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 28)) + scale_color_manual(values = colors)+ scale_shape_manual(values = shapes)
mss

cor.test(noW$ScaleMAD, log(noW$SSInd), method = c("pearson"))



#ClockCpG graphs
rclocks<- ggplot(noW, aes(x=as.numeric(ClockCPGs), y=as.numeric(Rsq)))+geom_point(size=8, aes(shape=Method, col=Method),, alpha=0.8) + theme(axis.line = element_line(colour = "gray", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+ xlab(NULL) + scale_x_continuous(trans='log10')+ ylab(NULL)+ ylim(0.5,1)+geom_smooth(method=lm,se=FALSE,fullrange=TRUE,linetype="dashed",color="black")+ theme(plot.margin = unit(c(3,0,0,0), "lines"))+  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 28))+ scale_color_manual(values = colors)+ scale_shape_manual(values = shapes)
rclocks

cor.test(noW$Rsq, log(noW$ClockCPGs), method = c("pearson"))


mclocks<- ggplot(noW, aes(x=as.numeric(ClockCPGs), y=as.numeric(ScaleMAD)))+geom_point(size=8, aes(shape=Method, col=Method), alpha=0.8) + theme(axis.line = element_line(colour = "grey", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+ xlab("Number of CpG sites used for clock") + scale_x_continuous(trans='log10')+ ylab(NULL)+ ylim(-2,-0.4)+geom_smooth(method=lm,se=FALSE,fullrange=TRUE,linetype="solid",color="black")+ theme(plot.margin = unit(c(3,0,0,0), "lines"))+theme(axis.title.x=element_text(size=22.3))+  theme(legend.title = element_text(size = 22.3), legend.text = element_text(size = 28))+ scale_color_manual(values = colors)+ scale_shape_manual(values = shapes)
mclocks
cor.test(noW$ScaleMAD, log(noW$ClockCPGs), method = c("pearson"))

#Age Range graphs

arr<- ggplot(data=noW, aes(x=as.numeric(AgeRangeN), y=Rsq))+geom_point(size=8, aes(shape=Method, col=Method), alpha=0.8)+ theme(axis.line = element_line(colour = "gray", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+ xlab(NULL) + scale_x_continuous(trans='log10')+ ylab(NULL)+ ylim(0.5, 1)+geom_smooth(method=lm,se=FALSE,fullrange=TRUE,linetype="dashed",color="black")+ theme(plot.margin = unit(c(3,0,0,0), "lines"))+  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 28))+ scale_color_manual(values = colors)+ scale_shape_manual(values = shapes)
arr

arm<- ggplot(data=noW,aes(x=as.numeric(AgeRangeN), y=ScaleMAD))+geom_point(size=8, aes(shape=Method, col=Method), alpha=0.8)+ theme(axis.line = element_line(colour = "gray", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+ xlab("Age range (years)") + scale_x_continuous(trans='log10')+ ylab(NULL)+ ylim(-2,-0.4) +geom_smooth(method=lm,se=FALSE,fullrange=TRUE,linetype="dashed",color="black")+ theme(plot.margin = unit(c(3,0,0,0), "lines"))+  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 28))+ scale_color_manual(values = colors)+ scale_shape_manual(values = shapes)

arm
noW$AgeRangeN<- as.numeric(noW$AgeRangeN)
cor.test(noW$ScaleMAD, log(noW$AgeRangeN), method=c("pearson"))

#Save Figure 7

tiff("RsqMADPanel.tiff", units="in", width=18, height=14, res=300)
f7<- ggarrange(rss, rclocks, arr, mss, mclocks,arm, nrow=2, ncol=3, align="v",  common.legend = TRUE, legend = "top",labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), label.x = c(0.1,0.1,0.1, 0.1, 0.1, 0.1), label.y = c(1,1,1,1, 1, 1), font.label = list(size=26))+ 
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"))
f7
dev.off()

ggsave(file="Figure_7.pdf",f7,  width=6000, height=4000, units="px")
ggsave(file="Figure_7.png",f7,  width=6000, height=4000, units="px")

#Figure 6
#Clock CpGs vs. SSInd
css<- ggplot(data=SumRev, aes(x=SSInd, y=ClockCPGs))+geom_point(size=3)+ theme(axis.line = element_line(colour = "gray", linetype = "solid"), axis.ticks = element_line(linetype = "twodash"), panel.background = element_rect(fill = NA),plot.background = element_rect(fill = "white"))+  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+geom_smooth(method=lm,se=T,fullrange=TRUE,linetype="solid",color="black", fill="#006666")+ylab("Number of CpG sites used for clock construction")+xlab("Number of individuals sampled")+ scale_x_continuous(trans='log10')+ scale_y_continuous(trans='log10')+theme(aspect.ratio = 1/1, panel.border = element_rect(colour = "black", fill=NA, size=1))
css

tiff("ClockCpGsSSIng.tiff", units="in", width=16, height=12, res=300)
ggscatter(SumRev, x="SSInd", y="ClockCPGs",xscale="log10", ylim=c(-10,600),
add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "gray"), # Customize reg. line
conf.int = TRUE, 
cor.coef = TRUE, 
cor.coeff.args = list(method = "spearman", label.x = 1, label.y=500,  label.sep = "\n"),
xlab="Number of Individuals sampled", ylab="Number of CpG sites used for \n clock construction")
dev.off()

#Save Figure 8
tiff("ClockCpGsSSIng.tiff", units="in", width=8, height=8, res=300)
css
dev.off()

ggsave(file="Figure_6.pdf",css,  width=2000, height=2000, units="px")
ggsave(file="Figure_6.png",css,  width=2000, height=2000, units="px")

#Supplementary figure

tiff("MADRsq.tiff", units="in", width=8, height=6, res=300)
S1<- ggscatter(noW, x="ScaleMAD", y="Rsq", size=3,shape="Method",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "gray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = -1.5, label.y=0.7, label.sep = "\n"), xlab="Scaled MAD", ylab="Coefficient of determination (R-sq)")
dev.off()

ggsave(file="Figure_S1.pdf",S1,  width=2000, height=2000, units="px")

cor.test(noW$Rsq,noW$ScaleMAD , method=c( "spearman"))


#model Clock CpGs SsInd
cor.test(log(SumRev$ClockCPGs), log(SumRev$SSInd))


