library(dplyr)
library(agricolae)
library(cowplot)
library(egg)
library(car)
library(ggplot2)
library(ggsignif)
library(forcats)

setwd(dir="C:/Users/Nans/Desktop/PVC_MPs_G.fossarum/Data/(H3) Shredding activities")
AE_FR<-read.csv("AE_FR.csv", header = T, sep = ";")
AE_FR$Size <- fct_relevel(AE_FR$Size, c("Small", "Large"))
AE_FR$Concentration<- fct_relevel(AE_FR$Concentration, c("Control","Low", "High", "Extreme"))
AE_FR$Treatment <- fct_relevel(AE_FR$Treatment, c("GLS", "GMS", "GHS","GLL","GML","GHL","CG-S","CG-L"))
View(AE_FR)


#################Statistical analysis

####For the Assimilation Efficiency (AE)

##Two-Way ANOVA 

res.aov1 <- aov(AE ~ Size*Concentration, data = AE_FR)
leveneTest(AE ~ Size*Concentration, data = AE_FR)
plot(res.aov1, 2)
aov_residuals <- residuals(object = res.aov1)
shapiro.test(x = aov_residuals)
summary(res.aov1)
TukeyHSD(res.aov1, group=TRUE)

##Multiple comparison with Tukey HSD

#Comparison between MP sizes and concentrations

tx<- with(AE_FR,interaction(Size,Concentration))
amod <- aov(AE ~ tx, data=AE_FR)
HSD.test(amod, "tx", group=TRUE, console=TRUE)

#Comparison between MP Concentrations

tz<- with(AE_FR,interaction(Concentration))
amodz <- aov(AE~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)

#Comparison between MP Sizes

tz<- with(AE_FR,interaction(Size))
amodz <- aov(AE~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)

#Comparison between the different exposure treatments

tz<- with(AE_FR,interaction(Treatment))
amodz <- aov(AE~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)




####For the Feeding Rate (FR)

res.aov2 <- aov(FR ~ Size*Concentration, data = AE_FR)
leveneTest(FR ~ Size*Concentration, data = AE_FR)

# leveneTest (p=O.04) = Homogeneity of variance is not verified. Thus we transform the FR data :

AE_FR$FR <- log(AE_FR$FR)
res.aov2 <- aov(FR ~ Size*Concentration, data = AE_FR)
leveneTest(FR ~ Size*Concentration, data = AE_FR)
plot(res.aov2, 2)
aov_residuals2 <- residuals(object = res.aov2)
shapiro.test(x = aov_residuals2)
summary(res.aov2)
TukeyHSD(res.aov2)

##Multiple comparison with Tukey HSD

#Comparison between MP sizes and concentrations

tz<- with(AE_FR,interaction(Size,Concentration))
amodz <- aov(FR ~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)

#Comparison between MP concentrations

tz<- with(AE_FR,interaction(Concentration))
amodz <- aov(FR ~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)

#Comparison between MP Sizes

tz<- with(AE_FR,interaction(Size))
amodz <- aov(FR ~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)

#Comparison between the different exposure treatments

tz<- with(AE_FR,interaction(Treatment))
amodz <- aov(FR ~ tz, data=AE_FR)
HSD.test(amodz, "tz", group=TRUE, console=TRUE)


#################Graphical representations

####For the Assimilation Efficiency (AE)
  
AEgraph<-ggplot(data=AE_FR, aes(x=Size, y=AE, color=Size)) + theme_article()+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Assimilation efficiency") + ylab("Assimilation efficiency (%)") +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1, alpha=0.5)
AEgraph<- AEgraph + theme(legend.position = "none")

AEgraph<-AEgraph+scale_color_manual(values=c("royalblue4","springgreen4")) 

AEgraph_FINAL<-AEgraph + facet_grid(cols = vars(Concentration))+
  theme(strip.background = element_rect(colour="black", fill="white", size=0.5, linetype="solid"))+ 
  theme(panel.border = element_rect(fill = "transparent",color = "black", linewidth = 0.5)) 


AEgraph_FINAL

####For the Feeding Rate (FR)

FRgraph<-ggplot(data=AE_FR, aes(x=Size, y=FR, color=Size)) + theme_article()+
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Feeding rate") + ylab("Log(Feeding rate)") +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1, alpha=0.5)
FRgraph<- FRgraph + theme(legend.position = "none") 

FRgraph<-FRgraph+scale_color_manual(values=c("royalblue4","springgreen4"))

FRgraph_FINAL<-FRgraph + facet_grid(cols = vars(Concentration))+
  theme(strip.background = element_rect(colour="black", fill="white", size=0.5, linetype="solid"))+ 
  theme(panel.border = element_rect(fill = "transparent",color = "black", linewidth = 0.5)) 

FRgraph_FINAL

###They can then be merged together

GraphFR_AE<-plot_grid(FRgraph_FINAL, AEgraph_FINAL, ncol=2, nrow=1)
GraphFR_AE

