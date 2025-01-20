library(dplyr)
library(agricolae)
library(cowplot)
library(egg)
library(car)
library(ggplot2)
library(ggsignif)
library(vegan)
library(forcats)

setwd(dir="C:/Users/Nans/Desktop/PVC_MPs_G.fossarum/Data/(H2) Proteomic")


###################Global peptides comparison using PERMANOVA (only done on peptides with 100% quantification rate)

Prot100<-read.csv("Peptides_100%_Quantification.csv", header = T, sep = ";")
View(Prot100)
Peptides = Prot100[,6:ncol(Prot100)]
View(Peptides)
Peptides.env <- Prot100[,-13]
Peptides.env <- Peptides.env[,-12]
Peptides.env <- Peptides.env[,-11]
Peptides.env <- Peptides.env[,-10]
Peptides.env <- Peptides.env[,-9]
Peptides.env <- Peptides.env[,-8]
Peptides.env <- Peptides.env[,-7]
Peptides.env <- Peptides.env[,-6]
View(Peptides.env)

adonis2(Peptides~Treatment, data = Peptides.env, permutations = 999, method="bray")
adonis2(Peptides ~ Size, data = Prot100, permutations = 999, method="bray")
adonis2(Peptides ~ Concentration, data = Prot100, permutations = 999, method="bray")

Peptides.dist <- vegdist(Peptides, method="bray")
dispersion <- betadisper(Peptides.dist, group=Peptides.env$Treatment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) 

#No effects of the different treatments, MP size and concentration on the overall peptides expression rate

#We rerun the analysis on subsets of data representing each treatment

CG_GHL<-subset(Prot100, Treatment=="CG"|Treatment=="GHL")
View(CG_GHL)
CG_GLL<-subset(Prot100, Treatment=="CG"|Treatment=="GLL")
View(CG_GLL)
CG_GEL<-subset(Prot100, Treatment=="CG"|Treatment=="GEL")
View(CG_GEL)
CG_GHS<-subset(Prot100, Treatment=="CG"|Treatment=="GHS")
View(CG_GHS)
CG_GLS<-subset(Prot100, Treatment=="CG"|Treatment=="GLS")
View(CG_GLS)
CG_GES<-subset(Prot100, Treatment=="CG"|Treatment=="GES")
View(CG_GES)

adonis2(Peptides~Prot100$Treatment, data = CG_GHL, permutations = 999, method="bray")
adonis2(Peptides~Prot100$Treatment, data = CG_GLL, permutations = 999, method="bray")
adonis2(Peptides~Prot100$Treatment, data = CG_GEL, permutations = 999, method="bray")
adonis2(Peptides~Prot100$Treatment, data = CG_GHS, permutations = 999, method="bray")
adonis2(Peptides~Prot100$Treatment, data = CG_GLS, permutations = 999, method="bray")
adonis2(Peptides~Prot100$Treatment, data = CG_GES, permutations = 999, method="bray")

######Graphical representation using NMDS

Peptides = Prot100[,6:ncol(Prot100)]
View(Peptides)
Peptides.env <- Prot100[,1]
View(Peptides.env)
ord<-metaMDS(Peptides, distance="bray")
ord
ord$stress
plot(ord, type="t")
ordihull(ord,groups=Peptides.env,draw="polygon",col="grey90",label=T)
orditorp(ord,display="species",col="red",air=0.01)

###################Statistical analysis and graphical representation of each individual peptide with 100% quantification rate

Prot100<-read.csv("Peptides_100%_Quantification.csv", header = T, sep = ";")
Prot100$Treatment <-fct_relevel(Prot$Treatment, c("CG", "GLS","GHS","GES","GLL","GHL","GEL"))
View(Prot100)


####Comparison of peptide concentrations between the different PVC-MP exposure treatment

#ADP

shapiro.test(Prot$ADP)
bartlett.test(ADP~Treatment, data=Prot)
res.aov.ADP <- aov(ADP ~ Treatment, data = Prot)
summary(res.aov.ADP)
TukeyHSD(res.aov.ADP, conf.level=.95)
ADP_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=ADP, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("ADPALGQAIQER concentration (pmol/mg)")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
ADP_Concentration<-ADP_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
ADP_Concentration

#EFW

shapiro.test(Prot$EFW)
bartlett.test(EFW~Treatment, data=Prot)
res.aov.EFW <- aov(EFW ~ Treatment, data = Prot)
summary(res.aov.EFW)
TukeyHSD(res.aov.EFW, conf.level=.95)
EFW_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=EFW, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("EFWIATDHNEVR concentration (pmol/mg)")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
EFW_Concentration<-EFW_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
EFW_Concentration

#FVG

shapiro.test(Prot$FVG)
bartlett.test(FVG~Treatment, data=Prot)
res.aov.FVG <- aov(FVG ~ Treatment, data = Prot)
summary(res.aov.FVG)
TukeyHSD(res.aov.FVG, conf.level=.95)
FVG_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=FVG, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("FVGLISLIDPPR concentration (pmol/mg)")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
FVG_Concentration<-FVG_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
FVG_Concentration

#LIF

shapiro.test(Prot$LIF)
bartlett.test(LIF~Treatment, data=Prot)
res.aov.LIF <- aov(LIF ~ Treatment, data = Prot)
summary(res.aov.LIF)
TukeyHSD(res.aov.LIF, conf.level=.95)
LIF_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=LIF, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LIFDNLK concentration (pmol/mg)")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
LIF_Concentration<-LIF_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
LIF_Concentration

#LII

shapiro.test(Prot$LII)
bartlett.test(LII~Treatment, data=Prot)
kruskal(Prot$LII,Prot$Treatment, group=TRUE)$groups
LII_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=LII, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LIIVEGCQR concentration (pmol/mg)")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
LII_Concentration<-LII_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
LII_Concentration

#LVL

shapiro.test(Prot$LVL)
bartlett.test(LVL~Treatment, data=Prot)
res.aov.LVL <- aov(LVL ~ Treatment, data = Prot)
summary(res.aov.LVL)
TukeyHSD(res.aov.LVL, conf.level=.95)
LVL_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=LVL, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LVLGTATYGR concentration (pmol/mg)")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
LVL_Concentration<-LVL_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
LVL_Concentration

#SGQ 

shapiro.test(Prot$SGQ)
bartlett.test(SGQ~Treatment, data=Prot)
res.aov.SGQ <- aov(SGQ ~ Treatment, data = Prot)
summary(res.aov.SGQ)
TukeyHSD(res.aov.SGQ, conf.level=.95)
SGQ_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=SGQ, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("SGQDGVPILK concentration (pmol/mg")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
SGQ_Concentration <-SGQ_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
SGQ_Concentration

#VIM

shapiro.test(Prot$VIM)
bartlett.test(VIM~Treatment, data=Prot)
res.aov.VIM <- aov(VIM ~ Treatment, data = Prot)
summary(res.aov.VIM)
TukeyHSD(res.aov.VIM, conf.level=.95)
VIM_Concentration<-ggplot(data=Prot, aes(x=Treatment, y=VIM, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("VIMVTGDHPITAK concentration (pmol/mg")+ theme_article() +
  geom_boxplot(notch=FALSE, outliers=FALSE)+
  geom_jitter(color="black", size=1.5, alpha=0.5)
VIM_Concentration <-VIM_Concentration+ scale_color_manual(values=c("gray4", "springgreen4","royalblue4"))
VIM_Concentration

#The graphics can then be merged together (2x2)

Plot_prot_1<-plot_grid(FVG_Concentration, VIM_Concentration,SGQ_Concentration,LIF_Concentration, ncol=2, nrow=2)
Plot_prot_1

Plot_prot_2<-plot_grid(LII_Concentration, EFW_Concentration,LVL_Concentration,ADP_Concentration, ncol=2, nrow=2)
Plot_prot_2





###################Statistical analysis and graphical representation of each individual peptide with 50-99% quantification rate

###ELF

#Concentration

shapiro.test(Prot50$ELF)
bartlett.test(ELF~Concentration, data=Prot50)
kruskal(Prot50$ELF,Prot50$Concentration, group=TRUE)$groups
ELF_Concentration<-ggplot(data=Prot50, aes(x=Concentration, y=ELF, color=Concentration))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("ELFDFADAHR (Cellulase)")+ theme_article() +
  geom_point(aes(color=Concentration))
ELF_Concentration<-ELF_Concentration+ scale_color_manual(values=c("grey4", "#009E73", "#0072B2", "#D55E00"))
ELF_Concentration

#Size

shapiro.test(Prot50$ELF)
bartlett.test(ELF~Size, data=Prot50)
res.aov.ELF <- aov(ELF ~ Size, data = Prot50)
summary(res.aov.ELF)
TukeyHSD(res.aov.ELF, conf.level=.95)
ELF_Size<-ggplot(data=Prot50, aes(x=Size, y=ELF, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("ELFDFADAHR (Cellulase)")+ theme_article() +
  geom_point(aes(color=Size))
ELF_Size<-ELF_Size+ scale_color_manual(values=c("gray4", "royalblue4","springgreen4"))
ELF_Size


###EVN

#Concentration

shapiro.test(Prot50$EVN)
bartlett.test(EVN~Concentration, data=Prot50)
kruskal(Prot50$EVN,Prot50$Concentration, group=TRUE)$groups
EVN_Concentration<-ggplot(data=Prot50, aes(x=Concentration, y=EVN, color=Concentration))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("EVNGDASEAALLK (Na+/K+ ATPase)")+ theme_article() +
  geom_point(aes(color=Concentration))
EVN_Concentration<-EVN_Concentration+ scale_color_manual(values=c("grey4", "#009E73", "#0072B2", "#D55E00"))
EVN_Concentration

#Size

shapiro.test(Prot50$EVN)
bartlett.test(EVN~Size, data=Prot50)
kruskal(Prot50$EVN,Prot50$Size, group=TRUE)$groups
EVN_Size<-ggplot(data=Prot50, aes(x=Size, y=EVN, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("EVNGDASEAALLK (Na+/K+ ATPase)")+ theme_article() +
  geom_point(aes(color=Size))
EVN_Size<-EVN_Size+ scale_color_manual(values=c("gray4", "royalblue4","springgreen4"))
EVN_Size

##GTL

#Concentration

shapiro.test(Prot50$GTL)
bartlett.test(GTL~Concentration, data=Prot50)
kruskal(Prot50$GTL,Prot50$Concentration, group=TRUE)$groups
GTL_Concentration<-ggplot(data=Prot50, aes(x=Concentration, y=GTL, color=Concentration))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("GTLAVIPVQNR (Transglutaminase)")+ theme_article() +
  geom_point(aes(color=Concentration))
GTL_Concentration<-GTL_Concentration+ scale_color_manual(values=c("grey4", "#009E73", "#0072B2", "#D55E00"))
GTL_Concentration

#Size

shapiro.test(Prot50$GTL)
bartlett.test(GTL~Size, data=Prot50)
kruskal(Prot50$GTL,Prot50$Size, group=TRUE)$groups
GTL_Size<-ggplot(data=Prot50, aes(x=Size, y=GTL, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("GTLAVIPVQNR (Transglutaminase)")+ theme_article() +
  geom_point(aes(color=Size))
GTL_Size<-GTL_Size+ scale_color_manual(values=c("gray4", "royalblue4","springgreen4"))
GTL_Size


##LAD

#Concentration

shapiro.test(Prot50$LAD)
bartlett.test(LAD~Concentration, data=Prot50)
kruskal(Prot50$LAD,Prot50$Concentration, group=TRUE)$groups
LAD_Concentration<-ggplot(data=Prot50, aes(x=Concentration, y=LAD, color=Concentration))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LADNIAGHVINTQEFIR (Catalase)")+ theme_article() +
  geom_point(aes(color=Concentration))
LAD_Concentration<-LAD_Concentration+ scale_color_manual(values=c("grey4", "#009E73", "#0072B2", "#D55E00"))
LAD_Concentration

#Size

shapiro.test(Prot50$LAD)
bartlett.test(LAD~Size, data=Prot50)
kruskal(Prot50$LAD,Prot50$Size, group=TRUE)$groups
LAD_Size<-ggplot(data=Prot50, aes(x=Size, y=LAD, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LADNIAGHVINTQEFIR (Catalase)")+ theme_article() +
  geom_point(aes(color=Size))
LAD_Size<-LAD_Size+ scale_color_manual(values=c("gray4", "royalblue4","springgreen4"))
LAD_Size


##LGS

#Concentration

shapiro.test(Prot50$LGS)
bartlett.test(LGS~Concentration, data=Prot50)
res.aov.LGS<- aov(LGS ~ Concentration, data = Prot50)
summary(res.aov.LGS)
TukeyHSD(res.aov.LGS, conf.level=.95)
LGS_Concentration<-ggplot(data=Prot50, aes(x=Concentration, y=LGS, color=Concentration))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LGSNFLQIPVNCPYR (Catalase)")+ theme_article() +
  geom_point(aes(color=Concentration))
LGS_Concentration<-LGS_Concentration+ scale_color_manual(values=c("grey4", "#009E73", "#0072B2", "#D55E00"))
LGS_Concentration

#Size

shapiro.test(Prot50$LGS)
bartlett.test(LGS~Size, data=Prot50)
shapiro.test(Prot50$LGS)
bartlett.test(LGS~Size, data=Prot50)
LGS_Size<-ggplot(data=Prot50, aes(x=Size, y=LGS, color=Size))+ theme(plot.title = element_text(hjust = 2)) +
  ylab("LGSNFLQIPVNCPYR (Catalase)")+ theme_article() +
  geom_point(aes(color=Size))
LGS_Size<-LGS_Size+ scale_color_manual(values=c("gray4", "royalblue4","springgreen4"))
LGS_Size

###The graphics can then be merged together (2 by 2)

Plot_prot_50_ELF_EVN<-plot_grid(ELF_Concentration,ELF_Size,EVN_Concentration,EVN_Size, ncol=2, nrow=2)
Plot_prot_50_ELF_EVN

Plot_prot_50_GTL_LAD<-plot_grid(GTL_Concentration,GTL_Size,LAD_Concentration,LAD_Size, ncol=2, nrow=2)
Plot_prot_50_GTL_LAD

Plot_prot_50_LGS<-plot_grid(LGS_Concentration,LGS_Size,LGS_Concentration,LGS_Size, ncol=2, nrow=2)
Plot_prot_50_LGS

