library(ggplot2)
library(cowplot)
library(dplyr)
library(plyr)
library(egg)
library(car)
library(forcats)

setwd(dir="C:/Users/Nans/Desktop/PVC_MPs_G.fossarum/Data/Supplementary materials")

Size_Test<-read.csv("Size_Test.csv", header = T, sep = ";")  
Size_Test$Size.class <- fct_relevel(Size_Test$Size.class, c("20-63", "63-80","80-100","125-160"))
View(Size_Test)


#Comparison of the average MP size recovered in each size class


shapiro.test(Size_Test$Mean)
bartlett.test(Mean~Size.class, data=Size_Test)
pairwise.wilcox.test(Size_Test$Mean,Size_Test$Size.class, p.adj="bonferroni", correct = FALSE)

Mean_MPs_Size<-ggplot(data=Size_Test, aes(x=Size.class, y=Mean)) +
  ggtitle("Mean particles size") + 
  xlab("Size classes (µm)") + ylab("PVC particle size")+ theme_article() +
  geom_boxplot(notch = FALSE,fill="gray92", color="gray23")+ theme(plot.title = element_text(hjust = 0.5))

Mean_MPs_Size



#Graphical representation and comparison of the abundance of MPs recovered in each size class


SizeClassA<-subset(Size_Test, Size.class=="20-63")
shapiro.test(SizeClassA$Abundance)
#Since there are only 6 replicates per Size Class, we used a kruskal.wallis test even if shapiro.test is >0.05
kruskal.test(Abundance ~ Replicate, data = SizeClassA)

AbundanceA<-ggplot(data=SizeClassA, aes(x=Replicate, y=Abundance)) +
  ggtitle("PVC pipe 20-63µm") + 
  xlab("Replicates") + ylab("Abundance")+ theme_minimal() +
  geom_bar(stat="identity",fill="gray92", color="gray23", width=0.5)+
  geom_text(aes(label=Abundance), vjust=1.6, color="black", size=3.5)+
  theme(plot.title = element_text(hjust = 0.5))
AbundanceA  


SizeClassB<-subset(Size_Test, Size.class=="63-80")
shapiro.test(SizeClassB$Abundance)
#Since there are only 6 replicates per Size Class, we used a kruskal.wallis test even if shapiro.test is >0.05
kruskal.test(Abundance ~ Replicate, data = SizeClassB)

AbundanceB<-ggplot(data=SizeClassB, aes(x=Replicate, y=Abundance)) +
  ggtitle("PVC pipe 63-80µm") + 
  xlab("Replicates") + ylab("Abundance")+ theme_minimal() +
  geom_bar(stat="identity",fill="gray92", color="gray23", width=0.5)+
  geom_text(aes(label=Abundance), vjust=1.6, color="black", size=3.5)+
  theme(plot.title = element_text(hjust = 0.5))
AbundanceB

SizeClassC<-subset(Size_Test, Size.class=="80-100")
shapiro.test(SizeClassC$Abundance)
#Since there are only 6 replicates per Size Class, we used a kruskal.wallis test even if shapiro.test is >0.05
kruskal.test(Abundance ~ Replicate, data = SizeClassC)

AbundanceC<-ggplot(data=SizeClassC, aes(x=Replicate, y=Abundance)) +
  ggtitle("PVC pipe 80-100µm") + 
  xlab("Replicates") + ylab("Abundance")+ theme_minimal() +
  geom_bar(stat="identity",fill="gray92", color="gray23", width=0.5)+
  geom_text(aes(label=Abundance), vjust=1.6, color="black", size=3.5)+
  theme(plot.title = element_text(hjust = 0.5))
AbundanceC


SizeClassD<-subset(Size_Test, Size.class=="125-160")
shapiro.test(SizeClassD$Abundance)
#Since there are only 6 replicates per Size Class, we used a kruskal.wallis test even if shapiro.test is >0.05
kruskal.test(Abundance ~ Replicate, data = SizeClassD)

AbundanceD<-ggplot(data=SizeClassD, aes(x=Replicate, y=Abundance)) +
  ggtitle("PVC pipe 125-160µm") + 
  xlab("Replicates") + ylab("Abundance")+ theme_minimal() +
  geom_bar(stat="identity",fill="gray92", color="gray23", width=0.5)+
  geom_text(aes(label=Abundance), vjust=1.6, color="black", size=3.5)+
  theme(plot.title = element_text(hjust = 0.5))
AbundanceD

Abundance_Pipe<-plot_grid(AbundanceA, AbundanceB ,AbundanceC ,AbundanceD, ncol = 2, nrow = 2)
Abundance_Pipe

