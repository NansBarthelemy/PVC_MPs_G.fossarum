library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(forcats)
library(egg)

setwd(dir="C:/Users/Nans/Desktop/PVC_MPs_G.fossarum/Data/(H1) Mortality")
data<- read.delim("Mortality_Gammarus.txt")
View(data)
data$Size <- fct_relevel(data$Size.Class, c("S", "L", "Control"))
data$Concentration<- fct_relevel(data$concentration, c("L", "H", "E","Control"))

msurv=Surv(data$Time,data$status)
msurv

#########################Graphical representation of the survival curve of all treatments

mfitTreatment=survfit(msurv~data$Treatment)
plot(mfitTreatment, col=c("darkorange", "darkorange3", "darkorange4", "plum", "orchid4", "purple2", "gold"))
legend("bottomleft", c("GLS", "GHS","GES","GLL","GHL","GEL","CG"),lty=1, col=c("darkorange", "darkorange3", "darkorange4", "plum", "orchid4", "purple2", "gold"))
survdiff(msurv~data$Treatment)

####Graphical representation and comparison of survival curve between the MP treatments and the control

#GLS vs CG

mfitLS=survfit(msurv~data$LS.control)
plot(mfitLS, col=c("red", "blue"))
legend("bottomleft", c("GLS","CG"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LS.control)

#GHS vs CG

mfitMS=survfit(msurv~data$HS.control)
plot(mfitMS, col=c("red", "blue"))
legend("bottomleft", c("GHS","CG"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$HS.control)

#GES vs CG

mfitHS=survfit(msurv~data$ES.control)
plot(mfitHS, col=c("red", "blue"))
legend("bottomleft", c("GES","CG"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$ES.control)

#GLL vs CG

mfitLL=survfit(msurv~data$LL.control)
plot(mfitLL, col=c("red", "blue"))
legend("bottomleft", c("GLL","CG"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LL.control)

#GHL vs CG

mfitML=survfit(msurv~data$HL.control)
plot(mfitML, col=c("red", "blue"))
legend("bottomleft", c("GHL","CG"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$HL.control)

#GEL vs CG

mfitHL=survfit(msurv~data$EL.control)
plot(mfitHL, col=c("red", "blue"))
legend("bottomleft", c("GEL","CG"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$EL.control)


#########################Statistical analysis and graphical representation focusing onlu on MP sizes

###i.e. all treatments with similar MP size were regrouped

#Global comparison and graphical representation

mfitSize=survfit(msurv~data$Size)
plot(mfitSize, col=c("grey1", "firebrick3", "royalblue1"))
legend("bottomleft", c("Small MPs (20-63µm)", "Large MPs (125-160µm)","CG"),lty=1, col=c("grey1", "firebrick3", "royalblue1"))
survdiff(msurv~data$Size)


####Small MPs vs CG 

mfitSCG=survfit(msurv~data$Small.CG)
plot(mfitSCG, col=c("royalblue1", "firebrick3"), xlab = "Day", ylab = c("Survival probability"))
legend("bottomleft", c("Small MPs (20-63µm)", "CG"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Small.CG)


####Large MPs vs CG 

mfitLCG=survfit(msurv~data$Large.CG)
plot(mfitLCG, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("Large MPs (125-160µm)", "CG"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Large.CG)

####Small MPs vs Large MPs

mfitSL=survfit(msurv~data$Small.Large)
plot(mfitSL, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("Small MPs (20-63µm)", "Large MPs (125-160µm)"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Small.Large)

#########################Statistical analysis and graphical representation focusing only on MP concentrations

###i.e. all treatments with similar MP concentration were regrouped

#Global comparison and graphical representation

mfitConcentration=survfit(msurv~data$Concentration)
plot(mfitConcentration, col=c("firebrick3", "green4", "darkorange4", "royalblue1"))
legend("bottomleft", c("Low (500 MPs/L)", "High (5000 MPs/L)","Extreme (50000 MPs/L)","CG"),lty=1, col=c("firebrick3", "green4", "darkorange4", "royalblue1"))
survdiff(msurv~data$Concentration)

###Low concentration vs CG 

mfitLowCG=survfit(msurv~data$Low.CG)
plot(mfitLowCG, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("Low (500 MPs/L)", "CG"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Low.CG)

###High concentration vs CG 

mfitHighCG=survfit(msurv~data$High.CG)
plot(mfitHighCG, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("High (5000 MPs/L)", "CG"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$High.CG)

###Extreme concentration vs CG 

mfitExtremeCG=survfit(msurv~data$Extreme.CG)
plot(mfitExtremeCG, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("Extreme (50000 MPs/L)", "CG"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Extreme.CG)

#Low concentration vs High concentration

mfitLowHigh=survfit(msurv~data$Low.High)
LowHigh<-plot(mfitLowHigh, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("Low (500 MPs/L)", "High (5000 MPs/L)"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Low.High)

#Low concentration vs Extreme concentration 

mfitLowExtreme=survfit(msurv~data$Low.Extreme)
LowExtreme<-plot(mfitLowExtreme, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("Low (500 MPs/L)", "Extreme (50000 MPs/L)"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$Low.Extreme)

#High concentration vs Extreme concentration 

mfitHighExtreme=survfit(msurv~data$High.Extreme)
HighExtreme<-plot(mfitHighExtreme, col=c("royalblue1", "firebrick3"))
legend("bottomleft", c("High (5000 MPs/L)", "Extreme (50000 MPs/L)"),lty=1, col=c("royalblue1", "firebrick3"))
survdiff(msurv~data$High.Extreme)




######################### Statistical analysis and graphical representation focusing on MP size and concentration

######Comparison between MP abundance within small MP treatments

#Low concentration vs High concentration (20-63µm)

mfitLS.MS=survfit(msurv~data$LS.HS)
plot(mfitLS.MS, col=c("red", "blue"))
legend("bottomleft", c("GLS","GHS"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LS.HS)

#Low concentration vs Extreme concentration (20-63µm)

mfitLSHS=survfit(msurv~data$LS.ES)
plot(mfitLSHS, col=c("red", "blue"))
legend("bottomleft", c("GLS","GES"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LS.ES)

#High concentration vs Extreme concentration (20-63µm)

mfitMLHL=survfit(msurv~data$HS.ES)
plot(mfitMLHL, col=c("red", "blue"))
legend("bottomleft", c("GHS","GES"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$HS.ES)

######Comparison between MP abundance within large MP treatments

#Low concentration vs High concentration (125-160µm)

mfitLLML=survfit(msurv~data$LL.HL)
plot(mfitLLML, col=c("red", "blue"))
legend("bottomleft", c("GLL","GHL"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LL.HL)

#Low concentration vs Extreme concentration (125-160µm)

mfitLLHL=survfit(msurv~data$LL.EL)
plot(mfitLLHL, col=c("red", "blue"))
legend("bottomleft", c("GLL","GEL"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LL.EL)

#High concentration vs Extreme concentration (125-160µm)

mfitMLHL=survfit(msurv~data$HL.EL)
plot(mfitMLHL, col=c("red", "blue"))
legend("bottomleft", c("GHL","GEL"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$HL.EL)


######Comparison between MP size within the abundance treatments


#at Low concentration : Small MPs (20-63µm) vs Large MPs (125-160µm)

mfitLSLL=survfit(msurv~data$LS.LL)
plot(mfitLSLL, col=c("red", "blue"))
legend("bottomleft", c("GLS","GLL"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$LS.LL)

#at High concentration : Small MPs (20-63µm) vs Large MPs (125-160µm)

mfitMSML=survfit(msurv~data$HS.HL)
plot(mfitMSML, col=c("red", "blue"))
legend("bottomleft", c("GHS","GHL"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$HS.HL)

#at Extreme concentration : Small MPs (20-63µm) vs Large MPs (125-160µm)

mfitHSHL=survfit(msurv~data$ES.EL)
plot(mfitHSHL, col=c("red", "blue"))
legend("bottomleft", c("GES","GEL"),lty=1,col=c("red", "blue"))
survdiff(msurv~data$ES.EL)


#############Final graphical representations based on the different MP concentrations within both size range

dataSMALL<- read.delim("Mortality_Gammarus_Small.txt")
dataLARGE<- read.delim("Mortality_Gammarus_Large.txt")
View(dataSMALL)
View(dataLARGE)

#####Different concentrations of small MPs (20-63µm)

msurv=Surv(dataSMALL$Time,dataSMALL$status)
msurv

fit_Small <-survfit(Surv(Time,status)~Treatment, data=dataSMALL)
Graph_Small<-ggsurvplot(fit_Small, palette = c ("#009E73", "#0072B2", "#D55E00", "#000000"), data = data, risk.table = FALSE, xlim = c(0,28),
                        break.time.by = 7,legend.labs =c("GLS (BC)", "GHS (B)","GES (C)","CG (A)"),font.legend = c(14,"plain"), title = "Small MPs (20-63µm)", ggtheme = theme_article()+
                          theme(plot.title = element_text(hjust = 0.5)), legend = "bottom",
                        legend.title = "")
Graph_Small 

#####Different concentrations of Large MPs (125-160µm)

msurv=Surv(dataLARGE$Time,dataLARGE$status)
msurv

fit_Large <-survfit(Surv(Time,status)~Treatment, data=dataLARGE)
Graph_Large<-ggsurvplot(fit_Large, palette = c ("#009E73", "#0072B2", "#D55E00", "#000000"), data = data, risk.table = FALSE, xlim = c(0,28),
                        break.time.by = 7,legend.labs =c("GLL (B)", "GHL (B)","GEL (AB)","CG (A)"),font.legend = c(14,"plain"), title = "Large MPs (125-160µm)", ggtheme = theme_article()+
                          theme(plot.title = element_text(hjust = 0.5)), legend = "bottom",
                        legend.title = "")
Graph_Large 
Graph_Large 











