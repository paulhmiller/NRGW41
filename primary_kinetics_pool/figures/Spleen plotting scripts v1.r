
library(reshape2)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)


source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")

############################ Spleen plots ##################################

#read in data from each experiment:
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRG_W41/primary_kinetics_pool/figures")
dat<-read.csv("QR spleen kinetics.csv", colClasses=c(Sex="factor"))

#change directory to pooled_data directory
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRG_W41/primary_kinetics_pool/figures")
dir.create("Spleen")
setwd("./Spleen")

dat<-dat[,c(1:5,8:22)]  #Take just the necessary columns

#calculate absolute counts with bead numbers
factor <- (5*983)/(dat[,6]+dat[,7])
tmp <- dat[,8:14]*factor
tmp[1:6,] <- tmp[1:6,]*10 #Multiply R mice by 10 (because FACed just 10%) Both had 5uL beads.
tmp[11:24,] <- tmp[11:24,]*4 #Multiply R mice by 10 (because FACed just 10%) Both had 5uL beads.
dat <- dat[,c(1:5,15:20)]
dat <- cbind(dat,tmp)

dat[,6:11][dat[,6:11] < 0.01] <- 0.01  #change low values to detection threshold
dat[,12:18][dat[,12:18] < 10] <- 10  #change low values to detection threshold
datlt <- dat  #log10 transform
datlt[6:ncol(dat)] <- log10(dat[6:ncol(dat)])

tall<-melt(datlt,id.vars=c("Exp","Week","Strain","Sex","Irradiation.Dose"),na.rm=TRUE)  #id.vars refer to the things you want to keep linked
#merge male and female
#levels(tall$Sex)[levels(tall$Sex)=="M"] <- "both"
#levels(tall$Sex)[levels(tall$Sex)=="F"] <- "both"
names(tall)[names(tall) == 'variable'] <- 'Lineage'	

#Change/merge 150 and 900 both to irradiated, so that there are only 2 groups (unirradiated and irradiated).
levels(tall$Irradiation.Dose)[levels(tall$Irradiation.Dose)=="0 Rad"] <- "Non-irradiated"
levels(tall$Irradiation.Dose)[levels(tall$Irradiation.Dose)=="150 Rad"] <- "Irradiated" 
levels(tall$Irradiation.Dose)[levels(tall$Irradiation.Dose)=="900 Rad"] <- "Irradiated"

##** toggles
#tall <- tall[tall$Irradiation.Dose=="Irradiated",]

datst <- ddply(tall, c("Week","Strain","Sex","Irradiation.Dose","Lineage"), summarise,
               N    = sum(!is.na(value)),
               mean = mean(value, na.rm=TRUE),
               sd   = sd(value, na.rm=TRUE),
               se   = sd / sqrt(N)
)

# THIS STOPPED WORKING. USE ddply INSTEAD. datst <- summarySE(tall, measurevar="value", groupvars=c("Week","Strain","Sex","Irradiation.Dose","Lineage"),na.rm=TRUE)   # summarySE provides the std dev, SEM, and (default 95%) confidence interval. #measurevar is the x-axis


tp<-datst[datst$Week=='30',]   #Take a single timepoint


#subset by lineage so that each can be plotted separately 
CD45<-datst[datst$Lineage == "CD45",]
GPA<-datst[datst$Lineage == "GPA",]
GM<-datst[datst$Lineage == "CD33.15",]
B<-datst[datst$Lineage == "CD19",]
T<-datst[datst$Lineage == "CD3",]
CD34<-datst[datst$Lineage == "CD34",]



plot(CD45$mean~CD45$Sex)
plot(CD45$mean~CD45$Irradiation.Dose)
plot(CD45$mean~CD45$Strain)
















#########  BM kinetics line plot ###########

plot<-ggplot(CD45, aes(x=Week,y=value,colour=Strain,shape=Sex,lty=Irradiation.Dose)) + 
  ggtitle("CD45") +
  coord_cartesian(xlim=c(1,31),ylim=c(-1,2)) +
#  scale_y_continuous(breaks=c(-2,-1,0,1,2,3,4),labels = c(0.01, 0.1, 1, 10, 100, 1000, 10000)) + 
  scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3, 4, 5, 6), labels = c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) + 
  #	geom_hline(yintercept=-0.69897, linetype="dotted", size=1) +	
  scale_colour_manual(values=paletteA) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se, lty=Irradiation.Dose),  width=0.5, size=0.75) +	
  geom_line(size=1) +		
  scale_linetype_manual(name="Strain", values=c("Non-irradiated"='dotted',"Irradiated"='solid')) + 
  #    scale_colour_hue(l=45) +
  geom_point(size = 3) +
  #    scale_shape_manual(values=c("Non-irradiated"=15,"Irradiated"=19)) +
  xlab("weeks post-transplant") +
  ylab("% total BM cells") +
  scale_x_continuous(breaks=c(3,6,10,20, 30)) +
  annotation_logticks(sides = "l") +
  #facet_wrap( ~Irradiation.Dose, ncol=2) +
#  guides(lty=F) +  #(colour=F,lty=F,shape=F) +
  themePM1()
plot
ggsave(filename="./BM/AMR CB20K38K BM CD45 kinetics.tiff",width=12,height=12,units="cm",dpi=600, compression="lzw")


#plot +
#	coord_cartesian(xlim=c(1,31),ylim=c(-1,2)) 
#ggsave(filename="150606 BM100K BM CD45 axis01 kinetics.tiff",width=8,height=8.1,units="cm",dpi=600, compression="lzw")

plot %+% GPA +
  ggtitle("GPA") 
ggsave(filename="./BM/AMR CB20K38K BM GPA kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% GM +
  ggtitle("CD15/33") 
ggsave(filename="./BM/AMR CB38K BM GM kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")
#ggsave(filename="150505 BM100K BM GM kinetics.wmf",width=12,height=12,units="cm",dpi=1200)

plot %+% B +
  ggtitle("CD19") 
ggsave(filename="./BM/AMR CB38K BM CD19 kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% T +
	ggtitle("CD3") 
ggsave(filename="./BM/AMR CB38K BM CD3 kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% CD34 +
  ggtitle("CD34") 
ggsave(filename="./BM/AMR CB38K BM CD34 kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")





#########  BM bar plots ###########

#tmp <- datst
datst <- tmp


#select timepoints to plot
datst <- datst[datst$Week==30,]


#subset by lineage so that each can be plotted separately 
CD45<-datst[datst$Lineage == "CD45+",]
GPA<-datst[datst$Lineage == "GPA+",]
GM<-datst[datst$Lineage == "CD33+15+",]
B<-datst[datst$Lineage == "CD19+",]
T<-datst[datst$Lineage == "CD3+",]
CD34<-datst[datst$Lineage == "CD34+",]



plot<-ggplot(CD45,aes(x=Strain, y=value, fill=Sex)) + 
#  scale_fill_manual(values=paletteA) +
  geom_bar(stat="identity", position="dodge", width=0.8, colour="black") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), position=position_dodge(width=0.8), width=0.4) +
  ylab("% total BM cells") +
  ggtitle("CD45") +
  coord_cartesian(ylim=c(0,2)) +
  scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3, 4, 5, 6), labels = c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) + 
#  scale_y_continuous(breaks=c(-1,0,1,2),labels = c(0.1, 1, 10, 100)) + 
#  #offset by 3 logs:  scale_y_continuous(breaks=c(1,2,3,4,5),labels = c(0.01, 0.1, 1, 10, 100)) + #offset by 3 logs
#  scale_x_discrete("", labels=c("Irradiated","Non-irradiated")) +
  annotation_logticks(sides = "l") +
#  guides(fill=FALSE,lty=FALSE) +
  themePM1()  
plot
ggsave(filename="./BM/AMR CB20K38K BM CD45 w30 bar.tiff",width=10,height=10,units="cm",dpi=600, compression="lzw")

























############################ CB20K - BM kinetic line plots ##################################


#read in data from each experiment:
setwd("C:/Users/paulm/CRC Paul/Experiments/")
Qdat<-read.csv("./2014/140825Q W41 txpt4/figures/QCB20K BM kinetics.csv", colClasses=c("Sex"="factor"))
Rdat<-read.csv("./2014/141118R W41 again/figures/R_BM kinetics.csv")

#change directory to pooled_data directory
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRG_W41/primary_kinetics_pool/figures")

#add exp label column
Qdat<-cbind("140825Q",Qdat)
Rdat<-cbind("141118R",Rdat)
names(Qdat)[1] <- "Exp"
names(Rdat)[1] <- "Exp"

#delete non-common columns
Qdat<-Qdat[,c(1:13)]
Rdat<-Rdat[,c(1:13)]

#merge
dat <- rbind(Qdat, Rdat)

#just take CB20K
dat<-dat[dat$Input=="CB.20K",] 

#write a csv file for record keeping
#write.table(dat, file= "QR CB20K BM kinetics.csv", sep=",", row.names=F)


dat<-dat[,c(1:5,8:13)]  #Take just the necessary columns


names(dat)<-c("Exp","Week","Strain","Sex","Irradiation.Dose","CD45+","GPA+","CD33+15+","CD19+","CD3+","CD34+")  #change column names so it looks better on plots
dat[,6:11][dat[,6:11] < 0.01] <- 0.01  #change low values in leukocyte columns to detection threshold
datlt <- dat  #log10 transform
datlt[6:11] <- log10(dat[6:11])
tall<-melt(datlt,id.vars=c("Exp","Week","Strain","Sex","Irradiation.Dose"),na.rm=TRUE)  #id.vars refer to the things you want to keep linked
#merge male and female
#levels(tall$Sex)[levels(tall$Sex)=="M"] <- "both"
#levels(tall$Sex)[levels(tall$Sex)=="F"] <- "both"
names(tall)[names(tall) == 'variable'] <- 'Lineage'	
datst <- summarySE(tall, measurevar="value", groupvars=c("Week","Strain","Sex","Irradiation.Dose","Lineage"),na.rm=TRUE)   # summarySE provides the std dev, SEM, and (default 95%) confidence interval. #measurevar is the x-axis
#tp<-datst[datst$Week=='10',]   #Take a single timepoint

#remove NRG 150 Rad
#datst <- datst[!(datst$Irradiation.Dose=="150 Rad" & datst$Strain=="NRG"),]

#Change/merge 150 and 900 both to irradiated, so that there are only 2 groups (unirradiated and irradiated).
levels(datst$Irradiation.Dose)[levels(datst$Irradiation.Dose)=="0 Rad"] <- "Non-irradiated"
levels(datst$Irradiation.Dose)[levels(datst$Irradiation.Dose)=="150 Rad"] <- "Irradiated" 
levels(datst$Irradiation.Dose)[levels(datst$Irradiation.Dose)=="900 Rad"] <- "Irradiated"

##** toggles
#datst <- datst[datst$Irradiation.Dose=="Irradiated",]


#subset by lineage so that each can be plotted separately 
CD45<-datst[datst$Lineage == "CD45+",]
GPA<-datst[datst$Lineage == "GPA+",]
GM<-datst[datst$Lineage == "CD33+15+",]
B<-datst[datst$Lineage == "CD19+",]
T<-datst[datst$Lineage == "CD3+",]
CD34<-datst[datst$Lineage == "CD34+",]


#########  BM kinetics line plot ###########

plot<-ggplot(CD45, aes(x=Week,y=value,colour=Strain,shape=Sex,lty=Irradiation.Dose)) + 
  ggtitle("CD45") +
  coord_cartesian(xlim=c(1,21),ylim=c(-1,2)) +
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3,4),labels = c(0.01, 0.1, 1, 10, 100, 1000, 10000)) + 
  #	geom_hline(yintercept=-0.69897, linetype="dotted", size=1) +	
  scale_colour_manual(values=paletteA) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se, lty=Irradiation.Dose),  width=0.5, size=0.75) +	
  geom_line(size=1) +		
  scale_linetype_manual(name="Strain", values=c("Non-irradiated"='dotted',"Irradiated"='solid')) + 
  #    scale_colour_hue(l=45) +
  geom_point(size = 3) +
  #    scale_shape_manual(values=c("Non-irradiated"=15,"Irradiated"=19)) +
  xlab("weeks post-transplant") +
  ylab("percent of total BM") +
  scale_x_continuous(breaks=c(3,6,10,20, 30)) +
  annotation_logticks(sides = "l") +
  #facet_wrap( ~Irradiation.Dose, ncol=2) +
  #  guides(colour=F,lty=F,shape=F) +
  themePM1()
plot
ggsave(filename="QR CB20K BM CD45 kinetics irrad.tiff",width=12,height=12,units="cm",dpi=600, compression="lzw")


#plot +
#	coord_cartesian(xlim=c(1,31),ylim=c(-1,2)) 
#ggsave(filename="150606 BM100K BM CD45 axis01 kinetics.tiff",width=8,height=8.1,units="cm",dpi=600, compression="lzw")

plot %+% GPA +
  ggtitle("GPA") 
ggsave(filename="QR CB20K BM GPA kinetics irrad.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% GM +
  ggtitle("CD15/33") 
ggsave(filename="QR CB20K BM GM kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")
#ggsave(filename="150505 BM100K BM GM kinetics.wmf",width=8,height=8.1,units="cm",dpi=1200)

plot %+% B +
  ggtitle("CD19") 
ggsave(filename="QR CB20K BM B kinetics.tiff",width=8,height=8.1,units='cm',dpi=600,compression="lzw")

#plot %+% T +
#	ggtitle("CD3") 
#ggsave(filename="QR CB20K BM T kinetics.tiff",width=8,height=8.1,units='cm',dpi=600,compression="lzw")

plot %+% CD34 +
  ggtitle("CD34") 
ggsave(filename="QR CB20K BM CD34 kinetics.tiff",width=8,height=8.1,units='cm',dpi=600,compression="lzw")




