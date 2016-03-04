
library(reshape2)
library(ggplot2)
library(grid)
library(plyr)


source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")


# 151109: modified to include CB20K (Qdat, etc)



############################ CB38K - PB kinetic line plots ##################################



#read in data from each experiment:
setwd("C:/Users/paulm/CRC Paul/Experiments/2014")
Adat<-read.csv("./140310A B6NRGW41 txpt2/figures/PB kinetics.csv")
Mdat<-read.csv("./140708M B6NRGW41 txpt3 & teratoma into W41/figures/MQCB38K PB Kinetics.csv")
Rdat<-read.csv("./141118R W41 again/figures/R_PB kinetics.csv")
Qdat<-read.csv("./140825Q W41 txpt4/figures/QCB20K PB Kinetics.csv", colClasses=c(Sex="factor"))


#change directory to pooled_data directory
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRG_W41/primary_kinetics_pool/figures")
dir.create("PB")

#delete non-common and un-needed columns
Adat<-Adat[,c(1:5,10:14)]
Mdat<-Mdat[,c(1:5,10,12:14,16)]
Rdat<-Rdat[,c(1:5,10,12:14,16)]
Qdat<-Qdat[,c(1:5,10,12:14,16)]

#add exp label column
Adat<-cbind("140310A",Adat)
Mdat<-cbind("140708M",Mdat)
Rdat<-cbind("141118R",Rdat)
Qdat<-cbind("140825Q",Qdat)

names(Adat)[1] <- "Exp"
names(Mdat)[1] <- "Exp"
names(Rdat)[1] <- "Exp"
names(Qdat)[1] <- "Exp"


#merge
dat <- rbind(Adat, Mdat, Rdat, Qdat)

#just take CB38K
#dat<-dat[dat$Input!="CB.20K",] 

#write a csv file for record keeping
write.table(dat, file= "AMR CB20K40K PB kinetics.csv", sep=",", row.names=F)

#fix time-points so that data can be pooled:
dat <- dat[!(dat$Week==26),] #remove week 26

dat<-dat[,c(1:5,7:11)]  #Take just the necessary columns

names(dat)<-c("Exp","Week","Strain","Sex","Irradiation.Dose","CD45+","CD33+15+","CD19+","CD3+","Platelets")  #change column names so it looks better on plots
dat[,6:9][dat[,6:9] < 0.2] <- 0.2  #change low values in leukocyte columns to detection threshold
dat[,10][dat[,10] < 25] <- 25  #change low values in platelet column to detection threshold
datlt <- dat  #log10 transform
datlt[6:10] <- log10(dat[6:10])
tall<-melt(datlt,id.vars=c("Exp","Week","Strain","Sex","Irradiation.Dose"),na.rm=TRUE)  #id.vars refer to the things you want to keep linked
#merge male and female
#levels(tall$Sex)[levels(tall$Sex)=="M"] <- "both"
#levels(tall$Sex)[levels(tall$Sex)=="F"] <- "both"
names(tall)[names(tall) == 'variable'] <- 'Lineage'	
datst <- summarySE(tall, measurevar="value", groupvars=c("Week","Strain","Sex","Irradiation.Dose","Lineage"),na.rm=TRUE)   # summarySE provides the std dev, SEM, and (default 95%) confidence interval. #measurevar is the x-axis
#tp<-datst[datst$Week=='10',]   #Take a single timepoint

#remove NRG 150 Rad
datst <- datst[!(datst$Irradiation.Dose=="150 Rad" & datst$Strain=="NRG"),]

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
PLT<-datst[datst$Lineage == "Platelets",]



#########  PB kinetics line plot ###########

plot<-ggplot(CD45, aes(x=Week,y=value,colour=Strain,shape=Sex,lty=Irradiation.Dose)) + 
  ggtitle("CD45") +
  coord_cartesian(xlim=c(2,31),ylim=c(-0.69897,4)) +
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
  ylab(expression(paste(cells~x~10^3,"/mL"))) +
  scale_x_continuous(breaks=c(3,6,10,20, 30)) +
  annotation_logticks(sides = "l") +
  #facet_wrap( ~Irradiation.Dose, ncol=2) +
#  guides(colour=F,lty=F,shape=F) +
  themePM1()
plot
ggsave(filename="./PB/AMR CB20K38K PB CD45 kinetics.tiff",width=12,height=12,units="cm",dpi=600, compression="lzw")


#plot +
#	coord_cartesian(xlim=c(1,31),ylim=c(-1,2)) 
#ggsave(filename="150606 BM100K BM CD45 axis01 kinetics.tiff",width=8,height=8.1,units="cm",dpi=600, compression="lzw")

plot %+% PLT+ 
  ggtitle("Platelets") +
  coord_cartesian(xlim=c(0,31),ylim=c(1.39794,4)) +
  scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3, 4, 5, 6), labels = c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) 
#  scale_y_continuous(breaks=c(1.39794,2,3,4),labels = c(25, 100, 1000, 10000)) 
ggsave(filename="./PB/AMR CB20K38K PB platelets kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% GM +
  ggtitle("CD15/33") 
ggsave(filename="./PB/AMR CB20K38K PB GM kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")
#ggsave(filename="150505 BM100K PB GM kinetics.wmf",width=8,height=8.1,units="cm",dpi=1200)

plot %+% B +
  ggtitle("CD19") 
ggsave(filename="./PB/AMR CB20K38K PB B kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% T +
  ggtitle("CD3") 
ggsave(filename="./PB/AMR CB20K38K PB T kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

















############################ CB20K - PB kinetic line plots ##################################


#read in data from each experiment:
setwd("C:/Users/paulm/CRC Paul/Experiments/")
Qdat<-read.csv("./2014/140825Q W41 txpt4/figures/QCB20K PB Kinetics.csv", colClasses=c("Sex"="factor"))
Rdat<-read.csv("./2014/141118R W41 again/figures/R_PB kinetics.csv")

#change directory to pooled_data directory
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRG_W41/primary_kinetics_pool/figures")

#add exp label column
Qdat<-cbind("140825Q",Qdat)
Rdat<-cbind("141118R",Rdat)
names(Qdat)[1] <- "Exp"
names(Rdat)[1] <- "Exp"

#merge
dat <- rbind(Qdat, Rdat)

#just take CB20K
dat<-dat[dat$Input=="CB.20K",] 

#write a csv file for record keeping
#write.table(dat, file= "QR CB20K BM kinetics.csv", sep=",", row.names=F)


dat<-dat[,c(1:5,11,13:15,17)]  #Take just the necessary columns


#names(dat)<-c("Exp","Week","Strain","Sex","Irradiation.Dose","CD45+","GPA+","CD33+15+","CD19+","CD3+","CD34+")  #change column names so it looks better on plots
dat[,6:9][dat[,6:9] < 0.2] <- 0.2  #change low values in leukocyte columns to detection threshold
dat[,10][dat[,10] < 25] <- 25  #change low values in platelet column to detection threshold
datlt <- dat  #log10 transform
datlt[6:10] <- log10(dat[6:10])
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
CD45<-datst[datst$Lineage == "CD45",]
GPA<-datst[datst$Lineage == "GPA",]
GM<-datst[datst$Lineage == "CD33.15",]
B<-datst[datst$Lineage == "CD19",]
T<-datst[datst$Lineage == "CD3",]
PLT<-datst[datst$Lineage == "Platelets",]


#########  PB kinetics line plot ###########

plot<-ggplot(CD45, aes(x=Week,y=value,colour=Strain,shape=Sex,lty=Irradiation.Dose)) + 
  ggtitle("CD45") +
  coord_cartesian(xlim=c(2,21),ylim=c(-0.69897,4)) +
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
  ylab(expression(paste(cells~x~10^3,"/mL"))) +
  scale_x_continuous(breaks=c(3,6,10,20, 30)) +
  annotation_logticks(sides = "l") +
  #facet_wrap( ~Irradiation.Dose, ncol=2) +
  #  guides(colour=F,lty=F,shape=F) +
  themePM1()
plot
ggsave(filename="./PB/QR CB20K PB CD45 kinetics.tiff",width=12,height=12,units="cm",dpi=600, compression="lzw")

#plot +
#	coord_cartesian(xlim=c(1,31),ylim=c(-1,2)) 
#ggsave(filename="150606 BM100K BM CD45 axis01 kinetics.tiff",width=8,height=8.1,units="cm",dpi=600, compression="lzw")

plot %+% GM +
  ggtitle("CD15/33") 
ggsave(filename="./PB/QR CB20K PB GM kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")
#ggsave(filename="150505 BM100K BM GM kinetics.wmf",width=8,height=8.1,units="cm",dpi=1200)

plot %+% B +
  ggtitle("CD19") 
ggsave(filename="./PB/QR CB20K PB CD19 kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% T +
	ggtitle("CD3") 
ggsave(filename="./PB/QR CB20K PB CD3 kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")

plot %+% PLT + 
  ggtitle("Platelets") +
  coord_cartesian(xlim=c(0,21),ylim=c(1.39794,4)) +
  scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3, 4, 5, 6), labels = c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) 
#  scale_y_continuous(breaks=c(1.39794,2,3,4),labels = c(25, 100, 1000, 10000)) 
ggsave(filename="./PB/QR CB20K PB PLT kinetics.tiff",width=12,height=12,units='cm',dpi=600,compression="lzw")




