# Scripts to merge CSV files from 4 experiments into one. 
# Makes separate CSV files for each of BM and PB data. 


# NRG-W41 versus NRG experiments
# Bone Marrow

#read in data from each experiment:
setwd("C:/Users/paulm/CRC Paul/Experiments/2014")
Adat<-read.csv("./140310A B6NRGW41 txpt2/figures/BM kinetics bc.csv")
Mdat<-read.csv("./140708M B6NRGW41 txpt3 & teratoma into W41/figures/MQCB38K BM Kinetics.csv")
Rdat<-read.csv("./141118R W41 again/figures/R_BM kinetics.csv")
Qdat<-read.csv("./140825Q W41 txpt4/figures/QCB20K BM kinetics.csv", colClasses=c(Sex="factor"))
Gdat<-read.csv("../2015/151217G nonirradiated W41/G_BM kinetics.csv", colClasses=c(Sex="factor"))

#change directory to pooled_data directory
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/primary_kinetics_pool/")

#add exp label column
Adat<-cbind("140310A",Adat)
Mdat<-cbind("140708M",Mdat)
Rdat<-cbind("141118R",Rdat)
Qdat<-cbind("140825Q",Qdat)
Gdat<-cbind("151217G",Gdat)

names(Adat)[1] <- "Exp"
names(Mdat)[1] <- "Exp"
names(Rdat)[1] <- "Exp"
names(Qdat)[1] <- "Exp"
names(Gdat)[1] <- "Exp"

#delete non-common columns
Adat<-Adat[,c(1:13)]
Mdat<-Mdat[,c(1:13)]
Rdat<-Rdat[,c(1:13)]
Qdat<-Qdat[,c(1:13)]
Gdat<-Gdat

#merge
dat <- rbind(Adat, Mdat, Rdat, Qdat, Gdat)

#just take CB38K
#dat<-dat[dat$Input!="CB.20K",] 

#adjust time-points so that data can be pooled:
dat$Week[dat$Week==24] <- 20 #change 24 to 20
dat <- dat[!(dat$Week==18),] #remove week 18
dat <- dat[!(dat$Week==26),] #remove week 26

#write a csv file 
write.table(dat, file= "AMQRG CB20K40K BM kinetics.csv", sep=",", row.names=F)


# Peripheral Blood

#read in data from each experiment:
setwd("C:/Users/paulm/CRC Paul/Experiments/2014")
Adat<-read.csv("./140310A B6NRGW41 txpt2/figures/PB kinetics.csv")
Mdat<-read.csv("./140708M B6NRGW41 txpt3 & teratoma into W41/figures/MQCB38K PB Kinetics.csv")
Rdat<-read.csv("./141118R W41 again/figures/R_PB kinetics.csv")
Qdat<-read.csv("./140825Q W41 txpt4/figures/QCB20K PB Kinetics.csv", colClasses=c(Sex="factor"))
Gdat<-read.csv("../2015/151217G nonirradiated W41/G_PB kinetics.csv", colClasses=c(Sex="factor"))


#change directory to pooled_data directory
setwd("C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/primary_kinetics_pool/")

#delete non-common and un-needed columns
Adat<-Adat[,c(1:6,10:14)]
Mdat<-Mdat[,c(1:6,10,12:14,16)]
Rdat<-Rdat[,c(1:6,10,12:14,16)]
Qdat<-Qdat[,c(1:6,10,12:14,16)]
Gdat<-Gdat[,c(2:7,8:12)]

#add exp label column
Adat<-cbind("140310A",Adat)
Mdat<-cbind("140708M",Mdat)
Rdat<-cbind("141118R",Rdat)
Qdat<-cbind("140825Q",Qdat)
Gdat<-cbind("151217G",Gdat)

names(Adat)[1] <- "Exp"
names(Mdat)[1] <- "Exp"
names(Rdat)[1] <- "Exp"
names(Qdat)[1] <- "Exp"
names(Gdat)[1] <- "Exp"


#merge
dat <- rbind(Adat, Mdat, Rdat, Qdat, Gdat)

#just take CB38K
#dat<-dat[dat$Input!="CB.20K",] 

#fix time-points so that data can be pooled:
dat <- dat[!(dat$Week==26),] #remove week 26

#write a csv file for record keeping
write.table(dat, file= "AMQRG CB20K40K PB kinetics.csv", sep=",", row.names=F)



