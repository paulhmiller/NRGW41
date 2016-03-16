source("C:/Users/paulm/Documents/R/source/functions.R")
library(reshape2); library(plyr)

# Get barcode sort info
setwd('C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/barcoding')
bc<-read.csv("barcode samples.csv",header=T,check.names=F ) #,colClasses=c("Experiment"="factor")) #,row.names=1)
names(bc)[names(bc) == 'Timepoint'] <- 'Week' #rename column so it's the same as other input dataframe
bc$Week[bc$Week==24] <- 20 #change 24 to 20

# Get gfp percent info
setwd('C:/Users/paulm/CRC Paul/Experiments/2014/140310A B6NRGW41 txpt2/figures')
A <- read.csv("BM kinetics bc.csv",header=T,check.names=F)
A <- A[,c(1:14)] #get rid of excess columns for merge
# Adjust time-points so data can be pooled:
A$Week[A$Week==24] <- 20 #change 24 to 20
A <- A[!(A$Week==18),] #remove week 18
A <- A[!(A$Week==26),] #remove week 26
setwd('C:/Users/paulm/CRC Paul/Experiments/2014/140708M B6NRGW41 txpt3 & teratoma into W41/figures')
MQ <- read.csv("MQCB38K BM Kinetics.csv",header=T,check.names=F) #,colClasses=c("Experiment"="factor")) #,row.names=1)
setwd('C:/Users/paulm/CRC Paul/Experiments/2014/141118R W41 again/figures')
R <- read.csv("R_BM kinetics.csv",header=T,check.names=F) #,colClasses=c("Experiment"="factor")) #,row.names=1)
R <- R[,c(1:12,14:15)] #get rid of excess columns for merge
gfp <- rbind(A,MQ,R)

# Build new dataframe that merges CD45, CD45.GFP.Percent, GPA & GPA.GFP.Percent with bc counts
dat <- merge(gfp, bc, by= c("Mouse", "Week"))

# Drop out duplicate columns (tried to confirm values were the same but an extra empty 
#  factor level in each prevented this; i checked in Excel instead.):
dat <- dat[, c(1:15,18:ncol(dat))]
#rename columns
names(dat)[names(dat) == 'Strain.x'] <- 'Strain' 
names(dat)[names(dat) == 'Sex.x'] <- 'Sex' 

#make new columns with adjusted barcode counts
bc.E.adj <- round(dat$bc.E * dat$`GPA GFP Percent` / 100, 0)
bc.G.adj <- round(dat$bc.GM * dat$`CD45 GFP Percent` / 100, 0)
bc.B.adj <- round(dat$bc.B * dat$`CD45 GFP Percent` / 100, 0)
bc.T.adj <- round(dat$bc.T * dat$`CD45 GFP Percent` / 100, 0)
dat <- cbind(dat, bc.E.adj, bc.G.adj, bc.B.adj, bc.T.adj)

#write over bc values for 140310A, because their sort was gated on GFP+ cells
dat[dat$Experiment=="140310A", c(24:27)] <- dat[dat$Experiment=="140310A", c(18:21)]

#write csv file
setwd('C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/barcoding')
#write.table(dat, file="out.csv", sep=",", col.names=NA)

#reduce to fewer columns:
tmp <- dat[,c(1:4, 16, 23:27)]
#convert any low barcode sort counts to NA:
tmp[,7:ncol(tmp)][tmp[,7:ncol(tmp)] < 21] <- NA
#get rid of rows where no counts above threshold (NA)
tmp <- tmp[rowSums(is.na(tmp)) != 4,]
#exclude NRG 150 Rad
tmp <- tmp[tmp$Strain!="NRG" | tmp$Irradiation!="150 Rad",]

#write csv file
#write.table(tmp, file="barcodes_with_threshold.csv", sep=",", col.names=NA)


## Remove mice and samples to get numbers down:
#tmp <- tmp[tmp$Week!=6, ]
#tmp <- tmp[tmp$Week!=30, ]
# Convert to single row per collection tube, and remove NA's
tmp<- melt(tmp,id.vars=c("Mouse","Week","Strain","Sex","Irradiation","BOX"),na.rm=TRUE)


# Write csv file
setwd('C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/barcoding')
#write.table(tmp, file="barcodes_with_threshold_tubes2.csv", sep=",", col.names=NA)


### At this point, tubes were ordered and labelled with new unique ID numbers. 
### These new numbers are in the csv file as "tubeID2'
### New boxes were made to store potentially usable tubes. 
### MS Excel was used to unput the above information as the tubes were labelled.
### Below are scripts to further reduce tube numbers. 

setwd('C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/barcoding')
dat <- read.csv("barcodes_send_160314.csv", header=TRUE)


dat$send <- "send"
# General discard conditions:

dat$send[dat$Week==6] <- "discard"
dat$send[is.na(dat$tubeID2)] <- "discard"
dat$send[dat$value < 50] <- "discard"

# Specific discard decisions:

# There's only one week 30 non-irrad NRG tube (also, it's only 57 cells):
dat$send[dat$Week==30 & dat$Strain=="NRG" & dat$Irradiation=="0 Rad"] <- "discard"
# Excessive NRG F 900 Rad mice. Choose 2 at random and remove them. 
set.seed(1)
discard <- sample(unique(dat$Mouse[dat$Week==3 & dat$Strain=="NRG" & dat$Sex=="F" & dat$Irradiation=="900 Rad"]), 2)
dat$send[dat$Mouse==discard[1] & dat$Week==3 & dat$Strain=="NRG" & dat$Sex=="F" & dat$Irradiation=="900 Rad"] <- "discard"
dat$send[dat$Mouse==discard[2] & dat$Week==3 & dat$Strain=="NRG" & dat$Sex=="F" & dat$Irradiation=="900 Rad"] <- "discard"

write.table(dat, file="barcodes_send_160314_v2.csv", sep=",", col.names=NA)


nrow(dat[dat$send=="send", ])
tmp <- dat[order(dat$Week, dat$Irradiation, dat$Strain, dat$Sex,  dat$Mouse),]
tmp[tmp$send=="send", ]












### Sequencing logistics ###

# 100 per lane




 