# NRG-W41 Project
# Scripts to Make Figures for Mouse Chimerism Data
# Paul Miller [paulhmiller@gmail.com]


.pardefault <- par(no.readonly=T)  # Stores current par settings
par(.pardefault)  # Reloads stored par settings
dev.off()

library(reshape2)
library(magicaxis)
library(dplyr) # Call this last. Provides: filter, select, do, piping, group_by


## Plotting Functions

# Std Error of the Mean (excludes NAs)
se <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))

  
GroupBarplot <- function(dataframe, week, valueCol, ylab, ytitle, ylim, cols="black", title){
#  Makes a grouped barplot, with strains plotted together and sexes separated.
#  Takes plot title, name of dataframe, time-point to plot, and column number 
#  that contains the values to plot. 
  tmp <- dataframe[dataframe$Week==week, ]
  if(ncol(tmp)==13){
    tmp[8:13] <- log10(tmp[8:13])
  } else if(ncol(tmp)==12){
    tmp[8:12] <- log10(tmp[8:12])
  } else {
    stop("error: function only handles 13 (for BM) or 12 (for PB) columns")
  }
  means <- tapply(tmp[, valueCol], list(tmp$Sex, tmp$Strain), mean, na.rm=TRUE)
  means <- means+5 # this is to enable subzero plotting on log transformed data
  SEMs <- tapply(tmp[, valueCol], list(tmp$Sex, tmp$Strain), se)
  bp <- barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
          xlab="", ylab=ylab, mgp=c(2.2,0.5,0), xpd=FALSE)
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  axis(2, las=2, mgp=c(3,0.6,0), at=c(-2+5, -1+5, 0+5, 1+5, 2+5, 3+5, 4+5, 5+5),
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),
            expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  title(main=title, line=0.5)
  arrows(bp, means+SEMs, bp, means, lwd = 1.5, angle = 90, code = 3, length = 0.05)
  #text(x = bp, y = GM$mean+GM$se+pdist, labels=GM$star , cex=0.7) #paste("p=",round(PLT$p.value,2))
  box()
}
#legend(locator(1),rownames(dat),fill=c("#ee7700","#3333ff"))


KineticsPlot <- function(dataframe, lineage="CD45.Percent", ylab, ylim, xlim=c(2,31), 
                         lty=1, cols="black", pcex=1, lcex=1, title){
  #  Makes a grouped lineplot, with strains plotted together and sexes separated.
  #  Takes plot title, name of dataframe, and column number with values. 
  #  Data is log10 transformed for mean and SEM calculation. 
  #  Requires dplyr and reshape. 
  tmp <- melt(dataframe, id.vars=c("Exp", "Week", "Strain", "Sex", "Irradiation.Dose", "Input", "Mouse"))
  tmp[9] <- log10(tmp[9])
  tmp <- tmp[tmp$variable==lineage, ]
  tmp <- dplyr::summarise(group_by(tmp, Week, Strain, Sex, variable), mean=mean(value, na.rm=TRUE), se=se(value))
  plot(tmp$mean ~ tmp$Week, type="n", axes=F,  ylim=log10(ylim), xlim=xlim, 
       xlab="weeks post-transplant", ylab=ylab, mgp=c(axtitledist,0.5,0))
  sep1 <- tmp[tmp$Strain=="NRG" & tmp$Sex=="M", ]
  sep2 <- tmp[tmp$Strain=="NRG" & tmp$Sex=="F", ]
  sep3 <- tmp[tmp$Strain=="NRG-W41" & tmp$Sex=="M", ]
  sep4 <- tmp[tmp$Strain=="NRG-W41" & tmp$Sex=="F", ]
  points(sep1$mean ~ sep1$Week, cex=pcex, pch=pchs1[1], col=cols[1]) 
  points(sep2$mean ~ sep2$Week, cex=pcex, pch=pchs1[2], col=cols[1]) 
  points(sep3$mean ~ sep3$Week, cex=pcex, pch=pchs1[1], col=cols[2])
  points(sep4$mean ~ sep4$Week, cex=pcex, pch=pchs1[2], col=cols[2])
  lines(sep1$mean ~ sep1$Week, cex=lcex, lty=lty, col=cols[1])
  lines(sep2$mean ~ sep2$Week, cex=lcex, lty=lty, col=cols[1])
  lines(sep3$mean ~ sep3$Week, cex=lcex, lty=lty, col=cols[2])
  lines(sep4$mean ~ sep4$Week, cex=lcex, lty=lty, col=cols[2])
  arrows(tmp$Week, tmp$mean+tmp$se, tmp$Week, tmp$mean-tmp$se, col=c(cols[1],cols[1],cols[2],cols[2]), 
         lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  axis(side=1, at=c(3, 6, 10, 20, 30),  mgp=c(0.8,0.4,0), cex=0.8, tck=-0.05)
  axis(2, las=2, mgp=c(2.5,0.2,0), tck=-0.01, at=c(-2, -1, 0, 1, 2, 3, 4, 5), 
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),expression(10^1),
              expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  title(main=title, line=0.5)
  #text(x = bp, y = GM$mean+GM$se+pdist, labels=GM$star , cex=0.7) #paste("p=",round(PLT$p.value,2))
  box()
}


setwd('C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/figures')


##### CB38K+20K - BM #####

# Read in BM. Uses a csv file that contains BM from multiple experiments
# Each AMQR letter refer to each experiment
BM <- read.csv("../primary_kinetics_pool/AMQR CB20K40K BM kinetics.csv")  #, colClasses=c(Mouse="character"))
# Adjust time-points so that BMa can be pooled:
BM$Week[BM$Week==24] <- 20 # change 24 to 20
BM <- BM[!(BM$Week==18),] # remove week 18. Does not have lineage info. 
BM <- BM[!(BM$Week==26),] # remove week 26. These were 'A' mice, some of which were treated with clodronate at week 25
# Change low values in leukocyte columns to detection threshold
BM[,8:13][BM[,8:13] < 0.01] <- 0.01  
# Remove NRG 150 Rad
BM <- BM[!(BM$Irradiation.Dose=="150 Rad" & BM$Strain=="NRG"),]
# Change/merge 150 and 900 both to irradiated, so that there are only 2 groups (unirradiated and irradiated).
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="0 Rad"] <- "Non-irradiated"
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="150 Rad"] <- "Irradiated" 
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="900 Rad"] <- "Irradiated"
# Check if any na's. If so, need to decide what to do about them
#apply(BM, 2, function(x) sum(is.na(x)))



##### CB38K+20K - PB #####

# Read in PB. Uses a csv file that contains PB from multiple experiments
# Each AMQR letter refer to each experiment
PB <- read.csv("../primary_kinetics_pool/AMQR CB20K40K PB kinetics.csv")  #, colClasses=c(Mouse="character"))
# Adjust time-points so that PBa can be pooled:
PB <- PB[!(PB$Week==26),] # remove week 26. These were 'A' mice, some of which were treated with clodronate at week 25
# Change low values in leukocyte & platelet columns to detection threshold
PB[,8:11][PB[,8:11] < 0.2] <- 0.2 
PB[,12][PB[,12] < 25] <- 25 
# Remove NRG 150 Rad
PB <- PB[!(PB$Irradiation.Dose=="150 Rad" & PB$Strain=="NRG"),]
# Change/merge 150 and 900 both to irradiated, so that there are only 2 groups (unirradiated and irradiated).
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="0 Rad"] <- "Non-irradiated"
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="150 Rad"] <- "Irradiated" 
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="900 Rad"] <- "Irradiated"
# Check if any na's. If so, need to decide what to do about them
#apply(PB, 2, function(x) sum(is.na(x)))



# Select Data
BMi <- BM[BM$Irradiation.Dose=="Irradiated" ,]
PBi <- PB[PB$Irradiation.Dose=="Irradiated" ,]



# Plotting 
ppi <- 300
pchs1 <- c(17,16)
BMylab <- "% total BM cells"
PBylab <- expression(paste(cells~x~10^3,"/mL"))
cols3 <- c("#ee7700","#3333ff")  # colours for bar plots

# Week 20 barplots
#pdf(file="./week20.pdf", width=11.5/2.54, height=8/2.54) #, family='Calibri')
#png("plot.png", width=(11.5*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=10)
par(mfrow=c(2,4), mar=c(2.1, 3.5, 2.1, 1.1), cex=0.7, mgp=c(2,0.6,0))
# BMi
GroupBarplot(BMi, week=20, valueCol=8,  cols=cols3, ylab=BMylab, ylim=c(1, 100), title="CD45")
GroupBarplot(BMi, week=20, valueCol=10, cols=cols3, ylab=BMylab, ylim=c(0.1, 100), title="GM")
GroupBarplot(BMi, week=20, valueCol=11, cols=cols3, ylab=BMylab, ylim=c(0.1, 100), title="B lymphoid")
GroupBarplot(BMi, week=20, valueCol=9,  cols=cols3, ylab=BMylab, ylim=c(0.1, 100), title="GPA")
# PBi
GroupBarplot(PBi, week=20, valueCol=8,  cols=cols3, ylab=PBylab, ylim=c(1, 1000), title="CD45")
GroupBarplot(PBi, week=20, valueCol=9,  cols=cols3, ylab=PBylab, ylim=c(1, 1000), title="GM")
GroupBarplot(PBi, week=20, valueCol=10, cols=cols3, ylab=PBylab, ylim=c(1, 1000), title="B lymphoid")
GroupBarplot(PBi, week=20, valueCol=12, cols=cols3, ylab=PBylab, ylim=c(10, 10000), title="Platelets")
#dev.off()


# Kinetics plots
xlim <- c(1.5, 30.5)   # X-axis range
lcols <- c("#000000","#CD0000")  # line colors
lty <- 1    # linetype (1=solid, 2=dash, 3=dotted)
pcex <- 1.6   # point scaling factor
lcex <- 1   # line scaling factor
axtitledist <- 1.5  # adjusts distance of x and y axis titles

png("plot.png", width=(11.5*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=8)
par(mfrow=c(2,4), mar=c(3, 2.9, 2, 0.3), cex=0.7, mgp=c(2,0.6,0))

# BMi
KineticsPlot(BMi, lineage="CD45.Percent",    ylab=BMylab, ylim=c(1, 100), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot(BMi, lineage="CD33.15.Percent", ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot(BMi, lineage="CD19.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
KineticsPlot(BMi, lineage="GPA.Percent",     ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GPA")
# PBi
KineticsPlot(PBi, lineage="CD45",      ylab=PBylab, ylim=c(1, 1000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot(PBi, lineage="CD33.15",   ylab=PBylab, ylim=c(1, 1000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot(PBi, lineage="CD19",      ylab=PBylab, ylim=c(1, 1000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
KineticsPlot(PBi, lineage="Platelets", ylab=PBylab, ylim=c(25, 10000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="Platelets")
dev.off()














# Statistics
# T.Test
datS2 <- NULL
for (i in c(0,4,24)){
  tmp <- (t.test(dat[dat$time==0, 2], dat[dat$time==i, 2],var.equal = FALSE))
  tmp <- c(i, tmp$p.value)
  datS2 <- rbind(datS2, tmp)
}	

datS2 <- data.frame(datS2)
#datS2[, 2] <- as.numeric(levels(datS2[,2]))[datS2[,2]]

# Add columns for asterisks. (Symbol meanings: . <= 0.10; * <= 0.05; ** <= 0.01; *** <= 0.001)
for (i in 1:nrow(datS2)){
  if (datS2[i, 2] <= 0.001){
    datS2[i ,3] <- "***"
  } else if (datS2[i, 2] <= 0.01){
    datS2[i ,3] <- "**"
  } else if (datS2[i, 2] <= 0.05){
    datS2[i ,3] <- "*"
  } else if (datS2[i, 2] <= 0.10){
    datS2[i ,3] <- "."
  } else{
    datS2[i ,3] <- ""
  }
}
colnames(datS2) <- c("time", "p.value", "star")		








