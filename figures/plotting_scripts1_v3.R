# NRG-W41 Project
# Scripts to Make Figures for Mouse Chimerism Data
# Paul Miller [paulhmiller@gmail.com]

## info:
# Elsevier: 1 column, 90 mm; 1.5 column, 140 mm; and 2 column, 190 mm (the full width of the page)
# Blood: 1 column, 80 mm; 1.5 column, 115 mm

.pardefault <- par(no.readonly=T)  # Stores current par settings
par(.pardefault)  # Reloads stored par settings
dev.off()

library(reshape2)
library(magicaxis)
library(reshape2)
library(dplyr) # Call this last. Provides: filter, select, do, pipe, group_by


setwd('C:/Users/paulm/CRC Paul/PROJECTS/NRGW41/figures')

## Functions
# Std Error of the Mean (excludes NAs)
se <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
# Makes a table of significance stars from an input table of p-values
p2stars <- function(table){
    stars <- table
    stars[stars<=0.001]              <- "***"
    stars[stars<=0.01 & stars>0.001] <- "**"
    stars[stars<=0.05 & stars>0.01]  <- "*"
    stars[stars<=0.10 & stars>0.05]  <- "."
    stars[stars> 0.10]               <- "-"
    return(stars)
}


# Global Plotting Parameters
# Standard widths for figures: 1 column, 85 mm; 1.5 column, 114 mm; and 2 column, 174 mm (the full width of the page).
ppi <- 300
pchs1 <- c(17,16)
BMylab <- "% total BM cells"
PBylab <- expression(paste(cells~x~10^3,"/mL"))

xlim <- c(1.5, 30.5)   # X-axis range on line plots
lcols <- c("#000000","#CD0000")  # line colors
lty <- 1    # linetype (1=solid, 2=dash, 3=dotted)
pcex <- 1.8   # point scaling factor
lcex <- 1   # line scaling factor
axtitledist <- 1.6  # adjusts distance of x and y axis titles on kinetics plots
vAdj <- 0.0  # Jitter-like effect on x-axis, use values between 0-2





# Function for Old Age plots (not in use):
KineticsPlot2 <- function(dataframe, lineage="Human.Percent", ylab, ylim, xlim=c(2,31), 
                          lty=1, cols="black", pcex=1, lcex=1, title){
  #  Makes a grouped lineplot, with strains plotted together and sexes separated.
  #  Takes plot title, name of dataframe, and column number with values. 
  #  Data is log10 transformed for mean and SEM calculation. 
  #  Requires dplyr and reshape. 
  tmp <- melt(dataframe, id.vars=variables)
  tmp[9] <- log10(tmp[9])
  tmp <- tmp[tmp$variable==lineage, ]
  tmp <- dplyr::summarise(group_by(tmp, Week, Sex, variable), mean=mean(value, na.rm=TRUE), se=se(value))
  plot(tmp$mean ~ tmp$Week, type="n", axes=F,  ylim=log10(ylim), xlim=xlim, col=tmp$Strain,
       xlab="weeks post-transplant", ylab=ylab, mgp=c(axtitledist,0.5,0))
  sep1 <- tmp[tmp$Sex=="F", ]
  sep1$Week <- sep1$Week + (vAdj/1)
  sep2 <- tmp[tmp$Sex=="M", ]
  sep2$Week <- sep2$Week + (vAdj/2)
  print (sep1)
  arrows(sep1$Week, sep1$mean+sep1$se, sep1$Week, sep1$mean-sep1$se, col=c(cols[1]), 
         lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  arrows(sep2$Week, sep2$mean+sep2$se, sep2$Week, sep2$mean-sep2$se, col=c(cols[1]), 
         lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  points(sep1$mean ~ sep1$Week, cex=pcex, pch=pchs1[1], col=cols[1]) 
  points(sep2$mean ~ sep2$Week, cex=pcex, pch=pchs1[2], col=cols[1]) 
  lines(sep1$mean ~ sep1$Week, cex=lcex, lty=lty, col=cols[1])
  lines(sep2$mean ~ sep2$Week, cex=lcex, lty=lty, col=cols[1])
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  axis(side=1, at=c(3, 6, 10, 20, 30),  mgp=c(0.8,0.4,0), cex=0.8, tck=-0.03)
  axis(2, las=2, mgp=c(2.5,1.7,0), tck=-0.01, hadj=0, at=c(-2, -1, 0, 1, 2, 3, 4, 5), 
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),expression(10^1),
                expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  title(main=title, line=0.5)
  #text(x = bp, y = GM$mean+GM$se+pdist, labels=GM$star , cex=0.7) #paste("p=",round(PLT$p.value,2))
  box()
}
#variables <- names(BM[0:7])
#KineticsPlot2(BM, lineage="Human.Percent",    ylab=BMylab, ylim=c(1, 100), xlim=xlim, 
#              cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")











##### NRG versus NRG-W41 #####

##### BM #####
# Read in BM. Uses a csv file that contains BM from multiple experiments
# Each AMQR letter refer to each experiment
BM <- read.csv("../primary_kinetics_pool/AMQRG CB20K40K BM kinetics.csv")  #, colClasses=c(Mouse="character"))
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


# Kinetics plots
# Function for NRG verus NRG-W41 plots:
KineticsPlot1 <- function(dataframe, lineage="CD45.Percent", ylab, ylim, xlim=c(2,31), 
                         lty=1, cols="black", pcex=1, lcex=1, title){
  #  Makes a grouped lineplot, with strains plotted together and sexes separated.
  #  Takes plot title, name of dataframe, and column number with values. 
  #  Data is log10 transformed for mean and SEM calculation. 
  #  Requires dplyr and reshape. 
  tmp <- melt(dataframe, id.vars=c("Exp", "Week", "Strain", "Sex", "Irradiation.Dose", "Input", "Mouse"))
  tmp[9] <- log10(tmp[9])
  tmp <- tmp[tmp$variable==lineage, ]
  tmp <- dplyr::summarise(group_by(tmp, Week, Strain, Sex, variable), mean=mean(value, na.rm=TRUE), se=se(value))
  plot(tmp$mean ~ tmp$Week, type="n", axes=F,  ylim=log10(ylim), xlim=xlim, col=tmp$Strain,
       xlab="weeks post-transplant", ylab=ylab, mgp=c(axtitledist,0.5,0))
  text(x = weeks, y=sig, labels=stars[ ,lineage] , cex=0.9) #paste("p=",round(PLT$p.value,2))
  sep1 <- tmp[tmp$Strain=="NRG" & tmp$Sex=="M", ]
  sep1$Week <- sep1$Week + (vAdj/1)
  sep2 <- tmp[tmp$Strain=="NRG" & tmp$Sex=="F", ]
  sep2$Week <- sep2$Week + (vAdj/2)
  sep3 <- tmp[tmp$Strain=="NRG-W41" & tmp$Sex=="M", ]
  sep3$Week <- sep3$Week + (-vAdj/2)
  sep4 <- tmp[tmp$Strain=="NRG-W41" & tmp$Sex=="F", ]
  sep4$Week <- sep4$Week + (-vAdj/1) 
  arrows(sep1$Week, sep1$mean+sep1$se, sep1$Week, sep1$mean-sep1$se, col=c(cols[1]), 
        lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  arrows(sep2$Week, sep2$mean+sep2$se, sep2$Week, sep2$mean-sep2$se, col=c(cols[1]), 
        lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  arrows(sep3$Week, sep3$mean+sep3$se, sep3$Week, sep3$mean-sep3$se, col=c(cols[2]), 
        lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  arrows(sep4$Week, sep4$mean+sep4$se, sep4$Week, sep4$mean-sep4$se, col=c(cols[2]), 
        lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  points(sep1$mean ~ sep1$Week, cex=pcex, pch=pchs1[1], col=cols[1]) 
  points(sep2$mean ~ sep2$Week, cex=pcex, pch=pchs1[2], col=cols[1]) 
  points(sep3$mean ~ sep3$Week, cex=pcex, pch=pchs1[1], col=cols[2])
  points(sep4$mean ~ sep4$Week, cex=pcex, pch=pchs1[2], col=cols[2])
  lines(sep1$mean ~ sep1$Week, cex=lcex, lty=lty, col=cols[1])
  lines(sep2$mean ~ sep2$Week, cex=lcex, lty=lty, col=cols[1])
  lines(sep3$mean ~ sep3$Week, cex=lcex, lty=lty, col=cols[2])
  lines(sep4$mean ~ sep4$Week, cex=lcex, lty=lty, col=cols[2])
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  axis(side=1, at=c(3, 6, 10, 20, 30),  mgp=c(0.8,0.4,0), cex=0.8, tck=-0.03)
  axis(2, las=2, mgp=c(2.5,1.7,0), tck=-0.01, hadj=0, at=c(-2, -1, 0, 1, 2, 3, 4, 5), 
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),expression(10^1),
              expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  title(main=title, line=0.5)
  #text(x = bp, y = GM$mean+GM$se+pdist, labels=GM$star , cex=0.7) #paste("p=",round(PLT$p.value,2))
  box()
}

lcols <- c("#000000","#CD0000")  # line colors
lty <- 1    # linetype (1=solid, 2=dash, 3=dotted)
vAdj <- 0.0  # Jitter-like effect, use values between 0-2




# BM
png("figX_W41_BM_kinetics_1.5col.png", width=(14.0*ppi)/2.54, height=(16*ppi)/2.54, res=ppi, pointsize=8)
par(mfcol=c(4,4), mar=c(3.2, 2.9, 2, 0.8), cex=0.7, mgp=c(2,0.6,0))
sig <- 2

# f, i
tmp <- BM[BM$Irradiation.Dose=="Irradiated" & BM$Sex=="F",]
# Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:4, 6)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out)
stars <- p2stars(out)
# Plots
KineticsPlot1(tmp, lineage="CD45.Percent",    ylab=BMylab, ylim=c(0.1, 100), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15.Percent", ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
KineticsPlot1(tmp, lineage="GPA.Percent",     ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GPA")

# m, i
tmp <- BM[BM$Irradiation.Dose=="Irradiated" & BM$Sex=="M",]
# Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:4, 6)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
KineticsPlot1(tmp, lineage="CD45.Percent",    ylab=BMylab, ylim=c(0.1, 100), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15.Percent", ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
KineticsPlot1(tmp, lineage="GPA.Percent",     ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GPA")

# f, ni
tmp <- BM[BM$Irradiation.Dose=="Non-irradiated"  & BM$Sex=="F",]
# Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20)    #,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:4, 6)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
KineticsPlot1(tmp, lineage="CD45.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15.Percent", ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
KineticsPlot1(tmp, lineage="GPA.Percent",     ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GPA")

# m, ni
tmp <- BM[BM$Irradiation.Dose=="Non-irradiated"  & BM$Sex=="M",]
# Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20)   #,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:4, 6)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
KineticsPlot1(tmp, lineage="CD45.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15.Percent", ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
KineticsPlot1(tmp, lineage="GPA.Percent",     ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GPA")
dev.off()





# PB
png("figX_W41_PB_kinetics_1.5col.png", width=(14.0*ppi)/2.54, height=(16*ppi)/2.54, res=ppi, pointsize=8)
par(mfcol=c(4,4), mar=c(3.2, 2.9, 2, 0.8), cex=0.7, mgp=c(2,0.6,0))
ylims1 <- c(0.3,1400)
ylims2 <- c(30,12000)
sig1 <- log10(1400)
sig2 <- log10(12000)

# f, i
tmp <- PB[PB$Irradiation.Dose=="Irradiated" & PB$Sex=="F",]
## Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:3, 5)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
sig <- sig1
KineticsPlot1(tmp, lineage="CD45",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15",   ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
sig <- sig2
KineticsPlot1(tmp, lineage="Platelets", ylab=PBylab, ylim=ylims2, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="Platelets")

# m, i
tmp <- PB[PB$Irradiation.Dose=="Irradiated" & PB$Sex=="M",]
## Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20)  #,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:3, 5)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
sig <- sig1
KineticsPlot1(tmp, lineage="CD45",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15",   ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
sig <- sig2
KineticsPlot1(tmp, lineage="Platelets", ylab=PBylab, ylim=ylims2, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="Platelets")

# f, ni
tmp <- PB[PB$Irradiation.Dose=="Non-irradiated"  & PB$Sex=="F",]
## Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(10)   #3,6,10,20,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:3, 5)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
sig <- sig1
KineticsPlot1(tmp, lineage="CD45",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15",   ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
sig <- sig2
KineticsPlot1(tmp, lineage="Platelets", ylab=PBylab, ylim=ylims2, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="Platelets")

# m ni
tmp <- PB[PB$Irradiation.Dose=="Non-irradiated"  & PB$Sex=="M",]
## Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(10)   #3,6,10,20,30)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:3, 5)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
sig <- sig1
KineticsPlot1(tmp, lineage="CD45",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45")
KineticsPlot1(tmp, lineage="CD33.15",   ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM")
KineticsPlot1(tmp, lineage="CD19",      ylab=PBylab, ylim=ylims1, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid")
sig <- sig2
KineticsPlot1(tmp, lineage="Platelets", ylab=PBylab, ylim=ylims2, xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="Platelets")
dev.off()





# T cells



png("figX_W41_Tcells_1.5col.png", width=(14.0*ppi)/2.54, height=(16*ppi)/2.54, res=ppi, pointsize=8)
par(mfcol=c(4,4), mar=c(3.2, 2.9, 2, 0.8), cex=0.7, mgp=c(2,0.6,0))
sig <- 2
## BM# f, i
tmp <- BM[BM$Irradiation.Dose=="Irradiated" & BM$Sex=="F",]
# Statistics
dat <- tmp
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20,30)
lineages <- names(dat[8:ncol(dat)])
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp2 <- (t.test(dat[dat$Strain=="NRG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG-W41" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp2$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
stars <- p2stars(out)
# Plots
KineticsPlot1(tmp, lineage="CD3.Percent",    ylab=BMylab, ylim=c(0.1, 100), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD3")


tmp <- BM[BM$Week==30, ]
tmp <- tmp[, c(3:5,12)]
tmp[4] <- log10(tmp[4])




# Couldn't figure out how to separate the dataframe base on multiple factors,
# so am doing it manually:
tmp1 <- dplyr::filter(tmp, Strain=="NRG-W41" & Sex=="F" & Irradiation.Dose=="Irradiated")
tmp2 <- dplyr::filter(tmp, Strain=="NRG"     & Sex=="F" & Irradiation.Dose=="Irradiated")
tmp3 <- dplyr::filter(tmp, Strain=="NRG-W41" & Sex=="M" & Irradiation.Dose=="Irradiated")
tmp4 <- dplyr::filter(tmp, Strain=="NRG"     & Sex=="M" & Irradiation.Dose=="Irradiated")
tmp5 <- dplyr::filter(tmp, Strain=="NRG-W41" & Sex=="F" & Irradiation.Dose=="Non-irradiated")
tmp6 <- dplyr::filter(tmp, Strain=="NRG"     & Sex=="F" & Irradiation.Dose=="Non-irradiated")
tmp7 <- dplyr::filter(tmp, Strain=="NRG-W41" & Sex=="M" & Irradiation.Dose=="Non-irradiated")
tmp8 <- dplyr::filter(tmp, Strain=="NRG"     & Sex=="M" & Irradiation.Dose=="Non-irradiated")
plot(seq(1:8), type="n", axes=F) 
points(cfcNSG$percent ~ jitter(as.numeric(cfcNSG$Input),factor=0.5), cex=pcex, pch=pch1, col=col1)

#myFun <- function(x){return(x)}
#tmp2 <- tapply(tmp[,4], list(tmp$Strain, tmp$Sex, tmp$Irradiation.Dose), myFun)
#means <- tapply(tmp[, valueCol], list(tmp$Age, tmp$Sex), mean, na.rm=TRUE)
#  means <- means+5 # this is to enable subzero plotting on log transformed data
#SEMs <- tapply(tmp[, valueCol], list(tmp$Age, tmp$Sex), se)
plot <- barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
                xlab="", ylab=ylab, mgp=c(axtitledist,0.5,0), xpd=FALSE)
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  par(new=TRUE)  # enables replotting over the yaxis ticks, using the next line 
  barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
          xlab="recipient sex", ylab="", mgp=c(axtitledist,0.5,0), xpd=FALSE)
  par(new=FALSE)
  axis(2, las=2, mgp=c(3,1.7,0), tck=-0.02, hadj=0, at=c(-2+5, -1+5, 0+5, 1+5, 2+5, 3+5, 4+5, 5+5),
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),
                expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  axis(1, las=2, tck=-0.02, at=c(2,5), labels=c("",""))
  title(main=title, line=0.5)
  arrows(bp, means+SEMs, bp, means-SEMs, lwd = 1.0, angle = 90, code = 3, length = 0.025)
  text(x=2, y=sig1+5, labels=starsF[, valueCol], cex=0.9) # add significance stars
  text(x=5, y=sig1+5, labels=starsM[, valueCol], cex=0.9) # add significance stars
  box()
}
#legend(locator(1),rownames(dat),fill=c("#ee7700","#3333ff"))





#tmp <- dplyr::summarise(group_by(cd3, Strain, Sex, Irradiation.Dose), 
#         mean=mean(CD3.Percent, na.rm=TRUE), se=se(CD3.Percent))
# Make plot
#png("figX_cd3.png", width=(9.0*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=8)
plot(tmp, yaxt="n", xaxt="n", ylim=c(0,1),
     xlab="weeks post-irradiation", ylab="fraction alive", 
     mgp=c(axtitledist+0.4,0.5,0), xpd=FALSE)
axis(2, las=2, mgp=c(0,1.7,0), tck=-0.02, hadj=0, at=seq(0,1,0.2) )
axis(side=1, at=seq(0,6,1),  mgp=c(0.4,0.4,0), cex=0.8, tck=-0.03)
xplace <- 5.5
toplab <- 1.02
step <- 0.06
units <- "cGy"
text(x=xplace, y=toplab-1*step, cex=0.9, labels=paste(doses[1],units))
text(x=xplace, y=toplab-2*step, cex=0.9, labels=paste(doses[2],units))
text(x=xplace, y=toplab-3*step, cex=0.9, labels=paste(doses[3],units))
text(x=xplace, y=toplab-4*step, cex=0.9, labels=paste(doses[5],units))
text(x=3.5, y=0.1, cex=0.9, labels=paste(doses[4],units))
#dev.off()



#cfcS <- merge(cfcS, data.frame(cfcstats), by = "Input")
#cfcS$Input<-factor(cfcS$Input, levels=c("5GF", "3GF", "2%FBS", "fresh"))
#cfcS <- cfcS[order(cfcS$Input), ]

# Separate out arms
cfcScntl <- cfcS[cfcS$Input=="fresh" ,]

# Create Chart
plot(seq(1:3), type="n", axes=F, xlim=xlims, ylim=ylim1, ylab=ylab1, xlab=xlab1, mgp=mpg1) 
#cfcS$mean ~ cfcS$Input, 
title(main="", line = titledist)
rect(-1, 100-cfcScntl$se, 10, 100+cfcScntl$se, col="gray", border="gray")
points(cfcNSG$percent ~ jitter(as.numeric(cfcNSG$Input),factor=0.5), cex=pcex, pch=pch1, col=col1) 
points(cfc3GS$percent ~ jitter(as.numeric(cfc3GS$Input),factor=0.5), cex=pcex, pch=pch2, col=col2) 
#axis(side=1, at=c(1:3), labels=rep("",3), mgp=c(3,0.5,0)) # makes tick marks
text(cex=1, x=c(1:3)+0.3, y=-15, labels=xlab2, xpd=TRUE, srt=45, pos=2)
axis(side=2, las=2, mgp=c(3,0.6,0))
axis(side=1, at=seq(1:3), labels=rep(NA,3), mgp=c(2,1,0))  # Axis ticks
# Error bars:
segments(x0=c(1:3)-0.4, y0=cfcS$mean, x1=c(1:3)+0.4, y1=cfcS$mean, lwd = 1.5, col=col3)
#segments(x0=c(1:3), y0=cfcS$mean-cfcS$se, x1=c(1:3), y1=cfcS$mean+cfcS$se, lwd = 1.5, col=col3)
arrows(x0=c(1:3), y0=cfcS$mean- cfcS$se, x1=c(1:3), y1=cfcS$mean + cfcS$se, lwd = 1.5, angle = 90, code = 3, length = 0.05, col=col3)
abline(h=100, lty=3, lwd = 1.2) 
text(x = c(1:3), y = cfcS$mean+cfcS$se+pdist, labels =cfcS$star, cex=0.7) #paste("p=",round(cfcstats$p.value,2))
box()




# Unused W41 plot funciton:
# Week 20 barplot
# Function for NRG versus NRG-W41 Barplot
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
  print(means)
  means <- means+5 # this is to enable subzero plotting on log transformed data
  SEMs <- tapply(tmp[, valueCol], list(tmp$Sex, tmp$Strain), se)
  bp <- barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
          xlab="", ylab=ylab, mgp=c(2.2,0.5,0), xpd=FALSE)
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  par(new=TRUE)  # enables replotting over the yaxis ticks, using the next line
  barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
          xlab="", ylab="", mgp=c(2.2,0.5,0), xpd=FALSE) 
  par(new=FALSE)
  axis(2, las=2, mgp=c(3,0.6,0), hadj=0, at=c(-2+5, -1+5, 0+5, 1+5, 2+5, 3+5, 4+5, 5+5),
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),
            expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  title(main=title, line=0.5)
  arrows(bp, means+SEMs, bp, means, lwd = 1.5, angle = 90, code = 3, length = 0.05)
  #text(x = bp, y = GM$mean+GM$se+pdist, labels=GM$star , cex=0.7) #paste("p=",round(PLT$p.value,2))
  box()
}
#legend(locator(1),rownames(dat),fill=c("#ee7700","#3333ff"))


cols3 <- c("#ee7700","#3333ff")  # colours for bar plots

#pdf(file="./week20.pdf", width=11.5/2.54, height=8/2.54) #, family='Calibri')
#png("plot1.png", width=(11.5*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=10)
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
# end of unused section





##### NSG versus NRG #####
##### BM #####

BM <- read.csv("../NSGvsNRG/C_BM kinetics.csv", colClasses=c(Sex="factor"))
# Adjust time-points so that BMa can be pooled:
BM$Week[BM$Week==24] <- 20
BM$Week[BM$Week==8] <- 10 
BM <- BM[!(BM$Week==17),] #remove week 17
# Change low values in leukocyte columns to detection threshold
BM[,8:13][BM[,8:13] < 0.01] <- 0.01  
#Change irradiation doses to enable merge
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="250 Rad"] <- "Irradiated" 
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="800 Rad"] <- "Irradiated"
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="315 Rad"] <- "Irradiated" 
levels(BM$Irradiation.Dose)[levels(BM$Irradiation.Dose)=="900 Rad"] <- "Irradiated"

##### PB #####

PB <- read.csv("../NSGvsNRG/C_PB kinetics.csv", colClasses=c(Sex="factor"))
# Adjust time-points so that BMa can be pooled:
PB$Week[PB$Week==24] <- 20
PB$Week[PB$Week==8] <- 10 
PB <- PB[!(PB$Week==17),] #remove week 17
# Change low values in leukocyte & platelet columns to detection threshold
PB[,8:11][PB[,8:11] < 0.2] <- 0.2 
PB[,12][PB[,12] < 25] <- 25 
#Change irradiation doses to enable merge
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="250 Rad"] <- "Irradiated" 
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="800 Rad"] <- "Irradiated"
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="315 Rad"] <- "Irradiated" 
levels(PB$Irradiation.Dose)[levels(PB$Irradiation.Dose)=="900 Rad"] <- "Irradiated"

#Statistics
dat <- BM
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:4, 6)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp <- (t.test(dat[dat$Strain=="NSG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out)
BMstars <- p2stars(out)

dat <- PB
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
weeks <- c(3,6,10,20)
lineages <- names(dat[8:ncol(dat)])
lineages <- lineages[c(1:3, 5)]
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp <- (t.test(dat[dat$Strain=="NSG" & dat$Week==w, l], 
                   dat[dat$Strain=="NRG" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out) 
PBstars <- p2stars(out)

# Do a separate t.test for PB CD19 w3 to check things
tmp <- PB[PB$Week==3 ,c(3,10)]
tmp[,2] <- log10(tmp[,2])
NSG <- tmp[tmp$Strain=="NSG", 2]
NRG <- tmp[tmp$Strain=="NRG", 2]
t.test(NSG, NRG, var.equal = FALSE)

# Kinetics plots
# Function for NSG verus NRG plots:
KineticsPlot2 <- function(dataframe, lineage="CD45.Percent", ylab, ylim, xlim=c(2,31), 
                          lty=1, cols="black", pcex=1, lcex=1, title, sig){
  #  Makes a grouped lineplot, with strains plotted together and sexes separated.
  #  Takes plot title, name of dataframe, and column number with values. 
  #  Data is log10 transformed for mean and SEM calculation. 
  #  Requires dplyr and reshape. 
  tmp <- melt(dataframe, id.vars=c("Exp", "Week", "Strain", "Sex", "Irradiation.Dose", "Input", "Mouse"))
  tmp[9] <- log10(tmp[9])
  tmp <- tmp[tmp$variable==lineage, ]
  tmp <- dplyr::summarise(group_by(tmp, Week, Strain, variable), mean=mean(value, na.rm=TRUE), se=se(value))
  plot(tmp$mean ~ tmp$Week, type="n", axes=F,  ylim=log10(ylim), xlim=xlim, col=tmp$Strain,
       xlab="weeks post-transplant", ylab=ylab, mgp=c(axtitledist,0.5,0))
  text(x = unique(tmp$Week), y=sig, labels=stars[ ,lineage] , cex=0.9) #paste("p=",round(PLT$p.value,2))
  sep1 <- tmp[tmp$Strain=="NRG", ]
  sep1$Week <- sep1$Week + (vAdj/1)
  sep2 <- tmp[tmp$Strain=="NSG", ]
  sep2$Week <- sep2$Week + (vAdj/2)
  arrows(sep1$Week, sep1$mean+sep1$se, sep1$Week, sep1$mean-sep1$se, col=c(cols[1]), 
         lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  arrows(sep2$Week, sep2$mean+sep2$se, sep2$Week, sep2$mean-sep2$se, col=c(cols[1]), 
         lwd = 1.5, angle = 90, code = 3, length = 0.02) 
  points(sep1$mean ~ sep1$Week, cex=pcex, pch=pchs1[1], col=cols[1]) 
  points(sep2$mean ~ sep2$Week, cex=pcex, pch=pchs1[2], col=cols[1]) 
  lines(sep1$mean ~ sep1$Week, cex=lcex, lty=lty, col=cols[1])
  lines(sep2$mean ~ sep2$Week, cex=lcex, lty=lty, col=cols[1])
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  axis(side=1, at=c(3, 6, 10, 20, 30),  mgp=c(0.8,0.4,0), cex=0.8, tck=-0.03)
  axis(2, las=2, mgp=c(2.5,1.7,0), tck=-0.01, hadj=0, at=c(-2, -1, 0, 1, 2, 3, 4, 5), 
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),expression(10^1),
                expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  title(main=title, line=0.5)
  box()
}

xlim <- c(1.5, 20.5)   # X-axis range
lcols <- c("#000000","#CD0000")  # line colors
lty <- 1    # linetype (1=solid, 2=dash, 3=dotted)
pchs1 <- c(16,1)

png("figX_NSG_NRG_1.5col.png", width=(14.0*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=8)
par(mfrow=c(2,4), mar=c(3.2, 2.9, 2, 0.8), cex=0.7, mgp=c(2,0.6,0))

# BM
stars <- BMstars
KineticsPlot2(BM, lineage="CD45.Percent",    ylab=BMylab, ylim=c(1, 100), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45", sig=2)
KineticsPlot2(BM, lineage="CD33.15.Percent", ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM", sig=2)
KineticsPlot2(BM, lineage="CD19.Percent",    ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid", sig=2)
KineticsPlot2(BM, lineage="GPA.Percent",     ylab=BMylab, ylim=c(0.1, 100),  xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GPA", sig=2)
# PB
stars <- PBstars
KineticsPlot2(PB, lineage="CD45",      ylab=PBylab, ylim=c(1, 1000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="CD45", sig=3)
KineticsPlot2(PB, lineage="CD33.15",   ylab=PBylab, ylim=c(1, 1000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="GM", sig=3)
KineticsPlot2(PB, lineage="CD19",      ylab=PBylab, ylim=c(1, 1000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="B lymphoid", sig=3)
KineticsPlot2(PB, lineage="Platelets", ylab=PBylab, ylim=c(30, 10000), xlim=xlim, 
             cols=lcols, pcex=pcex, lcex=lcex, lty=lty, title="Platelets", sig=4)
dev.off()














##### Old Age #####
##### BM #####
tmp <- read.csv("../old_age/OldAge1_BM.csv", strip.white=TRUE)  #, colClasses=c(Mouse="character"))
Exp <- rep("OldAge1", nrow(tmp))
tmp <- cbind(Exp, tmp)
BM1 <- tmp[c(1,6, 7:11, 22:27)]

tmp <- read.csv("../old_age/OldAge2_BM.csv", strip.white=TRUE)  #, colClasses=c(Mouse="character"))
Exp <- rep("OldAge2", nrow(tmp))
tmp <- cbind(Exp, tmp)
BM2 <- tmp[c(1,6, 7:11, 22:27)]

# Adjust time-points so that BMa can be pooled:
BM2$Week <- round(BM2$Week, digits=0)
BM2$Week[BM2$Week==21] <- 20
BM2$Week[BM2$Week==25] <- 20


# Merge experiment 1 & 2:
BM <- rbind(BM1, BM2)

# Make codes the same:
BM$Age[BM$Age=="O"] <- "Old"
BM$Age[BM$Age=="Y"] <- "Young"
levels(BM$Age)[levels(BM$Age)=="Old"] <- "old"
levels(BM$Age)[levels(BM$Age)=="Young"] <- "young"
BM$Carriers[BM$Carriers=="N"] <- "no"
BM$Carriers[BM$Carriers=="Y"] <- "yes"
BM <- droplevels(BM)
# reorder factor levels:
BM$Age <- factor(BM$Age,levels(BM$Age)[c(2,1)])
#print(levels(BM$Age))

# Change low values in leukocyte columns to detection threshold
BM[,8:13][BM[,8:13] < 0.01] <- 0.01  

# Select Data
BM34 <- BM[BM$Cell.type=="CD34+" ,]
BM34 <- BM34[BM34$Week==20, ]
#BM49f50 <- BM[BM$Cell.type=="50 49f+" ,]
#BM49f60 <- BM[BM$Cell.type=="60 49f+" ,]
#BM49f <- rbind(BM49f50, BM49f60)

# Function for barplot of Old versus Young, subdivided by Sex.
GroupBarplot2 <- function(dataframe, week, valueCol, ylab, ytitle, ylim, cols="black", title){
  #  Makes a grouped barplot, time-points plotted together and sexes separated.
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
  means <- tapply(tmp[, valueCol], list(tmp$Age, tmp$Sex), mean, na.rm=TRUE)
  print(means)
  means <- means+5 # this is to enable subzero plotting on log transformed data
  SEMs <- tapply(tmp[, valueCol], list(tmp$Age, tmp$Sex), se)
  bp <- barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
                xlab="", ylab=ylab, mgp=c(axtitledist,0.5,0), xpd=FALSE)
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  par(new=TRUE)  # enables replotting over the yaxis ticks, using the next line 
  barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
          xlab="recipient sex", ylab="", mgp=c(axtitledist,0.5,0), xpd=FALSE)
  par(new=FALSE)
  axis(2, las=2, mgp=c(3,1.7,0), tck=-0.02, hadj=0, at=c(-2+5, -1+5, 0+5, 1+5, 2+5, 3+5, 4+5, 5+5),
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),
                expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  axis(1, las=2, tck=-0.02, at=c(2,5), labels=c("",""))
  title(main=title, line=0.5)
  arrows(bp, means+SEMs, bp, means-SEMs, lwd = 1.0, angle = 90, code = 3, length = 0.025)
  text(x=2, y=sig1+5, labels=starsF[, valueCol], cex=0.9) # add significance stars
  text(x=5, y=sig1+5, labels=starsM[, valueCol], cex=0.9) # add significance stars
  box()
}
#legend(locator(1),rownames(dat),fill=c("#ee7700","#3333ff"))



# Function for barplot of +/- carriers, subdivided by sex. 
GroupBarplot4 <- function(dataframe, week, valueCol, ylab, ytitle, ylim, cols="black", title){
  #  Makes a grouped barplot, time-points plotted together and sexes separated.
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
  means <- tapply(tmp[, valueCol], list(tmp$Carriers, tmp$Sex), mean, na.rm=TRUE)
  means <- means+5 # this is to enable subzero plotting on log transformed data
  SEMs <- tapply(tmp[, valueCol], list(tmp$Carriers, tmp$Sex), se)
  bp <- barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
                xlab="", ylab=ylab, mgp=c(axtitledist,0.5,0), xpd=FALSE)
  magaxis(side=2, las=2, mgp=c(3.0, 0.6, 0.0), labels=FALSE, unlog=TRUE)  # magaxis provides easy log ticks
  par(new=TRUE)  # enables replotting over the yaxis ticks, using the next line 
  barplot(means, beside=T, yaxt="n", col=cols, ylim=5+log10(ylim), 
          xlab="recipient sex", ylab="", mgp=c(axtitledist,0.5,0), xpd=FALSE)
  par(new=FALSE)
  axis(2, las=2, mgp=c(3,1.7,0), tck=-0.02, hadj=0, at=c(-2+5, -1+5, 0+5, 1+5, 2+5, 3+5, 4+5, 5+5),
       labels=c(expression(10^-2), expression(10^-1),expression(10^0),
                expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
  axis(1, las=2, tck=-0.02, at=c(2,5), labels=c("",""))
  title(main=title, line=0.5)
  arrows(bp, means+SEMs, bp, means-SEMs, lwd = 1.0, angle = 90, code = 3, length = 0.025)
  text(x=2, y=sig1+5, labels=starsF[, valueCol], cex=0.9) # add significance stars
  text(x=5, y=sig1+5, labels=starsM[, valueCol], cex=0.9) # add significance stars
  box()
}
#legend(locator(1),rownames(dat),fill=c("#ee7700","#3333ff"))

#pdf(file="./week20.pdf", width=11.5/2.54, height=8/2.54) #, family='Calibri')
png("figX_age_carrier_1col.png", width=(9.0*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=8)
par(mfrow=c(2,3), mar=c(3.0, 3.0, 2.1, 0.8), cex=0.7, mgp=c(2,0.6,0))
axtitledist <- 1.6  # adjusts distance of x and y axis titles 
sig1 <- log10(85)
sig2 <- log10(85)
cols3 <- c("#ee7700","#3333ff")  # colours for bar plots
cols4 <- c("orange","yellow")  # colours for bar plots
cols5 <- c("purple4","green4")  # colours for bar plots 
cols6 <- c("#1b9e77","#7570b3")  # colours for bar plots 

# BM
weeks <- c(20)

# Age
## stats
### females
dat <- BM34
dat <- dat[dat$Sex=="F", ]
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
lineages <- names(dat[8:ncol(dat)])
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp <- (t.test(dat[dat$Age=="young" & dat$Week==w, l], 
                   dat[dat$Age=="old" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out) 
stars <- p2stars(out)
# Workaround: adds 7 columns to the right so that my previously made valueCol script works: 
starsF <- cbind(matrix(NA,1,7), stars)
### males
dat <- BM34
dat <- dat[dat$Sex=="M", ]
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
lineages <- names(dat[8:ncol(dat)])
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp <- (t.test(dat[dat$Age=="young" & dat$Week==w, l], 
                   dat[dat$Age=="old" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out) 
stars <- p2stars(out)
# Workaround: adds 7 columns to the right so that my previously made valueCol script works: 
starsM <- cbind(matrix(NA,1,7), stars)
## plots
GroupBarplot2(BM34, week=weeks, valueCol=8,  cols=cols3, ylab=BMylab, ylim=c(0.1, 100), title="CD45")
GroupBarplot2(BM34, week=weeks, valueCol=9,  cols=cols3, ylab=BMylab, ylim=c(0.1, 100), title="GM")
GroupBarplot2(BM34, week=weeks, valueCol=10, cols=cols3, ylab=BMylab, ylim=c(0.1, 100), title="B Lymphoid")


# Carriers
## stats
### females
dat <- BM34
dat <- dat[dat$Sex=="F", ]
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
lineages <- names(dat[8:ncol(dat)])
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp <- (t.test(dat[dat$Carriers=="no" & dat$Week==w, l], 
                   dat[dat$Carriers=="yes" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out) 
stars <- p2stars(out)
# Workaround: adds 7 columns to the right so that my previously made valueCol script works: 
starsF <- cbind(matrix(NA,1,7), stars)
### males
dat <- BM34
dat <- dat[dat$Sex=="M", ]
dat[8:ncol(dat)] <- log10(dat[8:ncol(dat)])
lineages <- names(dat[8:ncol(dat)])
out <- NULL
for (l in lineages){
  lin <- NULL
  for (w in weeks){
    tmp <- (t.test(dat[dat$Carriers=="no" & dat$Week==w, l], 
                   dat[dat$Carriers=="yes" & dat$Week==w, l], var.equal = FALSE))
    lin <- c(lin, tmp$p.value)
  }
  out <- cbind(out, lin) 
}
colnames(out) <- lineages
rownames(out) <- weeks
print(out) 
stars <- p2stars(out)
# Workaround: adds 7 columns to the right so that my previously made valueCol script works: 
starsM <- cbind(matrix(NA,1,7), stars)
## plots
GroupBarplot4(BM34, week=weeks, valueCol=8,  cols=cols5, ylab=BMylab, ylim=c(0.1, 100), title="CD45")
GroupBarplot4(BM34, week=weeks, valueCol=9,  cols=cols5, ylab=BMylab, ylim=c(0.1, 100), title="GM")
GroupBarplot4(BM34, week=weeks, valueCol=10, cols=cols5, ylab=BMylab, ylim=c(0.1, 100), title="B Lymphoid")

dev.off()





# Radiation Survival Curve (of NRG)
library(survival)
data(package = "survival")
data(lung)


dat <- read.csv("../radiation_sensitivity/survival.csv")  #, colClasses=c(Mouse="character"))
# Add survival object and do calculation
dat$SurvObj <- with(dat, Surv(Survival, Censor == "N"))
tmp<- survfit(SurvObj ~ Dose, data = dat, conf.type = "log-log")
# Irradiation Doses
doses <- unique(dat$Dose)
# Make plot
png("figX_irradiation_1col.png", width=(9.0*ppi)/2.54, height=(8*ppi)/2.54, res=ppi, pointsize=8)
plot(tmp, yaxt="n", xaxt="n", ylim=c(0,1),
     xlab="weeks post-irradiation", ylab="fraction alive", 
     mgp=c(axtitledist+0.4,0.5,0), xpd=FALSE)
axis(2, las=2, mgp=c(0,1.7,0), tck=-0.02, hadj=0, at=seq(0,1,0.2) )
axis(side=1, at=seq(0,6,1),  mgp=c(0.4,0.4,0), cex=0.8, tck=-0.03)
xplace <- 5.5
toplab <- 1.02
step <- 0.06
units <- "cGy"
text(x=xplace, y=toplab-1*step, cex=0.9, labels=paste(doses[1],units))
text(x=xplace, y=toplab-2*step, cex=0.9, labels=paste(doses[2],units))
text(x=xplace, y=toplab-3*step, cex=0.9, labels=paste(doses[3],units))
text(x=xplace, y=toplab-4*step, cex=0.9, labels=paste(doses[5],units))
text(x=3.5, y=0.1, cex=0.9, labels=paste(doses[4],units))
dev.off()





