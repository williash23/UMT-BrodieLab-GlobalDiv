

library(ggplot2)
library(dplyr)
library(sf)
library(sp)
library(cowplot)


# desktop
#setwd('D:/Box Sync/Projects/Projects (active)/Functional diversity/Analysis')
# laptop
#setwd('C:/Users/JB/Box Sync/Projects/Projects (active)/Functional diversity/Analysis')



############## DATA ######################################################################
rm(list = ls())
dat0 <- data.frame(read.csv("FuncPhylo_mammals190816.csv", header=T))
#  Output from main analysis


############## SOME SIMPLE CORRELATIONS ########################################################
dat1 <- cbind.data.frame(SR.0=dat0$SR0, FD.0=dat0$FD0, PD.0=dat0$PD0)
dat1 <- dat1[complete.cases(dat1),]
cor(dat1)
cor.test(dat1$SR.0, dat1$FD.0)
cor.test(dat1$SR.0, dat1$PD.0)
cor.test(dat1$FD.0, dat1$PD.0)



############## HOW THE DIFFERENT THREATS CHANGE FD AND PD ########################################
#----------- CREATES THE BOX PLOTS USED IN FIG 1 OF THE MANUSCRIPT ---------------------------------


#------- Change in FD
dat1 <- cbind.data.frame(FD.change.habitat = dat0$FD.habitat1-1, FD.change.hunting = dat0$FD.hunting1-1,
	FD.change.climate = dat0$FD.climate1-1, FD.change.pollution = dat0$FD.pollution1-1, FD.change.prey = dat0$FD.prey1-1, FD.change.conflict = dat0$FD.conflict1-1, 
	FD.change.nonnatives = dat0$FD.nonNatives1-1, FD.change.hybrid = dat0$FD.hybrid1-1, FD.change.disease = dat0$FD.disease1-1, FD.change.inbreeding = dat0$FD.inbreeding1-1)
datstack <- stack(dat1)
names(datstack) <- c("diversityloss", "threat")

# one-way ANOVA
resaov <- aov(diversityloss ~ threat, data=datstack) 
summary(resaov)
TukeyHSD(resaov)

# box plot
plotfd <- ggplot(datstack, aes(x=threat, y=diversityloss)) + geom_boxplot(outlier.shape=NA, fill='#A4A4A4', color="black") + 
	ylim(-0.25, 0) + 
	theme_classic() +
	geom_hline(yintercept=0, linetype="dashed", color="black") + 
	geom_point(aes())

# how many data points per category?
dattmp <- na.omit(datstack) 
summary(dattmp)


#------- Change in PD
dat1 <- cbind.data.frame(PD.change.habitat = dat0$PD.habitat1-1, PD.change.hunting = dat0$PD.hunting1-1,
	PD.change.climate = dat0$PD.climate1-1, PD.change.pollution = dat0$PD.pollution1-1, PD.change.prey = dat0$PD.prey1-1, PD.change.conflict = dat0$PD.conflict1-1, 
	PD.change.nonnatives = dat0$PD.nonNatives1-1, PD.change.hybrid = dat0$PD.hybrid1-1, PD.change.disease = dat0$PD.disease1-1, PD.change.inbreeding = dat0$PD.inbreeding1-1)
datstack <- stack(dat1)
names(datstack) <- c("diversityloss", "threat")

# one-way ANOVA
resaov <- aov(diversityloss ~ threat, data=datstack) 
summary(resaov)
TukeyHSD(resaov)

# box plot
plotpd <- ggplot(datstack, aes(x=threat, y=diversityloss)) + geom_boxplot(outlier.shape=NA, fill='#A4A4A4', color="black") + 
	ylim(-0.25, 0) + 
	theme_classic() +
	geom_hline(yintercept=0, linetype="dashed", color="black") + 
	geom_point(aes())


# how many data points per category?
dattmp <- na.omit(datstack) 
summary(dattmp)


#------ Multi-panel figure
plot_grid(plotfd, plotpd, labels = "AUTO", ncol = 1, nrow = 2, label_size = 11)
ggsave("Results/Graphs/ThreatImpacts_A_fd_B_pd.jpg")


#------ Some associated stats
# "...on average, the impacts of habitat loss and harvest... exceed those of climate change by..."

# factor by which habitat loss impacts on FD exceed those of climate change
(1 - mean(dat0$FD.habitat1, na.rm=T)) / (1 - mean(dat0$FD.climate1, na.rm=T))

# factor by which harvest impacts on FD exceed those of climate change
(1 - mean(dat0$FD.hunting1, na.rm=T)) / (1 - mean(dat0$FD.climate1, na.rm=T))

# factor by which hunting impacts on PD exceed those of climate change
(1 - mean(dat0$PD.habitat1, na.rm=T)) / (1 - mean(dat0$PD.climate1, na.rm=T))

# impacts on PD do not differ significantly between climate and harvest







############## CHANGES IN WHICH TRAITS ARE DRIVING DECLINES IN FD? ########################################
# "...Across countries, reductions in mammal functional diversity driven by habitat loss and harvest were correlated with nationwide 
# declines in frugivores..."

#--- habitat
DivChange <- (1 - dat0$FD.habitat1)
traitcols <- seq(149, 160, by=1) # columns of dat0 for trait changes
dat1 <- dat0[traitcols]
dat1 <- cbind.data.frame(DivChange, dat1)
dat1 <- dat1[complete.cases(dat1),]
corrs <- cor(dat1)
corrs[,1]

#--- harvest
DivChange <- (1 - dat0$FD.hunting1)
traitcols <- seq(161, 172, by=1) # columns of dat0 for trait changes
dat1 <- dat0[traitcols]
dat1 <- cbind.data.frame(DivChange, dat1)
dat1 <- dat1[complete.cases(dat1),]
corrs <- cor(dat1)
corrs[,1]






############## COMPARE PROPORTIONAL LOSS RATES OF SR, FD, & PD ######################################################################
# "...average losses of taxonomic, functional, and phylogenetic diversity across the world from all threats combined were..."


#----------------- FD vs SR ----------------------------------------

#--- habitat ---
dat1 <- dat0
dat1$SR.loss.prop <- (dat1$SR0 - dat1$SR.habitat1) / dat1$SR0
dat1$loss.prop <- (1 - dat1$FD.habitat1)
dat2 <- cbind.data.frame(country = dat1$Country, area = dat1$Area_sq_km, SR0 = dat1$SR0, 
	SR.loss.prop = dat1$SR.loss.prop, loss.prop = dat1$loss.prop)
dat2$diff <- dat2$SR.loss.prop - dat2$loss.prop
#dat2 <- subset(dat2, area >= 50000)

# paired t test
#dat3 <- cbind.data.frame(SR.loss.prop = dat2$SR.loss.prop, loss.prop = dat2$loss.prop)
#t.test(dat3$SR.loss.prop, dat3$loss.prop, paired = TRUE, alternative = "two.sided")

# regression of standing diversity vs the diff between SR.loss.prop and loss.prop
m1 <- glm(diff ~ SR0, data=dat2)
summary(m1)

# regression of countries' areas vs the diff between SR.loss.prop and loss.prop
m2 <- glm(diff ~ area, data=dat2)
summary(m2)

# plotting
newdat <- data.frame(SR0 = seq(50, 700, by=50))
newdat$Ypred <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$fit
newdat$YpredSE <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$se.fit
newdat$YpredLO <- newdat$Ypred - (1.96 * newdat$YpredSE)
newdat$YpredHI <- newdat$Ypred + (1.96 * newdat$YpredSE)
newdat


#--- hunting ---
dat1 <- dat0
dat1$SR.loss.prop <- (dat1$SR0 - dat1$SR.hunting1) / dat1$SR0
dat1$loss.prop <- (1 - dat1$FD.hunting1)
dat2 <- cbind.data.frame(country = dat1$Country, area = dat1$Area_sq_km, SR0 = dat1$SR0, 
	SR.loss.prop = dat1$SR.loss.prop, loss.prop = dat1$loss.prop)
dat2$diff <- dat2$SR.loss.prop - dat2$loss.prop
#dat2 <- subset(dat2, area >= 50000)

# paired t test
dat3 <- cbind.data.frame(SR.loss.prop = dat2$SR.loss.prop, loss.prop = dat2$loss.prop)
t.test(dat3$SR.loss.prop, dat3$loss.prop, paired = TRUE, alternative = "two.sided")

# regression of standing diversity vs the diff between SR.loss.prop and loss.prop
m1 <- glm(diff ~ SR0, data=dat2)
summary(m1)

# regression of countries' areas vs the diff between SR.loss.prop and loss.prop
m2 <- glm(diff ~ area, data=dat2)
summary(m2)

# plotting
newdat <- data.frame(SR0 = seq(50, 700, by=50))
newdat$Ypred <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$fit
newdat$YpredSE <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$se.fit
newdat$YpredLO <- newdat$Ypred - (1.96 * newdat$YpredSE)
newdat$YpredHI <- newdat$Ypred + (1.96 * newdat$YpredSE)
newdat



#----------------- PD vs SR----------------------------------------

#--- habitat ---
dat1 <- dat0
dat1$SR.loss.prop <- (dat1$SR0 - dat1$SR.habitat1) / dat1$SR0
dat1$loss.prop <- (1 - dat1$PD.habitat1)
dat2 <- cbind.data.frame(country = dat1$Country, area = dat1$Area_sq_km, SR0 = dat1$SR0, 
	SR.loss.prop = dat1$SR.loss.prop, loss.prop = dat1$loss.prop)
dat2$diff <- dat2$SR.loss.prop - dat2$loss.prop
#dat2 <- subset(dat2, area >= 50000)

# paired t test
#dat3 <- cbind.data.frame(SR.loss.prop = dat2$SR.loss.prop, loss.prop = dat2$loss.prop)
#t.test(dat3$SR.loss.prop, dat3$loss.prop, paired = TRUE, alternative = "two.sided")

# regression of standing diversity vs the diff between SR.loss.prop and loss.prop
m1 <- glm(diff ~ SR0, data=dat2)
summary(m1)

# regression of countries' areas vs the diff between SR.loss.prop and loss.prop
m2 <- glm(diff ~ area, data=dat2)
summary(m2)

# plotting
newdat <- data.frame(SR0 = seq(50, 700, by=50))
newdat$Ypred <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$fit
newdat$YpredSE <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$se.fit
newdat$YpredLO <- newdat$Ypred - (1.96 * newdat$YpredSE)
newdat$YpredHI <- newdat$Ypred + (1.96 * newdat$YpredSE)
newdat


#--- hunting ---
dat1 <- dat0
dat1$SR.loss.prop <- (dat1$SR0 - dat1$SR.hunting1) / dat1$SR0
dat1$loss.prop <- (1 - dat1$PD.hunting1)
dat2 <- cbind.data.frame(country = dat1$Country, area = dat1$Area_sq_km, SR0 = dat1$SR0, 
	SR.loss.prop = dat1$SR.loss.prop, loss.prop = dat1$loss.prop)
dat2$diff <- dat2$SR.loss.prop - dat2$loss.prop
#dat2 <- subset(dat2, area >= 50000)

# paired t test
#dat3 <- cbind.data.frame(SR.loss.prop = dat2$SR.loss.prop, loss.prop = dat2$loss.prop)
#t.test(dat3$SR.loss.prop, dat3$loss.prop, paired = TRUE, alternative = "two.sided")

# regression of standing diversity vs the diff between SR.loss.prop and loss.prop
m1 <- glm(diff ~ SR0, data=dat2)
summary(m1)

# regression of countries' areas vs the diff between SR.loss.prop and loss.prop
m2 <- glm(diff ~ area, data=dat2)
summary(m2)

# plotting
newdat <- data.frame(SR0 = seq(50, 700, by=50))
newdat$Ypred <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$fit
newdat$YpredSE <- predict(m1, newdata=newdat, type="response", se.fit=T, na.omit=T)$se.fit
newdat$YpredLO <- newdat$Ypred - (1.96 * newdat$YpredSE)
newdat$YpredHI <- newdat$Ypred + (1.96 * newdat$YpredSE)
newdat







