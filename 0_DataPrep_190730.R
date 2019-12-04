


#########################  Prepare inputs for analysis  ################################################



# Laptop
#setwd("C:/Users/JB/Box Sync/Projects/Projects (active)/Functional diversity/Analysis")
# Desktop
#setwd("D:/Box Sync/Projects/Projects (active)/Functional diversity/Analysis")
# Sara
#setwd("C:/Users/saraw/Desktop/FD")



library(dplyr)
#library(tidyr)
#library(sf)
#library(picante)
#library(raster)
#library(treeman)
library(ape)
library(phytools)


rm(list = ls())


#----------------------- LOAD AND CLEAN DATA -------------------------------------------------
Threats0 <- data.frame(read.csv("Data/Raw/Threats_mammals.csv", header=T))
Traits <- data.frame(read.csv("Data/Raw/Traits_mammals.csv", header=T)) # Wilman (2014) trait data (see "Articles/Measuring functional diversity" folder)
spptree0 <- read.tree("Data/Raw/Uyeda_etal_tetrapods.tree")
#  This tree was shared with J. Bordie privately and so has not been made publicly available here 


# Clean some of the data
names(Traits)[names(Traits) == 'Scientific'] <- 'binary'
Threats0$binary <- stringr::str_replace_all(Threats0$binary, c(" " = "_"))
Threats <- Threats0
dim(Threats) # 5674 species


# Remove marine species
Threats <- subset(Threats, order != "Sirenia")
Threats <- subset(Threats, binary != "Enhydra_lutris")
Threats <- subset(Threats, family != "Eschrichtiidae")
Threats <- subset(Threats, family != "Monodontidae")
Threats <- subset(Threats, family != "Iniidae")
Threats <- subset(Threats, family != "Delphinidae")
Threats <- subset(Threats, family != "Lipotidae")
Threats <- subset(Threats, family != "Phocoenidae")
Threats <- subset(Threats, family != "Pontoporiidae")
Threats <- subset(Threats, family != "Balaenidae")
Threats <- subset(Threats, family != "Balaenopteridae")
Threats <- subset(Threats, family != "Ziphiidae")
Threats <- subset(Threats, family != "Neobalaenidae")
Threats <- subset(Threats, family != "Physeteridae")
Threats <- subset(Threats, family != "Platanistidae")


# Add back in the river dolphins
riv01 <- subset(Threats0, binary == "Platanista_gangetica")
riv02 <- subset(Threats0, binary == "Inia_geoffrensis")
riv03 <- subset(Threats0, binary == "Inia_araguaiaensis")
riv04 <- subset(Threats0, binary == "Inia_boliviensis")
riv05 <- subset(Threats0, binary == "Pontoporia_blainvillei")
riv06 <- subset(Threats0, binary == "Lipotes_vexillifer")
rivdol <- rbind(riv01, riv02, riv03, riv04, riv05, riv06)
Threats <- rbind(Threats, rivdol)
dim(Threats) # 5583 species (removed 91 spp)


# Remove extinct species
Threats <- subset(Threats, risk != "Extinct")
dim(Threats) # 5501 species (removed 82 spp)



#----------------------- LINK TO SPECIES TRAITS -------------------------------------------

#----- Combine threats and traits
# Link the traits and threats databases
Dat0 <- suppressWarnings(left_join(Threats, Traits, by="binary"))


# Combine vertebrate feeding categories
Dat0$Vert <- Dat0$Diet.Vend + Dat0$Diet.Vect + Dat0$Diet.Vunk # vertebrate endotherms, ectotherms, & vertebrate unknown


# Foraging strata
Dat0$foraging <- Dat0$ForStrat.Value
levels(Dat0$foraging)[levels(Dat0$foraging)=="S"] <- "Ar" # lumping scansorial with arboreal
levels <- levels(Dat0$foraging)
levels[length(levels) + 1] <- "None"
Dat0$foraging <- factor(Dat0$foraging, levels = levels)
Dat0$foraging[is.na(Dat0$foraging)] <- "None"
levels(Dat0$foraging)[levels(Dat0$foraging)=="None"] <- "G"
Dat0$Ground <- ifelse(Dat0$foraging=="G", 1, 0)
Dat0$Climbing <- ifelse(Dat0$foraging=="Ar", 1, 0)
Dat0$Volant <- ifelse(Dat0$foraging=="A", 1, 0)


# Activity time
# Not currently used. It's not clear that activity time would strongly affect ecological function


# Data frame of species threats and traits
dat0 <- cbind.data.frame(rowID=Dat0$rowID, species=Dat0$binary, family=Dat0$family, order=Dat0$order,
	binary=Dat0$binary, trend=Dat0$trend, ExceptionsToDecline=Dat0$ExceptionsToDecline,
	habitat=Dat0$habitat, hunting=Dat0$hunting, conflict=Dat0$conflict, climate=Dat0$climate,
	NonNatives=Dat0$NonNatives, pollution=Dat0$pollution, hybrid=Dat0$hybrid, prey=Dat0$prey, 
	disease=Dat0$disease, inbreeding=Dat0$inbreeding,
	DietInv=Dat0$Diet.Inv, DietVert=Dat0$Vert, 
	DietFish=Dat0$Diet.Vfish, DietScav=Dat0$Diet.Scav, DietFruit=Dat0$Diet.Fruit, DietNect=Dat0$Diet.Nect,
	DietSeed=Dat0$Diet.Seed, DietHerb=Dat0$Diet.PlantO, BodyMass=Dat0$BodyMass.Value, Ground=Dat0$Ground,
	Climbing=Dat0$Climbing, Volant=Dat0$Volant)
dat0$species <- as.character(dat0$species)
dat0$DietInv <- as.numeric(dat0$DietInv)
dat0$DietVert <- as.numeric(dat0$DietVert)
dat0$DietFish <- as.numeric(dat0$DietFish)
dat0$DietScav <- as.numeric(dat0$DietScav)
dat0$DietFruit <- as.numeric(dat0$DietFruit)
dat0$DietNect <- as.numeric(dat0$DietNect)
dat0$DietSeed <- as.numeric(dat0$DietSeed)
dat0$DietHerb <- as.numeric(dat0$DietHerb)
dat1 <- cbind.data.frame(rowID=dat0$rowID, species=dat0$species, family=dat0$family, order=dat0$order,
	trend=dat0$trend, ExceptionsToDecline=dat0$ExceptionsToDecline,
	habitat=dat0$habitat, hunting=dat0$hunting, conflict=dat0$conflict,
	climate=dat0$climate, NonNatives=dat0$NonNatives, pollution=dat0$pollution, hybrid=dat0$hybrid,
	prey=dat0$prey, disease=dat0$disease, inbreeding=dat0$inbreeding, DietInv=dat0$DietInv, 
	DietVert=dat0$DietVert, DietFish=dat0$DietFish, DietScav=dat0$DietScav, DietFruit=dat0$DietFruit,
	DietNect=dat0$DietNect, DietSeed=dat0$DietSeed, DietHerb=dat0$DietHerb, BodyMass=dat0$BodyMass)
dat1$Ground <- dat0$Ground
dat1$Climbing <- dat0$Climbing
dat1$Volant <- dat0$Volant



#----- Deal with missing traits
dim(dat1) # We're at 5501 species


# Add a "genus" column to dat1
tmp1 <- data.frame(species = dat1$species)
tmp1$species <- as.character(tmp1$species)
tmp2 <- t(data.frame(strsplit(tmp1$species, "_")))
tmp2 <- data.frame(tmp2)
names(tmp2) <- c("genus", "epithet")
dat1$genus <- tmp2$genus


# Find out which species are missing trait data
MissingTraits <- dat1[rowSums(is.na(dat1)) > 0,] # rows of dat1 that are missing trait data
dim(MissingTraits) # Missing trait data for 745 species


# Remove those species from dat1
dat2 <- dplyr::anti_join(dat1, MissingTraits, by="species") 


# Use genus-average trait data for each missing species
NewTraits <- MissingTraits
for(i in 1:nrow(MissingTraits)){
	#i=745
	dattmp <- subset(dat2, genus==NewTraits$genus[i])
	NewTraits$DietInv[i] <- mean(dattmp$DietInv, na.omit=T)
	NewTraits$DietVert[i] <- mean(dattmp$DietVert, na.omit=T)
	NewTraits$DietFish[i] <- mean(dattmp$DietFish, na.omit=T)
	NewTraits$DietScav[i] <- mean(dattmp$DietScav, na.omit=T)
	NewTraits$DietFruit[i] <- mean(dattmp$DietFruit, na.omit=T)
	NewTraits$DietNect[i] <- mean(dattmp$DietNect, na.omit=T)
	NewTraits$DietSeed[i] <- mean(dattmp$DietSeed, na.omit=T)
	NewTraits$DietHerb[i] <- mean(dattmp$DietHerb, na.omit=T)
	NewTraits$BodyMass[i] <- mean(dattmp$BodyMass, na.omit=T)
	NewTraits$Ground[i] <- max(dattmp$Ground, na.rm=T)
	NewTraits$Climbing[i] <- max(dattmp$Climbing, na.rm=T)
	NewTraits$Volant[i] <- max(dattmp$Volant, na.rm=T)	}
summary(NewTraits) # 212 NA's
# Genus-level averaging worked for 745-212=533 species
# Still missing traits for 212 species. 


# Use family-level averaging for species that are still missing trait data
StillMissingTraits <- NewTraits[rowSums(is.na(NewTraits)) > 0,] # rows of NewTraits that are missing trait data
NewTraits <- dplyr::anti_join(NewTraits, StillMissingTraits, by="species") # remove those species from NewTraits 
NewTraits2 <- StillMissingTraits
for(i in 1:nrow(StillMissingTraits)){
	#i=1
	dattmp <- subset(dat2, family==NewTraits2$family[i])
	NewTraits2$DietInv[i] <- mean(dattmp$DietInv, na.omit=T)
	NewTraits2$DietVert[i] <- mean(dattmp$DietVert, na.omit=T)
	NewTraits2$DietFish[i] <- mean(dattmp$DietFish, na.omit=T)
	NewTraits2$DietScav[i] <- mean(dattmp$DietScav, na.omit=T)
	NewTraits2$DietFruit[i] <- mean(dattmp$DietFruit, na.omit=T)
	NewTraits2$DietNect[i] <- mean(dattmp$DietNect, na.omit=T)
	NewTraits2$DietSeed[i] <- mean(dattmp$DietSeed, na.omit=T)
	NewTraits2$DietHerb[i] <- mean(dattmp$DietHerb, na.omit=T)
	NewTraits2$BodyMass[i] <- mean(dattmp$BodyMass, na.omit=T)
	NewTraits2$Ground[i] <- max(dattmp$Ground, na.rm=T)
	NewTraits2$Climbing[i] <- max(dattmp$Climbing, na.rm=T)
	NewTraits2$Volant[i] <- max(dattmp$Volant, na.rm=T)	}
summary(NewTraits2) # 1 NA
# Genus-level averaging worked for 212-1=211 species
# Still missing traits for 1 species. 


# Find species that is *still* missing trait data
StillStillMissing <- NewTraits2[rowSums(is.na(NewTraits2)) > 0,] 
NewTraits2 <- dplyr::anti_join(NewTraits2, StillStillMissing, by="species") # remove those species from NewTraits 
NewTraits3 <- StillStillMissing
NewTraits3$species # Just one critter still missing: Laonastes_aenigmamus!


# Enter data for Laonastes manually, based on IUCN Red List (accessed 6 May 2019)
NewTraits3$DietInv <- 10
NewTraits3$DietVert <- 0
NewTraits3$DietFish <- 0
NewTraits3$DietScav <- 0
NewTraits3$DietFruit <- 0
NewTraits3$DietNect <- 0
NewTraits3$DietSeed <- 30
NewTraits3$DietHerb <- 60
rats <- subset(dat1, genus == "Rattus")
NewTraits3$BodyMass <- mean(rats$BodyMass, na.rm=T)
NewTraits3$Ground <- 1
NewTraits3$Climbing <- 1
NewTraits3$Volant <- 0


# Assemble the new data
dat3 <- rbind.data.frame(dat2, NewTraits, NewTraits2, NewTraits3)
# Now we have trait data for all 5501 species
#----- End of section on 'dealing with missing traits'


#----- Standardize data for continuous traits
dat1 <- dat3
dat1$DietInv <- as.numeric(scale(dat1$DietInv, scale=TRUE, center=TRUE))
dat1$DietVert <- as.numeric(scale(dat1$DietVert, scale=TRUE, center=TRUE))
dat1$DietFish <- as.numeric(scale(dat1$DietFish, scale=TRUE, center=TRUE))
dat1$DietScav <- as.numeric(scale(dat1$DietScav, scale=TRUE, center=TRUE))
dat1$DietFruit <- as.numeric(scale(dat1$DietFruit, scale=TRUE, center=TRUE))
dat1$DietNect <- as.numeric(scale(dat1$DietNect, scale=TRUE, center=TRUE))
dat1$DietSeed <- as.numeric(scale(dat1$DietSeed, scale=TRUE, center=TRUE))
dat1$DietHerb <- as.numeric(scale(dat1$DietHerb, scale=TRUE, center=TRUE))
dat1$BodyMass <- as.numeric(scale(dat1$BodyMass, scale=TRUE, center=TRUE))


#----- Save file
dat1 <- subset(dat1, select = -c(genus))
write.csv(dat1, "Data/Raw/mammal_threats_traits.csv", row.names=F)
#  This file is then used as input in the main processing script: FD_PD_190721.R








#----------------- PHYLOGENY -------------------------------------

ape::is.ultrametric(spptree0)
ape::is.rooted(spptree0)
ape::is.binary.tree(spptree0)


#---- Remove genera from the phylogeny that aren't in the threats database
#---- (this massively cuts down the size of the phylogeny file, while keeping genera we'll need later)

# add a "genus" column back to dat1
tmp1 <- data.frame(species = dat1$species)
tmp1$species <- as.character(tmp1$species)
tmp2 <- t(data.frame(strsplit(tmp1$species, "_")))
tmp2 <- data.frame(tmp2)
names(tmp2) <- c("genus", "epithet")
dat1$genus <- tmp2$genus


# genera in the phylogeny
tmp <- spptree0$tip.label
test <- sub("_[^_]+$", "", tmp)
test <- sub("_[^_]+$", "", test)
gentree <- data.frame(test)
gentree <- data.frame(genus=unique(gentree[,1]))
gentree$genus <- as.character(gentree$genus)


# list of genera to remove from phylogeny
genremove <- suppressWarnings(dplyr::anti_join(gentree, dat1, by = "genus")) 


# remove those genera from the phylogeny
spptree <- spptree0
for(i in 1:nrow(genremove)){
	spptree <- drop.tip(spptree, tip = grep(genremove[i,1], spptree$tip.label, value=T))	}


# check that that worked
tmp <- spptree$tip.label
test <- sub("_[^_]+$", "", tmp)
test <- sub("_[^_]+$", "", test)
gentree2 <- data.frame(test)
gentree2 <- data.frame(genus = unique(gentree2[,1]))
gentree2$genus <- as.character(gentree2$genus)
genremove2 <- suppressWarnings(dplyr::anti_join(gentree2, dat1, by = "genus"))
nrow(genremove2) # should be 0


# find out which species are in the threats/traits database but not the phylogeny
speciesthreats <- data.frame(species = unique(dat1$species))
speciesthreats$species <- as.character(speciesthreats$species)
speciestree <- data.frame(species = spptree$tip.label)
speciestree$species <- as.character(speciestree$species)
SppMissing <- suppressWarnings(dplyr::anti_join(speciesthreats, speciestree, by = "species"))
NumSppNotInPhylo <- nrow(SppMissing)
NumSppNotInPhylo # 1664 species


# add a 'genus' column to SppMissing
tmp <- SppMissing$species
test <- sub("_[^_]+$", "", tmp)
gen <- data.frame(test)
SppMissing$genus <- test


# remove species from SppMissing that don't have a genus in the phylogeny
SppMissing2 <- suppressWarnings(dplyr::inner_join(SppMissing, gentree2, by = "genus"))
nrow(SppMissing2) # able to add 1383 of the missing species to the phylogeny at the root node of their genus
NumSppCantBeAddedToPhylo <- NumSppNotInPhylo - nrow(SppMissing2)
NumSppCantBeAddedToPhylo # 281 species from the threats data can't be added to the phylogeny



#---- Make the tree ultrametric
#l <- seq(0, 100, by=1)
#LL.out <- NULL
#start <- Sys.time()
#for(i in 1:length(l)){
#	LL.out[i] <- attributes(chronos(spptree, lambda=l[i]))$ploglik	}
#Sys.time() - start
#write.csv(LL.out, "LLout.csv", row.names=F)


#--- Find a good lambda value
LL.out <- matrix(, nrow=6, ncol=2)
LL.out <- data.frame(LL.out)
names(LL.out) <- c("lambda", "loglik")
LL.out$lambda <- c(0, 0.5, 1, 2, 0.25, 0.75)

start <- Sys.time()
LL.out$loglik[1] <- attributes(chronos(spptree, lambda = 0))$ploglik 
Sys.time() - start
write.csv(LL.out, "LLout.csv", row.names=F)

start <- Sys.time()
LL.out$loglik[2] <- attributes(chronos(spptree, lambda = 0.5))$ploglik 
Sys.time() - start
write.csv(LL.out, "LLout.csv", row.names=F)

start <- Sys.time()
LL.out$loglik[3] <- attributes(chronos(spptree, lambda = 1.0))$ploglik 
Sys.time() - start
write.csv(LL.out, "LLout.csv", row.names=F)

start <- Sys.time()
LL.out$loglik[4] <- attributes(chronos(spptree, lambda = 2.0))$ploglik 
Sys.time() - start
write.csv(LL.out, "LLout.csv", row.names=F)

start <- Sys.time()
LL.out$loglik[5] <- attributes(chronos(spptree, lambda = 0.25))$ploglik 
Sys.time() - start
write.csv(LL.out, "LLout.csv", row.names=F)

start <- Sys.time()
LL.out$loglik[6] <- attributes(chronos(spptree, lambda = 0.75))$ploglik 
Sys.time() - start
write.csv(LL.out, "LLout.csv", row.names=F)


#--- Make the phylogeny ultrametric using best lambda value
# log likelihoods decline linearly and precipitously from 0 to 2.
spptree1 <- ape::chronos(spptree, lambda = 0) 


ape::is.ultrametric(spptree1)




#-------- Add species to the phylogeny (at the root node of the appropriate genus) that are in the threats database but not the phylogeny

spptree2 <- spptree1
start <- Sys.time()
for(i in 1:nrow(SppMissing2)){
   tryCatch({
	#i=9
	spptree2 <- phytools::add.species.to.genus(spptree2, SppMissing2[i,1], where = "root")
   }, error=function(e){}) # end of the tryCatch function
}
Sys.time() - start


# were we able to add all missing spp (that had genera in the phylogeny) to the phylogeny?
nrow(SppMissing2) == length(spptree2$tip.label) - length(spptree1$tip.label)
# if 'TRUE', we've added them all


# save phylogeny
write.tree(spptree2, file = "mammaltree.tree")






































