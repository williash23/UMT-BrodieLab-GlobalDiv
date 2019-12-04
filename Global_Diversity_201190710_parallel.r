


#setwd("D:/Box Sync/Projects/Projects (active)/Functional diversity/Analysis")
#setwd("C:/Users/sw223723e/Desktop/FD")

library(dplyr)
library(tidyr)
library(sf)
library(picante)
#library(raster)
library(ape)


rm(list = ls())



# ----------------------- USER-DEFINED PARAMETERS ----------------------------------------------------------------
# Load functions
source("FunctionalDiversity190528_funs.R")

st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))

# Desired patial projection
prj <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
	
# Number of randomization iterations to calculate bias in funct & phylogen loss
numrands <- 100

# Parameters for Petchy's functional diversity calculation
Distance.method <- "euclidean"
Cluster.method <- "average"


#----------------------- MAMMAL DATA -------------------------------------------------
Threats <- data.frame(read.csv("Data/Raw/Threats_mammals.csv", header=T))
Threats$ExceptionsToDecline[Threats$ExceptionsToDecline == "na"] <- NA
Traits <- data.frame(read.csv("Data/Raw/Function_mammals.csv", header=T)) # Wilman (2014) trait data (see "Articles" folder)
mammtree <- read.tree("Data/Raw/mammaltree.tree")
spptree <- mammtree

Threats <- subset(Threats, risk != "Extinct")

# Clean some of the data
names(Traits)[names(Traits) == 'Scientific'] <- 'binary'
Threats$binary <- stringr::str_replace_all(Threats$binary, c(" " = "_"))
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
#Not currently used. It's not clear that activity time would strongly affect ecological function

# Data frame of species threats and traits
dat0 <- cbind.data.frame(rowID=Dat0$rowID, species=Dat0$binary, family=Dat0$family, order=Dat0$order,
	binary=Dat0$binary, trend=Dat0$trend, 
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
	trend=dat0$trend, habitat=dat0$habitat, hunting=dat0$hunting, conflict=dat0$conflict,
	climate=dat0$climate, NonNatives=dat0$NonNatives, pollution=dat0$pollution, hybrid=dat0$hybrid,
	prey=dat0$prey, disease=dat0$disease, inbreeding=dat0$inbreeding)

# Standardize data for continuous traits
dat1$DietInv <- as.numeric(scale(dat0$DietInv, scale=TRUE, center=TRUE))
dat1$DietVert <- as.numeric(scale(dat0$DietVert, scale=TRUE, center=TRUE))
dat1$DietFish <- as.numeric(scale(dat0$DietFish, scale=TRUE, center=TRUE))
dat1$DietScav <- as.numeric(scale(dat0$DietScav, scale=TRUE, center=TRUE))
dat1$DietFruit <- as.numeric(scale(dat0$DietFruit, scale=TRUE, center=TRUE))
dat1$DietNect <- as.numeric(scale(dat0$DietNect, scale=TRUE, center=TRUE))
dat1$DietSeed <- as.numeric(scale(dat0$DietSeed, scale=TRUE, center=TRUE))
dat1$DietHerb <- as.numeric(scale(dat0$DietHerb, scale=TRUE, center=TRUE))
dat1$BodyMass <- as.numeric(scale(dat0$BodyMass, scale=TRUE, center=TRUE))
dat1$Ground <- dat0$Ground
dat1$Climbing <- dat0$Climbing
dat1$Volant <- dat0$Volant

#  Convert species names to character 
dat1$species <- as.character(dat1$species)



########################################################################################################################
### REDO THIS PART SO THAT WE GET A LIST OF SPECIES LISTS PER *COUNTRY* INSTEAD OF PER GRID CELL

#----------------------- COUNTRIES AND SPP RANGE PREP ---------------------------------

ctrs <- st_read("Data/Raw/UIA_World_Countries_Boundaries.shp") %>%
	dplyr::select(Country, geometry) %>%
	st_transform(st_crs(prj)) %>%
	mutate(ctry_sq_km = as.numeric(st_area(.))/1000000)  
ctrs_area <- ctrs %>%
	as.data.frame() %>%
	dplyr::select(Country, ctry_sq_km, - geometry) 
ctrs_df <- ctrs %>%	
	as.data.frame() %>%
	dplyr::select(-geometry) 
ctrs_df$Country <- as.character(ctrs_df$Country)


#  Run functions to prepare range polygons
#  Now removes mammals who are listed as extinct
mam <- prep_range_polys(fn = "Data/Raw/TERRESTRIAL_MAMMALS.shp", ctrs = ctrs)


#----------------------- SPECIES - COUNTRIES INTERSECTION ---------------------------------
polys_ctrs_int <- st_intersects(mam, ctrs, sparse = FALSE) 
#  Dimensions: nrow = number of spp (12112 for mammals), ncol = number of countries (254)

spp_int_ls <- list()

for(i in 1:nrow(ctrs_df)){
	polys_ctrs_int_vec <- which(polys_ctrs_int[,i] == TRUE)
	spp_list_ctr <- mam[polys_ctrs_int_vec, 1]
	st_geometry(spp_list_ctr) <- NULL
	spp_int_ls[[i]] <- spp_list_ctr	
	}
	

#----------------------- EXCEPTIONS TO DECLINE ---------------------------------
exc_tmp <- Dat0 %>%
	dplyr::select(binary, ExceptionsToDecline)
exc_tmp$ExceptionsToDecline <- as.character(exc_tmp$ExceptionsToDecline)


	spp <- exc_tmp$binary[1]
	tmp <- exc_tmp$ExceptionsToDecline[1]
	tmp_df <- as.data.frame(tmp) %>%
		tidyr::separate(tmp, paste("Country", 1:22, sep="_"), sep = ",")
	exc_df <- cbind(spp, tmp_df)

for(i in 1:nrow(exc_tmp)){
	spp <- exc_tmp$binary[i]
	tmp <- exc_tmp$ExceptionsToDecline[i]
	tmp_df <- as.data.frame(tmp) %>%
		tidyr::separate(tmp, paste("Country", 1:22, sep="_"), sep = ",")
	new_row <- cbind(spp, tmp_df)
	exc_df <- rbind(exc_df, new_row)
	}
	
exc_df <- exc_df[rowSums(is.na(exc_df[,2:23])) != 22,]
exc_df$spp <- as.character(exc_df$spp)
	
	
#----------------------- SAVE DATA FOR LATER USE ---------------------------------
#save(spp_int_ls, file="Data/spp_int_ls.Rdata")
#save(dat1, file="Data/dat1.Rdata")
#save(ctrs, file="Data/ctrs.Rdata")
#save(ctrs_df, file="Data/ctrs_df.Rdata")



 
################################################################################################################################
### THEN HERE WE'LL START THE 'FOR' LOOP, RUNNING THROUGH EACH COUNTRY ONE BY ONE





#----------------------- ANALYSIS W/ PARALLEL PROCESSING ------------------------------------
library(doSNOW)
library(parallel)
library(foreach)

# Number of countries to analyze
numcountries <- length(spp_int_ls)

nThreads <- 7
cluster = makeCluster(nThreads, type = "SOCK", outfile="Log.txt")

registerDoSNOW(cluster)
getDoParWorkers()

start <- Sys.time()
start

##  Parallel processing loop
#out <- foreach(i = 1:numcountries, .packages=c("dplyr", "sf", "picante"), .combine=rbind, .verbose = TRUE) %dopar% {
out <- foreach(i = 1:numcountries, .packages=c("dplyr", "sf", "picante"), .errorhandling = c("pass"), .verbose = TRUE) %dopar% {


	#  Exceptions to decline for country i
	country_name <- as.character(ctrs_df[i,1])
	keep_exc <- which(apply(exc_df, 1, function(r) any(r %in% country_name)))
	spp_exc <- exc_df[keep_exc,1]
	spp_exc_df <- as.data.frame(spp_exc)
	colnames(spp_exc_df) <- "species"
	
	# Species list in country i
	spp_list <- spp_int_ls[[i]]
	spp_list <- unique(spp_list) 

	datXY <- dplyr::left_join(spp_list, dat1, by="species") 
	datXY <- datXY[complete.cases(datXY),]
	all_spp <- as.data.frame(datXY$species)
	colnames(all_spp)[1] <- "species"
	SR0 <- nrow(datXY) # Taxonomic diversity (species richness)

	#----- Current diversity
  
	if(SR0 < 2){
		FD0 <- NA; FD_r <- rep(NA, 20); 
		FDrand_hab1 <- NA; FDrand_hab2 <- NA; FDrand_hunt1 <- NA; FDrand_hunt2 <- NA; FDrand_clim1 <- NA; 
		FDrand_clim2 <- NA; FDrand_con1 <- NA; FDrand_con2 <- NA; FDrand_nonnat1 <- NA; FDrand_nonnat2 <- NA; 
		FDrand_pol1 <- NA; FDrand_pol2 <- NA; FDrand_hyb1 <- NA; FDrand_hyb2 <- NA; FDrand_prey1 <- NA; 
		FDrand_prey2 <- NA; FDrand_dis1 <- NA; FDrand_dis2 <- NA; FDrand_inb1 <- NA; FDrand_inb2 <- NA;
		PD0 <- NA; PD_hab1 <- NA; PD_hab2 <- NA; PD_hunt1 <- NA; PD_hunt2 <- NA; 
		PD_clim1 <- NA; PD_clim2 <- NA; PD_con1 <- NA; PD_con2 <- NA; PD_nonnat1 <- NA; 
		PD_nonnat2 <- NA; PD_pol1 <- NA; PD_pol2 <- NA; PD_hyb1 <- NA; PD_hyb2 <- NA; PD_prey1 <- NA; 
		PD_prey2 <- NA; PD_dis1 <- NA; PD_dis2 <- NA; PD_inb1 <- NA; PD_inb2 <- NA;
		PDrand_hab1 <- NA; PDrand_hab2 <- NA; PDrand_hunt1 <- NA; PDrand_hunt2 <- NA; PDrand_clim1 <- NA; 
		PDrand_clim2 <- NA; PDrand_con1 <- NA; PDrand_con2 <- NA; PDrand_nonnat1 <- NA; PDrand_nonnat2 <- NA; 
		PDrand_pol1 <- NA; PDrand_pol2 <- NA; PDrand_hyb1 <- NA; PDrand_hyb2 <- NA; PDrand_prey1 <- NA; 
		PDrand_prey2 <- NA; PDrand_dis1 <- NA; PDrand_dis2 <- NA; PDrand_inb1 <- NA; PDrand_inb2 <- NA;

		hab1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & habitat == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hab1 <- nrow(hab1_spp)
	    
		hab2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & habitat == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & habitat == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hab2 <- nrow(hab2_spp)

		hunt1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & hunting == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hunt1 <- nrow(hunt1_spp)
		
		hunt2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & hunting == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & hunting == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hunt2 <- nrow(hunt2_spp)
		
		con1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & conflict == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_con1 <- nrow(con1_spp)
		
		con2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & conflict == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & conflict == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_con2 <- nrow(con2_spp)
	  
		clim1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & climate == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_clim1 <- nrow(clim1_spp)

		clim2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & climate == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & climate == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_clim2 <- nrow(clim2_spp)

		nonnat1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & NonNatives == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_nonnat1 <- nrow(nonnat1_spp)
		
		nonnat2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & NonNatives == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & NonNatives == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_nonnat2 <- nrow(nonnat2_spp)

		pol1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & pollution == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_pol1 <- nrow(pol1_spp)
	  
		pol2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & pollution == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & pollution == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_pol2 <- nrow(pol2_spp)
	  
		hyb1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & hybrid == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hyb1 <- nrow(hyb1_spp)
	  
		hyb2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & hybrid == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & hybrid == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hyb2 <- nrow(hyb2_spp)
		
		prey1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & prey == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_prey1 <- nrow(prey1_spp)
		
		prey2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & prey == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & prey == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_prey2 <- nrow(prey2_spp)
	  
		dis1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & disease == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_dis1 <- nrow(dis1_spp)
	  
		dis2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & disease == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & disease == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_dis2 <- nrow(dis2_spp)
		
		inb1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & inbreeding == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_inb1 <- nrow(inb1_spp)
	  
		inb2_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & inbreeding == 1)) %>%
			dplyr::filter(!(trend == "Decreasing" & inbreeding == 2)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_inb2 <- nrow(inb2_spp)
		
		} else {
			
			#  Functional diversity
			# Trait dendrogram for all species at the grid point
			####CHANGE: DietVend=datXY$DietVend no longer in dat1
			species.traits <- cbind.data.frame(species=datXY$species, DietInv=datXY$DietInv, 
			   DietVert=datXY$DietVert, DietFish=datXY$DietFish, 
			   DietScav=datXY$DietScav, DietFruit=datXY$DietFruit, DietNect=datXY$DietNect, 
			   DietSeed=datXY$DietSeed, DietHerb=datXY$DietHerb, BodyMass=datXY$BodyMass, 
			   Ground=datXY$Ground, Climbing=datXY$Climbing, Volant=datXY$Volant)
			colnames(species.traits)[1] <- "species"
			dimnames(species.traits) <- list("species"=as.character(species.traits[,1]),"traits"=dimnames(species.traits)[[2]])
			distances <- dist(species.traits[,-1], method=Distance.method)
			tree <- hclust(distances, method=Cluster.method)
			xtree <- Xtree(tree)
			i.prime <- ifelse(colSums(xtree$H1)>0, 1, 0)
			FD0 <- sum(i.prime * xtree$h2.prime)
			#  Phylogenetic diversity
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(datXY)))
			names(phymat) <- datXY$species
			PD0 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			
			#----- Spp impactes: Habitat loss
			# Removing species where the threat is major (1) and causing population decline
			hab1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & habitat == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hab1 <- nrow(hab1_spp)
			numspploss_hab1 <- SR0 - nrow(hab1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hab1_spp)))
			names(phymat) <- hab1_spp$species
			PD_hab1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_hab1 < 2){FDrand_hab1 <- NA; PDrand_hab1 <- NA}else
			  if(SR_hab1 == SR0){FDrand_hab1 <- 1; PDrand_hab1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_hab1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_hab1 <- mean(FDrand$FD.new)
				PDrand_hab1 <- mean(PDr)
				}
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			hab2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & habitat == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & habitat == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hab2 <- nrow(hab2_spp)
			numspploss_hab2 <- SR0 - nrow(hab2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hab2_spp)))
			names(phymat) <- hab2_spp$species
			PD_hab2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_hab2 < 2){FDrand_hab2 <- NA; PDrand_hab2 <- NA}else
			  if(SR_hab2 == SR0){FDrand_hab2 <- 1; PDrand_hab2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)		
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_hab2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_hab2 <- mean(FDrand$FD.new)
				PDrand_hab2 <- mean(PDr)
				}
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			#----- Spp impacts: Hunting
			# Removing species where the threat is major (1) and causing population decline
			hunt1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & hunting == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hunt1 <- nrow(hunt1_spp)
			numspploss_hunt1 <- SR0 - nrow(hunt1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hunt1_spp)))
			names(phymat) <- hunt1_spp$species
			PD_hunt1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_hunt1 < 2){FDrand_hunt1 <- NA; PDrand_hunt1 <- NA}else
			  if(SR_hunt1 == SR0){FDrand_hunt1 <- 1; PDrand_hunt1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_hunt1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_hunt1 <- mean(FDrand$FD.new)
				PDrand_hunt1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			hunt2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & hunting == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & hunting == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hunt2 <- nrow(hunt2_spp)
			numspploss_hunt2 <- SR0 - nrow(hunt2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hunt2_spp)))
			names(phymat) <- hunt2_spp$species
			PD_hunt2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_hunt2 < 2){FDrand_hunt2 <- NA; PDrand_hunt2 <- NA}else
			  if(SR_hunt2 == SR0){FDrand_hunt2 <- 1; PDrand_hunt2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_hunt2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_hunt2 <- mean(FDrand$FD.new)
				PDrand_hunt2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			#----- Spp impacts: Human-wildlife conflict
			# Removing species where the threat is major (1) and causing population decline
			con1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & conflict == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_con1 <- nrow(con1_spp)
			numspploss_con1 <- SR0 - nrow(con1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(con1_spp)))
			names(phymat) <- con1_spp$species
			PD_con1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_con1 < 2){FDrand_con1 <- NA; PDrand_con1 <- NA}else
			  if(SR_con1 == SR0){FDrand_con1 <- 1; PDrand_con1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_con1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_con1 <- mean(FDrand$FD.new)
				PDrand_con1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			con2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & conflict == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & conflict == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_con2 <- nrow(con2_spp)
			numspploss_con2 <- SR0 - nrow(con2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(con2_spp)))
			names(phymat) <- con2_spp$species
			PD_con2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_con2 < 2){FDrand_con2 <- NA; PDrand_con2 <- NA}else
			  if(SR_con2 == SR0){FDrand_con2 <- 1; PDrand_con2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_con2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_con2 <- mean(FDrand$FD.new)
				PDrand_con2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Climate change
			# Removing species where the threat is major (1) and causing population decline
			clim1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & climate == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_clim1 <- nrow(clim1_spp)
			numspploss_clim1 <- SR0 - nrow(clim1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(clim1_spp)))
			names(phymat) <- clim1_spp$species
			PD_clim1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_clim1 < 2){FDrand_clim1 <- NA; PDrand_clim1 <- NA}else
			  if(SR_clim1 == SR0){FDrand_clim1 <- 1; PDrand_clim1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_clim1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_clim1 <- mean(FDrand$FD.new)
				PDrand_clim1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			clim2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & climate == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & climate == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_clim2 <- nrow(clim2_spp)
			numspploss_clim2 <- SR0 - nrow(clim2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(clim2_spp)))
			names(phymat) <- clim2_spp$species
			PD_clim2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_clim2 < 2){FDrand_clim2 <- NA; PDrand_clim2 <- NA}else
			  if(SR_clim2 == SR0){FDrand_clim2 <- 1; PDrand_clim2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_clim2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_clim2 <- mean(FDrand$FD.new)
				PDrand_clim2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Non-natives
			# Removing species where the threat is major (1) and causing population decline
			nonnat1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & NonNatives == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_nonnat1 <- nrow(nonnat1_spp)
			numspploss_nonnat1 <- SR0 - nrow(nonnat1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(nonnat1_spp)))
			names(phymat) <- nonnat1_spp$species
			PD_nonnat1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_nonnat1 < 2){FDrand_nonnat1 <- NA; PDrand_nonnat1 <- NA}else
			  if(SR_nonnat1 == SR0){FDrand_nonnat1 <- 1; PDrand_nonnat1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_nonnat1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_nonnat1 <- mean(FDrand$FD.new)
				PDrand_nonnat1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			nonnat2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & NonNatives == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & NonNatives == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_nonnat2 <- nrow(nonnat2_spp)
			numspploss_nonnat2 <- SR0 - nrow(nonnat2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(nonnat2_spp)))
			names(phymat) <- nonnat2_spp$species
			PD_nonnat2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_nonnat2 < 2){FDrand_nonnat2 <- NA; PDrand_nonnat2 <- NA}else
			  if(SR_nonnat2 == SR0){FDrand_nonnat2 <- 1; PDrand_nonnat2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_nonnat2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_nonnat2 <- mean(FDrand$FD.new)
				PDrand_nonnat2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Pollution
			# Removing species where the threat is major (1) and causing population decline
			pol1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & pollution == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_pol1 <- nrow(pol1_spp)
			numspploss_pol1 <- SR0 - nrow(pol1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(pol1_spp)))
			names(phymat) <- pol1_spp$species
			PD_pol1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_pol1 < 2){FDrand_pol1 <- NA; PDrand_pol1 <- NA}else
			  if(SR_pol1 == SR0){FDrand_pol1 <- 1; PDrand_pol1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_pol1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_pol1 <- mean(FDrand$FD.new)
				PDrand_pol1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			pol2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & pollution == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & pollution == 2)) %>%
			  dplyr::select(species)	%>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_pol2 <- nrow(pol2_spp)
			numspploss_pol2 <- SR0 - nrow(pol2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(pol2_spp)))
			names(phymat) <- pol2_spp$species
			PD_pol2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_pol2 < 2){FDrand_pol2 <- NA; PDrand_pol2 <- NA}else
			  if(SR_pol2 == SR0){FDrand_pol2 <- 1; PDrand_pol2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_pol2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_pol2 <- mean(FDrand$FD.new)
				PDrand_pol2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Hybridization
			# Removing species where the threat is major (1) and causing population decline
			hyb1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & hybrid == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hyb1 <- nrow(hyb1_spp)
			numspploss_hyb1 <- SR0 - nrow(hyb1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hyb1_spp)))
			names(phymat) <- hyb1_spp$species
			PD_hyb1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_hyb1 < 2){FDrand_hyb1 <- NA; PDrand_hyb1 <- NA}else
			  if(SR_hyb1 == SR0){FDrand_hyb1 <- 1; PDrand_hyb1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_hyb1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_hyb1 <- mean(FDrand$FD.new)
				PDrand_hyb1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			hyb2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & hybrid == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & hybrid == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()	  
			SR_hyb2 <- nrow(hyb2_spp)
			numspploss_hyb2 <- SR0 - nrow(hyb2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hyb2_spp)))
			names(phymat) <- hyb2_spp$species
			PD_hyb2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_hyb2 < 2){FDrand_hyb2 <- NA; PDrand_hyb2 <- NA}else
			  if(SR_hyb2 == SR0){FDrand_hyb2 <- 1; PDrand_hyb2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_hyb2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_hyb2 <- mean(FDrand$FD.new)
				PDrand_hyb2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Prey depletion
			# Removing species where the threat is major (1) and causing population decline
			prey1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & prey == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_prey1 <- nrow(prey1_spp)
			numspploss_prey1 <- SR0 - nrow(prey1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(prey1_spp)))
			names(phymat) <- prey1_spp$species
			PD_prey1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_prey1 < 2){FDrand_prey1 <- NA; PDrand_prey1 <- NA}else
			  if(SR_prey1 == SR0){FDrand_prey1 <- 1; PDrand_prey1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_prey1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_prey1 <- mean(FDrand$FD.new)
				PDrand_prey1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			prey2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & prey == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & prey == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_prey2 <- nrow(prey2_spp)
			numspploss_prey2 <- SR0 - nrow(prey2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(prey2_spp)))
			names(phymat) <- prey2_spp$species
			PD_prey2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_prey2 < 2){FDrand_prey2 <- NA; PDrand_prey2 <- NA}else
			  if(SR_prey2 == SR0){FDrand_prey2 <- 1; PDrand_prey2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_prey2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_prey2 <- mean(FDrand$FD.new)
				PDrand_prey2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Disease
			# Removing species where the threat is major (1) and causing population decline
			dis1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & disease == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_dis1 <- nrow(dis1_spp)
			numspploss_dis1 <- SR0 - nrow(dis1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(dis1_spp)))
			names(phymat) <- dis1_spp$species
			PD_dis1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_dis1 < 2){FDrand_dis1 <- NA; PDrand_dis1 <- NA}else
			  if(SR_dis1 == SR0){FDrand_dis1 <- 1; PDrand_dis1 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_dis1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_dis1 <- mean(FDrand$FD.new)
				PDrand_dis1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			dis2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & disease == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & disease == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_dis2 <- nrow(dis2_spp)
			numspploss_dis2 <- SR0 - nrow(dis2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(dis2_spp)))
			names(phymat) <- dis2_spp$species
			PD_dis2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_dis2 < 2){FDrand_dis2 <- NA; PDrand_dis2 <- NA}else
			  if(SR_dis2 == SR0){FDrand_dis2 <- 1; PDrand_dis2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_dis2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_dis2 <- mean(FDrand$FD.new)
				PDrand_dis2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			#----- Spp impacts: Inbreeding
			# Removing species where the threat is major (1) and causing population decline
			inb1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & inbreeding == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_inb1 <- nrow(inb1_spp)
			numspploss_inb1 <- SR0 - nrow(inb1_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(inb1_spp)))
			names(phymat) <- inb1_spp$species
			PD_inb1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_inb1 < 2){FDrand_inb1 <- NA; PDrand_inb1 <- NA} else
			  if(SR_inb1 == SR0){FDrand_inb1 <- 1; PDrand_inb1 <- 1} else {	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_inb1)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_inb1 <- mean(FDrand$FD.new)
				PDrand_inb1 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			# Removing species where the threat is major or minor (1 or 2) and causing population decline
			inb2_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & inbreeding == 1)) %>%
			  dplyr::filter(!(trend == "Decreasing" & inbreeding == 2)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_inb2 <- nrow(inb2_spp)
			numspploss_inb2 <- SR0 - nrow(inb2_spp)
			
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(inb2_spp)))
			names(phymat) <- inb2_spp$species
			PD_inb2 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			if(SR_inb2 < 2){FDrand_inb2 <- NA; PDrand_inb2 <- NA}else
			  if(SR_inb2 == SR0){FDrand_inb2 <- 1; PDrand_inb2 <- 1}else{	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_inb2)),]) # random sample of spp
					names(sppr) <- names(all_spp)
					rand_comm <- as.data.frame(cbind(rep(r, nrow(sppr)), sppr))
					names(rand_comm)[1] <- "comm_idx"
					community.composition <- rbind(community.composition, rand_comm)
					# Phyl div in the community given the same amount of spp loss, but random spp
					phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(sppr)))
					names(phymat) <- sppr$species
					PDr[r,1] <- picante::pd(phymat, spptree, include.root=TRUE)$PD    
					}	
				colnames(community.composition)[1] <- "community"
				FDrand <- Getlength(xtree, community.composition)
				FDrand_inb2 <- mean(FDrand$FD.new)
				PDrand_inb2 <- mean(PDr)
			  }
			rm(phymat, FD_rand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
			
			# Any communities completely gone (0 spp)?
			FD_new_0s <- c(nrow(hab1_spp), nrow(hab2_spp),
				nrow(hunt1_spp), nrow(hunt2_spp), 
				nrow(con1_spp), nrow(con2_spp),
				nrow(clim1_spp), nrow(clim2_spp),
				nrow(nonnat1_spp), nrow(nonnat2_spp),
				nrow(pol1_spp), nrow(pol2_spp),
				nrow(hyb1_spp), nrow(hyb2_spp),
				nrow(prey1_spp), nrow(prey2_spp),
				nrow(dis1_spp), nrow(dis2_spp),
				nrow(inb1_spp), nrow(inb2_spp))
			FD_new_0s_vec <- which(FD_new_0s == 0)
			
			#----- New diversity
			# Generate community composition matrix for each scenario
			comm_comp_spp <- rbind(all_spp, 
				hab1_spp, hab2_spp, 
				hunt1_spp, hunt2_spp,
				con1_spp, con2_spp, 
				clim1_spp, clim2_spp, 
				nonnat1_spp, nonnat2_spp,
				pol1_spp, pol2_spp, 
				hyb1_spp, hyb2_spp, 
				prey1_spp, prey2_spp,
				dis1_spp, dis2_spp, 
				inb1_spp, inb2_spp)
			comm_comp_idx <- c(rep(1, nrow(datXY)), 
				rep(2, nrow(hab1_spp)), rep(3, nrow(hab2_spp)), 
				rep(4, nrow(hunt1_spp)), rep(5, nrow(hunt2_spp)),
				rep(6, nrow(con1_spp)), rep(7, nrow(con2_spp)), 
				rep(8, nrow(clim1_spp)), rep(9, nrow(clim2_spp)), 
				rep(10, nrow(nonnat1_spp)), rep(11, nrow(nonnat2_spp)),
				rep(12, nrow(pol1_spp)), rep(13, nrow(pol2_spp)), 
				rep(14, nrow(hyb1_spp)), rep(15, nrow(hyb2_spp)), 
				rep(16, nrow(prey1_spp)), rep(17, nrow(prey2_spp)),
				rep(18, nrow(dis1_spp)), rep(19, nrow(dis2_spp)),
				rep(20, nrow(inb1_spp)), rep(21, nrow(inb2_spp)))
			community.composition <- as.data.frame(cbind(comm_comp_idx, comm_comp_spp))
			colnames(community.composition)[1] <- "community"
			
			FD <- Getlength(xtree, community.composition)
			FD_c <- FD[2]
			FD_r_tmp1 <- t(FD_c)
			FD_r_tmp2 <- FD_r_tmp1[-1]
			
			if(any(FD_new_0s == 0) == TRUE){
				ins_vals <- rep(0, length(FD_new_0s_vec))
				FD_r <- R.utils::insert(FD_r_tmp2, ats = FD_new_0s_vec, values = ins_vals)
				} else {
					FD_r <- FD_r_tmp2
					}
				}	
  
  res_vec <- c(SR0, 
               SR_hab1, SR_hab2, 
			   SR_hunt1, SR_hunt2, 
			   SR_clim1, SR_clim2, 
			   SR_con1, SR_con2, 
               SR_nonnat1, SR_nonnat2, 
			   SR_pol1, SR_pol2, 
			   SR_hyb1, SR_hyb2, 
			   SR_prey1, SR_prey2, 
               SR_dis1, SR_dis2, 
			   SR_inb1, SR_inb2, 
               FD0, 
			   FD_r, #length 20
               PD0, 
               PD_hab1/PD0, PD_hab2/PD0, 
			   PD_hunt1/PD0, PD_hunt2/PD0, 
			   PD_clim1/PD0, PD_clim2/PD0, 
               PD_con1/PD0, PD_con2/PD0, 
			   PD_nonnat1/PD0, PD_nonnat2/PD0, 
			   PD_pol1/PD0, PD_pol2/PD0, 
               PD_hyb1/PD0, PD_hyb2/PD0, 
			   PD_prey1/PD0, PD_prey2/PD0, 
			   PD_dis1/PD0, PD_dis2/PD0, 
               PD_inb1/PD0, PD_inb2/PD0, 
               FDrand_hab1, FDrand_hab2, 
			   FDrand_hunt1, FDrand_hunt2, 
			   FDrand_clim1, FDrand_clim2, 
               FDrand_con1, FDrand_con2,
			   FDrand_nonnat1, FDrand_nonnat2, 
			   FDrand_pol1, FDrand_pol2, 
			   FDrand_hyb1, FDrand_hyb2, 
			   FDrand_prey1, FDrand_prey2, 
			   FDrand_dis1, FDrand_dis2, 
			   FDrand_inb1, FDrand_inb2,	
               PDrand_hab1/PD0, PDrand_hab2/PD0, 
			   PDrand_hunt1/PD0, PDrand_hunt2/PD0, 
               PDrand_clim1/PD0, PDrand_clim2/PD0, 
			   PDrand_con1/PD0, PDrand_con2/PD0, 
               PDrand_nonnat1/PD0, PDrand_nonnat2/PD0, 
			   PDrand_pol1/PD0, PDrand_pol2/PD0, 
               PDrand_hyb1/PD0, PDrand_hyb2/PD0, 
			   PDrand_prey1/PD0, PDrand_prey2/PD0, 
               PDrand_dis1/PD0, PDrand_dis2/PD0, 
			   PDrand_inb1/PD0, PDrand_inb2/PD0)
  
	if(res_vec[1] == 0){
		res_vec[] <- 0
		} else {
			res_vec <- res_vec
			}


return(res_vec)

	
	
rm("datXY", "all_spp",
	"country_name",
	"keep_exc", "spp_exc", "spp_exc_df", 
	"spp_list", 
	"SR0", 
	"hab1_spp", "hab2_spp", "hunt1_spp", "hunt2_spp", "hyb1_spp", "hyb2_spp", "inb1_spp" , "inb2_spp",
	"nonnat1_spp", "nonnat2_spp", "pol1_spp", "pol2_spp", "prey1_spp", "prey2_spp", "con1_spp", "con2_spp",
	"dis1_spp", "dis2_spp", "clim1_spp" , "clim2_spp",
	"numspploss_clim1", "numspploss_clim2", "numspploss_con1", "numspploss_con2", "numspploss_dis1", 
	"numspploss_dis2", "numspploss_hab1", "numspploss_hab2", "numspploss_hunt1",  "numspploss_hunt2", 
	"numspploss_hyb1", "numspploss_hyb2", "numspploss_inb1", "numspploss_inb2", "numspploss_nonnat1" ,
	"numspploss_nonnat2", "numspploss_pol1", "numspploss_pol2", "numspploss_prey1", "numspploss_prey2",	
	"SR_hab1", "SR_hab2", "SR_hunt1", "SR_hunt2", "SR_clim1", "SR_clim2", "SR_con1", "SR_con2", 
	"SR_nonnat1", "SR_nonnat2", "SR_pol1", "SR_pol2", "SR_hyb1", "SR_hyb2", "SR_prey1", "SR_prey2", 
	"SR_dis1", "SR_dis2", "SR_inb1", "SR_inb2" , 
	"i.prime", "tree", "xtree", "distances", 
	"FD0", "FDrand",
	"FD_r", "FD_r_tmp1", "FD_r_tmp2", "FD_c", "FD", "FD_new_0s", "FD_new_0s_vec",
	"comm_comp_idx" , "comm_comp_spp", "community.composition", "species.traits",
	"PD0", 
	"PD_hab1", "PD_hab2", "PD_hunt1", "PD_hunt2", "PD_clim1", "PD_clim2", 
	"PD_con1", "PD_con2", "PD_nonnat1", "PD_nonnat2", "PD_pol1", "PD_pol2", 
	"PD_hyb1", "PD_hyb2", "PD_prey1", "PD_prey2", "PD_dis1", "PD_dis2", 
	"PD_inb1", "PD_inb2", "FD_r_tmp1", "FD_r_tmp1", "FD_r", "FD", "FD_c", "FD_new_0s", "FD_new_0s_vec",
	"comm_comp_idx", "comm_comp_spp", "community.composition",
	"FDrand_hab1", "FDrand_hab2", "FDrand_hunt1", "FDrand_hunt2", "FDrand_clim1", "FDrand_clim2", 
	"FDrand_con1", "FDrand_con2", "FDrand_nonnat1", "FDrand_nonnat2", "FDrand_pol1", 
	"FDrand_pol2", "FDrand_hyb1", "FDrand_hyb2", "FDrand_prey1", "FDrand_prey2", "FDrand_dis1", 
	"FDrand_dis2", "FDrand_inb1", "FDrand_inb2",	
	"PDrand_hab1", "PDrand_hab2", "PDrand_hunt1", "PDrand_hunt2", 
	"PDrand_clim1", "PDrand_clim2", "PDrand_con1", "PDrand_con2", 
	"PDrand_nonnat1", "PDrand_nonnat2", "PDrand_pol1", "PDrand_pol2", 
	"PDrand_hyb1", "PDrand_hyb2", "PDrand_prey1", "PDrand_prey2", 
	"PDrand_dis1", "PDrand_dis2", "PDrand_inb1", "PDrand_inb2",
	"res_vec")
	

	}



stopCluster(cluster); print("Cluster stopped.")
registerDoSEQ()


end <- Sys.time()
elapsed <- end - start
elapsed

out_no_59 <- out
out_no_59[59] <- NA
out_df <- do.call(rbind, out_no_59)

out_with_59 <- out_df


## At this point rerun the loop a single time for i = 59. Replace row 59 of the results datafarme 9out_df)
##  with the output from that loop. The loop must be run manually so that you can replace PD_hunt2 with 0,
##   instead the error that occurs when trying to calculate PD_hunt2 because SR)hunt2 = 0.
## PD_hunt2 <- 0
##  This is the only country that this situation occurs!!!!

out_with_59[59,]  <- res_vec




#----------- ASSEMBLE AND SAVE THE RESULTS MATRIX -------------------------------
#  Save results as dataframe and CSV file
#res_df <- cbind.data.frame(ctrs_area, out_df)
res_df <- cbind.data.frame(ctrs_area, out_with_59)
names(res_df) <- c("Country", "Area_sq_km",
                   "SR0", "SR.habitat1", "SR.habitat2", "SR.hunting1", "SR.hunting2", "SR.climate1", "SR.climate2", "SR.conflict1", 
                   "SR.conflict2", "SR.nonNatives1", "SR.nonNatives2", "SR.pollution1", "SR.pollution2", "SR.hybrid1", "SR.hybrid2", 
                   "SR.prey1", "SR.prey2", "SR.disease1", "SR.disease2", "SR.inbreeding1", "SR.inbreeding2",
                   "FD0", "FD.habitat1", "FD.habitat2", "FD.hunting1", "FD.hunting2", "FD.climate1", "FD.climate2", "FD.conflict1", 
                   "FD.conflict2", "FD.nonNatives1", "FD.nonNatives2", "FD.pollution1", "FD.pollution2", "FD.hybrid1", "FD.hybrid2", 
                   "FD.prey1", "FD.prey2", "FD.disease1", "FD.disease2", "FD.inbreeding1", "FD.inbreeding2",
                   "PD0", "PD.habitat1", "PD.habitat2", "PD.hunting1", "PD.hunting2", "PD.climate1", "PD.climate2", "PD.conflict1", 
                   "PD.conflict2", "PD.nonNatives1", "PD.nonNatives2", "PD.pollution1", "PD.pollution2", "PD.hybrid1", "PD.hybrid2", 
                   "PD.prey1", "PD.prey2", "PD.disease1", "PD.disease2", "PD.inbreeding1", "PD.inbreeding2",
                   "FDrand.habitat1", "FDrand.habitat2", "FDrand.hunting1", "FDrand.hunting2", "FDrand.climate1", "FDrand.climate2", 
                   "FDrand.conflict1", "FDrand.conflict2", "FDrand.nonNatives1", "FDrand.nonNatives2", "FDrand.pollution1", "FDrand.pollution2", 
                   "FDrand.hybrid1", "FDrand.hybrid2", "FDrand.prey1", "FDrand.prey2", "FDrand.disease1", "FDrand.disease2", 
                   "FDrand.inbreeding1", "FDrand.inbreeding2",
                   "PDrand.habitat1", "PDrand.habitat2", "PDrand.hunting1", "PDrand.hunting2", "PDrand.climate1", "PDrand.climate2", 
                   "PDrand.conflict1", "PDrand.conflict2", "PDrand.nonNatives1", "PDrand.nonNatives2", "PDrand.pollution1", "PDrand.pollution2", 
                   "PDrand.hybrid1", "PDrand.hybrid2", "PDrand.prey1", "PDrand.prey2", "PDrand.disease1", "PDrand.disease2", 
                   "PDrand.inbreeding1", "PDrand.inbreeding2")
write.csv(res_df, "FuncPhylo_mammals190625.csv", row.names=F) 



















