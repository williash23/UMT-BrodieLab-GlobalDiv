#setwd("D:/Box Sync/Projects/Projects (active)/Functional diversity/Analysis")
setwd("C:/Users/JB/Box Sync/Projects/Projects (active)/Functional diversity/Analysis")
#setwd("C:/Users/saraw/Documents/FD")


library(dplyr)
library(tidyr)
library(sf)
library(picante)
library(raster)
library(ape)



rm(list = ls())



# ----------------------- USER-DEFINED PARAMETERS ----------------------------------------------------------------
# Load functions
source("FunctionalDiversity190528_funs.R")
st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))

# Desired patial projection
prj <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
	
# Number of randomization iterations to calculate bias in funct & phylogen loss
numrands <- 501  # 1 will be discarded

# Bonferroni-correction
Ntot <- 232 # number of coutries for which we calculate FD and PD (i.e. where SR0>1)
bc <- (0.05/Ntot)/2 # two-tailed; family-wise alpha = 0.05

# Parameters for Petchy's functional diversity calculation
Distance.method <- "euclidean"
Cluster.method <- "average"



#----------------------- SPECIES DATA -------------------------------------------------
dat1 <- data.frame(read.csv("Data/Raw/mammal_threats_traits.csv", header=T))
dat1$ExceptionsToDecline[dat1$ExceptionsToDecline == "na"] <- NA
spptree <- read.tree("Data/Raw/mammaltree.tree")
dat1$species <- as.character(dat1$species)




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

# Run functions to prepare range polygons
# removes mammals who are listed as extinct
mam <- prep_range_polys(fn = "Data/Raw/TERRESTRIAL_MAMMALS.shp", ctrs = ctrs)
#  This data set was obtained from the IUCN website: 
#   https://www.iucnredlist.org/resources/spatial-data-download
#  Downloaded on 2017-06-15
#  Processing function found in script: FunctionalDiversity190528_funs.R


#----------------------- SPECIES - COUNTRIES INTERSECTION ---------------------------------
polys_ctrs_int <- st_intersects(mam, ctrs, sparse = FALSE) 
# Dimensions: nrow = number of spp (12112 for mammals), ncol = number of countries (254)

spp_int_ls <- list()
for(i in 1:nrow(ctrs_df)){
	polys_ctrs_int_vec <- which(polys_ctrs_int[,i] == TRUE)
	spp_list_ctr <- mam[polys_ctrs_int_vec, 1]
	st_geometry(spp_list_ctr) <- NULL
	spp_int_ls[[i]] <- spp_list_ctr	}

	

#----------------------- MANUALLY ADD CERTAIN SPECIES ---------------------------------
# Species present in IUCN range descriptions but not IUCN range maps
# Maldives
add <- data.frame(species = "Myotis_blythii")
spp_int_ls[[59]] <- rbind(spp_int_ls[[59]], add)
# Turkey
add <- data.frame(species = "Miniopterus_schreibersii")
spp_int_ls[[16]] <- rbind(spp_int_ls[[16]], add)
# France
add <- data.frame(species = "Speothos_venaticus")
spp_int_ls[[34]] <- rbind(spp_int_ls[[34]], add)
# Sri Lanka
add <- data.frame(species = "Myotis_blythii")
spp_int_ls[[63]] <- rbind(spp_int_ls[[63]], add)
# Bangladesh
add <- data.frame(species = "Myotis_blythii")
spp_int_ls[[80]] <- rbind(spp_int_ls[[80]], add)
# Bhutan
add <- data.frame(species = "Myotis_blythii")
spp_int_ls[[81]] <- rbind(spp_int_ls[[81]], add)
# Israel
add <- data.frame(species = "Gazella_dorcas")
spp_int_ls[[224]] <- rbind(spp_int_ls[[224]], add)
# Croatia
add <- data.frame(species = "Miniopterus_schreibersii")
spp_int_ls[[248]] <- rbind(spp_int_ls[[248]], add)
# Bulgaria
add <- data.frame(species = "Miniopterus_schreibersii")
spp_int_ls[[254]] <- rbind(spp_int_ls[[254]], add)



#----------------------- EXCEPTIONS TO DECLINE ---------------------------------
exc_tmp <- dat1 %>%
	dplyr::select(species, ExceptionsToDecline)
exc_tmp$ExceptionsToDecline <- as.character(exc_tmp$ExceptionsToDecline)

#  Do first one to start data frame
spp <- exc_tmp$species[1]
tmp <- exc_tmp$ExceptionsToDecline[1]
tmp_df <- as.data.frame(tmp) %>%
	tidyr::separate(tmp, paste("Country", 1:22, sep="_"), sep = ", ") 
exc_df <- cbind(spp, tmp_df)

for(i in 1:nrow(exc_tmp)){
	spp <- exc_tmp$species[i]
	tmp <- exc_tmp$ExceptionsToDecline[i]
	tmp_df <- as.data.frame(tmp) %>%
		tidyr::separate(tmp, paste("Country", 1:22, sep="_"), sep = ", ") 
	new_row <- cbind(spp, tmp_df)
	exc_df <- rbind(exc_df, new_row)	}
	
exc_df <- exc_df[rowSums(is.na(exc_df[,2:23])) != 22,]
exc_df <- exc_df %>%
	mutate_if(is.factor, as.character)

dat1 <- subset(dat1, select = -c(ExceptionsToDecline))

	

#----------------------- SAVE DATA FOR LATER USE / LOAD PREVIOUSLY PREPPED DATA ----------------
#save(spp_int_ls, file="Data/spp_int_ls.Rdata")
#save(dat1, file="Data/dat1.Rdata")
#save(ctrs, file="Data/ctrs.Rdata")
#save(ctrs_df, file="Data/ctrs_df.Rdata")

#save.image("working")
#load("working")



#----------------------- ASSESS THREATS PER COUNTRY ------------------------------
### 'FOR' LOOP, RUNNING THROUGH EACH COUNTRY ONE BY ONE

# Number of countries to analyze
numcountries <- length(spp_int_ls)

out_mat <- matrix(NA, nrow = numcountries, ncol = 170)
start <- Sys.time()
start

for(i in 1:numcountries){
tryCatch({


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

	t0.01 <- mean(datXY$DietInv)
	t0.02 <- mean(datXY$DietVert)
	t0.03 <- mean(datXY$DietFish)
	t0.04 <- mean(datXY$DietScav)
	t0.05 <- mean(datXY$DietFruit)
	t0.06 <- mean(datXY$DietNect)
	t0.07 <- mean(datXY$DietSeed)
	t0.08 <- mean(datXY$DietHerb)
	t0.09 <- mean(datXY$BodyMass)
	t0.10 <- mean(datXY$Ground)
	t0.11 <- mean(datXY$Climbing)
	t0.12 <- mean(datXY$Volant)


	#----- Current diversity
  
	if(SR0 < 2){
		FD0 <- NA; FD_r <- rep(NA, 11); 
		FDrand_allthreats <- NA; FDrand_hab1 <- NA; FDrand_hunt1 <- NA; FDrand_clim1 <- NA; FDrand_con1 <- NA; 
		FDrand_nonnat1 <- NA; FDrand_pol1 <- NA; FDrand_hyb1 <- NA; FDrand_prey1 <- NA; FDrand_dis1 <- NA; FDrand_inb1 <- NA; 
		FDrandLO_allthreats <- NA; FDrandLO_hab1 <- NA; FDrandLO_hunt1 <- NA; FDrandLO_clim1 <- NA; FDrandLO_con1 <- NA; 
		FDrandLO_nonnat1 <- NA; FDrandLO_pol1 <- NA; FDrandLO_hyb1 <- NA; FDrandLO_prey1 <- NA; FDrandLO_dis1 <- NA; FDrandLO_inb1 <- NA; 
		FDrandHI_allthreats <- NA; FDrandHI_hab1 <- NA; FDrandHI_hunt1 <- NA; FDrandHI_clim1 <- NA; FDrandHI_con1 <- NA; 
		FDrandHI_nonnat1 <- NA; FDrandHI_pol1 <- NA; FDrandHI_hyb1 <- NA; FDrandHI_prey1 <- NA; FDrandHI_dis1 <- NA; FDrandHI_inb1 <- NA; 
		FDrandLObonferroni_allthreats <- NA; FDrandLObonferroni_hab1 <- NA; FDrandLObonferroni_hunt1 <- NA; FDrandLObonferroni_clim1 <- NA; 
		FDrandLObonferroni_con1 <- NA; FDrandLObonferroni_nonnat1 <- NA; FDrandLObonferroni_pol1 <- NA; FDrandLObonferroni_hyb1 <- NA; 
		FDrandLObonferroni_prey1 <- NA; FDrandLObonferroni_dis1 <- NA; FDrandLObonferroni_inb1 <- NA; 
		FDrandHIbonferroni_allthreats <- NA; FDrandHIbonferroni_hab1 <- NA; FDrandHIbonferroni_hunt1 <- NA; FDrandHIbonferroni_clim1 <- NA; 
		FDrandHIbonferroni_con1 <- NA; FDrandHIbonferroni_nonnat1 <- NA; FDrandHIbonferroni_pol1 <- NA; FDrandHIbonferroni_hyb1 <- NA; 
		FDrandHIbonferroni_prey1 <- NA; FDrandHIbonferroni_dis1 <- NA; FDrandHIbonferroni_inb1 <- NA; 
		PD0 <- NA; PD_allthreats <- NA; PD_hab1 <- NA; PD_hunt1 <- NA; PD_clim1 <- NA; PD_con1 <- NA; PD_nonnat1 <- NA; 
		PD_pol1 <- NA; PD_hyb1 <- NA; PD_prey1 <- NA; PD_dis1 <- NA; PD_inb1 <- NA; 
		PDrand_allthreats <- NA; PDrand_hab1 <- NA; PDrand_hunt1 <- NA; PDrand_clim1 <- NA; PDrand_con1 <- NA; 
		PDrand_nonnat1 <- NA; PDrand_pol1 <- NA; PDrand_hyb1 <- NA; PDrand_prey1 <- NA; PDrand_dis1 <- NA; PDrand_inb1 <- NA; 
		PDrandLO_allthreats <- NA; PDrandLO_hab1 <- NA; PDrandLO_hunt1 <- NA; PDrandLO_clim1 <- NA; PDrandLO_con1 <- NA; 
		PDrandLO_nonnat1 <- NA; PDrandLO_pol1 <- NA; PDrandLO_hyb1 <- NA; PDrandLO_prey1 <- NA; PDrandLO_dis1 <- NA; PDrandLO_inb1 <- NA; 
		PDrandHI_allthreats <- NA; PDrandHI_hab1 <- NA; PDrandHI_hunt1 <- NA; PDrandHI_clim1 <- NA; PDrandHI_con1 <- NA; 
		PDrandHI_nonnat1 <- NA; PDrandHI_pol1 <- NA; PDrandHI_hyb1 <- NA; PDrandHI_prey1 <- NA; PDrandHI_dis1 <- NA; PDrandHI_inb1 <- NA; 
		PDrandLObonferroni_allthreats <- NA; PDrandLObonferroni_hab1 <- NA; PDrandLObonferroni_hunt1 <- NA; PDrandLObonferroni_clim1 <- NA; 
		PDrandLObonferroni_con1 <- NA; PDrandLObonferroni_nonnat1 <- NA; PDrandLObonferroni_pol1 <- NA; PDrandLObonferroni_hyb1 <- NA; 
		PDrandLObonferroni_prey1 <- NA; PDrandLObonferroni_dis1 <- NA; PDrandLObonferroni_inb1 <- NA; 
		PDrandHIbonferroni_allthreats <- NA; PDrandHIbonferroni_hab1 <- NA; PDrandHIbonferroni_hunt1 <- NA; PDrandHIbonferroni_clim1 <- NA; 
		PDrandHIbonferroni_con1 <- NA; PDrandHIbonferroni_nonnat1 <- NA; PDrandHIbonferroni_pol1 <- NA; PDrandHIbonferroni_hyb1 <- NA; 
		PDrandHIbonferroni_prey1 <- NA; PDrandHIbonferroni_dis1 <- NA; PDrandHIbonferroni_inb1 <- NA; 

		allthreats_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing")) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_allthreats <- nrow(allthreats_spp)
	    
		hab1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & habitat == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hab1 <- nrow(hab1_spp)
	    
		hunt1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & hunting == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hunt1 <- nrow(hunt1_spp)
		
		clim1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & climate == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_clim1 <- nrow(clim1_spp)

		con1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & conflict == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_con1 <- nrow(con1_spp)
		
  		nonnat1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & NonNatives == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_nonnat1 <- nrow(nonnat1_spp)
		
		pol1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & pollution == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_pol1 <- nrow(pol1_spp)
	  
		hyb1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & hybrid == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_hyb1 <- nrow(hyb1_spp)
	  
		prey1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & prey == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_prey1 <- nrow(prey1_spp)
		
		dis1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & disease == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_dis1 <- nrow(dis1_spp)
	  
		inb1_spp <- datXY %>%
			dplyr::filter(!(trend == "Decreasing" & inbreeding == 1)) %>%
			dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
		SR_inb1 <- nrow(inb1_spp)
	  
		hab1_tmp <- NA; hunt1_tmp <- NA
		thab.01 <- NA; thab.02 <- NA; thab.03 <- NA; thab.04 <- NA; thab.05 <- NA; thab.06 <- NA 
		thab.07 <- NA; thab.08 <- NA; thab.09 <- NA; thab.10 <- NA; thab.11 <- NA; thab.12 <- NA 
		thunt.01 <- NA; thunt.02 <- NA; thunt.03 <- NA; thunt.04 <- NA; thunt.05 <- NA; thunt.06 <- NA 
		thunt.07 <- NA; thunt.08 <- NA; thunt.09 <- NA; thunt.10 <- NA; thunt.11 <- NA; thunt.12 <- NA 


		} else {

			
			#----- Standing diversity
			# Functional diversity
			# Trait dendrogram for all species at the grid point
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
			# Phylogenetic diversity
			phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(datXY)))
			names(phymat) <- datXY$species
			PD0 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
			
			

			#----- Spp impacts: All threats
			# Removing species undergoing population decline
			allthreats_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing")) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_allthreats <- nrow(allthreats_spp)
			numspploss_allthreats <- SR0 - nrow(allthreats_spp)
			
			if(SR_allthreats == 0){PD_allthreats <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(allthreats_spp)))
				names(phymat) <- allthreats_spp$species
				PD_allthreats <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_allthreats < 2){FDrand_allthreats <- NA; FDrandLO_allthreats <- NA; FDrandHI_allthreats <- NA; 
			  FDrandLObonferroni_allthreats <- NA; FDrandHIbonferroni_allthreats <- NA; PDrand_allthreats <- NA;
			  PDrandLO_allthreats <- NA; PDrandHI_allthreats <- NA; 
			  PDrandLObonferroni_allthreats <- NA; PDrandHIbonferroni_allthreats <- NA} else
			 
			 if(SR_allthreats == SR0){FDrand_allthreats <- 1; FDrandLO_allthreats <- 1; FDrandHI_allthreats <- 1; 
			  FDrandLObonferroni_allthreats <- 1; FDrandHIbonferroni_allthreats <- 1; PDrand_allthreats <- 1;
			  PDrandLO_allthreats <- 1; PDrandHI_allthreats <- 1; 
			  PDrandLObonferroni_allthreats <- 1; PDrandHIbonferroni_allthreats <- 1} else {	
				rep_comm_1 <- rep(1, nrow(datXY))
				community.composition <- as.data.frame(cbind(rep_comm_1, all_spp))
				names(community.composition)[1] <- "comm_idx"
				PDr <- matrix(, nrow=numrands, 1)	
				for(r in 1:numrands){	
					sppr <- data.frame(all_spp[sample(nrow(all_spp), (SR0-numspploss_allthreats)),]) # random sample of spp
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
				FDrand <- FDrand[-1,]
				FDrand_allthreats <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_allthreats <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_allthreats <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_allthreats <- mean(PDr, na.rm=T); 
				PDrandLO_allthreats <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_allthreats <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_allthreats <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_allthreats <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_allthreats <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_allthreats <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			

		
			#----- Spp impacts: Habitat loss
			# Removing species where the threat is major (1) and causing population decline
			hab1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & habitat == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hab1 <- nrow(hab1_spp)
			numspploss_hab1 <- SR0 - nrow(hab1_spp)
			
			hab1_tmp <- dplyr::left_join(hab1_spp, dat1, by="species") 
			hab1_tmp <- hab1_tmp[complete.cases(hab1_tmp),]
			thab.01 <- mean(hab1_tmp$DietInv)
			thab.02 <- mean(hab1_tmp$DietVert)
			thab.03 <- mean(hab1_tmp$DietFish)
			thab.04 <- mean(hab1_tmp$DietScav)
			thab.05 <- mean(hab1_tmp$DietFruit)
			thab.06 <- mean(hab1_tmp$DietNect)
			thab.07 <- mean(hab1_tmp$DietSeed)
			thab.08 <- mean(hab1_tmp$DietHerb)
			thab.09 <- mean(hab1_tmp$BodyMass)
			thab.10 <- mean(hab1_tmp$Ground)
			thab.11 <- mean(hab1_tmp$Climbing)
			thab.12 <- mean(hab1_tmp$Volant)

			if(SR_hab1 == 0){PD_hab1 <- 0}else{
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hab1_spp)))
				names(phymat) <- hab1_spp$species
				PD_hab1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
				
			if(SR_hab1 < 2){FDrand_hab1 <- NA; FDrandLO_hab1 <- NA; FDrandHI_hab1 <- NA; 
			  FDrandLObonferroni_hab1 <- NA; FDrandHIbonferroni_hab1 <- NA; PDrand_hab1 <- NA;
			  PDrandLO_hab1 <- NA; PDrandHI_hab1 <- NA; 
			  PDrandLObonferroni_hab1 <- NA; PDrandHIbonferroni_hab1 <- NA} else
			 
			 if(SR_hab1 == SR0){FDrand_hab1 <- 1; FDrandLO_hab1 <- 1; FDrandHI_hab1 <- 1; 
			  FDrandLObonferroni_hab1 <- 1; FDrandHIbonferroni_hab1 <- 1; PDrand_hab1 <- 1;
			  PDrandLO_hab1 <- 1; PDrandHI_hab1 <- 1; 
			  PDrandLObonferroni_hab1 <- 1; PDrandHIbonferroni_hab1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_hab1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_hab1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_hab1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_hab1 <- mean(PDr, na.rm=T); 
				PDrandLO_hab1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_hab1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_hab1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_hab1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_hab1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_hab1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
				}
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
	
		
		
			#----- Spp impacts: Hunting
			# Removing species where the threat is major (1) and causing population decline
			hunt1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & hunting == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hunt1 <- nrow(hunt1_spp)
			numspploss_hunt1 <- SR0 - nrow(hunt1_spp)
			
			hunt1_tmp <- dplyr::left_join(hunt1_spp, dat1, by="species") 
			hunt1_tmp <- hunt1_tmp[complete.cases(hunt1_tmp),]
			thunt.01 <- mean(hunt1_tmp$DietInv)
			thunt.02 <- mean(hunt1_tmp$DietVert)
			thunt.03 <- mean(hunt1_tmp$DietFish)
			thunt.04 <- mean(hunt1_tmp$DietScav)
			thunt.05 <- mean(hunt1_tmp$DietFruit)
			thunt.06 <- mean(hunt1_tmp$DietNect)
			thunt.07 <- mean(hunt1_tmp$DietSeed)
			thunt.08 <- mean(hunt1_tmp$DietHerb)
			thunt.09 <- mean(hunt1_tmp$BodyMass)
			thunt.10 <- mean(hunt1_tmp$Ground)
			thunt.11 <- mean(hunt1_tmp$Climbing)
			thunt.12 <- mean(hunt1_tmp$Volant)

			if(SR_hunt1 == 0){PD_hunt1 <- 0}else{
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hunt1_spp)))
				names(phymat) <- hunt1_spp$species
				PD_hunt1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_hunt1 < 2){FDrand_hunt1 <- NA; FDrandLO_hunt1 <- NA; FDrandHI_hunt1 <- NA; 
			  FDrandLObonferroni_hunt1 <- NA; FDrandHIbonferroni_hunt1 <- NA; PDrand_hunt1 <- NA;
			  PDrandLO_hunt1 <- NA; PDrandHI_hunt1 <- NA; 
			  PDrandLObonferroni_hunt1 <- NA; PDrandHIbonferroni_hunt1 <- NA} else
			 
			 if(SR_hunt1 == SR0){FDrand_hunt1 <- 1; FDrandLO_hunt1 <- 1; FDrandHI_hunt1 <- 1; 
			  FDrandLObonferroni_hunt1 <- 1; FDrandHIbonferroni_hunt1 <- 1; PDrand_hunt1 <- 1;
			  PDrandLO_hunt1 <- 1; PDrandHI_hunt1 <- 1; 
			  PDrandLObonferroni_hunt1 <- 1; PDrandHIbonferroni_hunt1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_hunt1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_hunt1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_hunt1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_hunt1 <- mean(PDr, na.rm=T); 
				PDrandLO_hunt1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_hunt1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_hunt1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_hunt1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_hunt1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_hunt1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
	
		
	
			#----- Spp impacts: Climate change
			# Removing species where the threat is major (1) and causing population decline
			clim1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & climate == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_clim1 <- nrow(clim1_spp)
			numspploss_clim1 <- SR0 - nrow(clim1_spp)
			
			if(SR_clim1 == 0){PD_clim1 <- 0}else{
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(clim1_spp)))
				names(phymat) <- clim1_spp$species
				PD_clim1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_clim1 < 2){FDrand_clim1 <- NA; FDrandLO_clim1 <- NA; FDrandHI_clim1 <- NA; 
			  FDrandLObonferroni_clim1 <- NA; FDrandHIbonferroni_clim1 <- NA; PDrand_clim1 <- NA;
			  PDrandLO_clim1 <- NA; PDrandHI_clim1 <- NA; 
			  PDrandLObonferroni_clim1 <- NA; PDrandHIbonferroni_clim1 <- NA} else
			 
			 if(SR_clim1 == SR0){FDrand_clim1 <- 1; FDrandLO_clim1 <- 1; FDrandHI_clim1 <- 1; 
			  FDrandLObonferroni_clim1 <- 1; FDrandHIbonferroni_clim1 <- 1; PDrand_clim1 <- 1;
			  PDrandLO_clim1 <- 1; PDrandHI_clim1 <- 1; 
			  PDrandLObonferroni_clim1 <- 1; PDrandHIbonferroni_clim1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_clim1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_clim1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_clim1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_clim1 <- mean(PDr, na.rm=T); 
				PDrandLO_clim1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_clim1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_clim1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_clim1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_clim1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_clim1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			

		
			#----- Spp impacts: Human-wildlife conflict
			# Removing species where the threat is major (1) and causing population decline
			con1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & conflict == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_con1 <- nrow(con1_spp)
			numspploss_con1 <- SR0 - nrow(con1_spp)
			
			if(SR_con1 == 0){PD_con1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(con1_spp)))
				names(phymat) <- con1_spp$species
				PD_con1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_con1 < 2){FDrand_con1 <- NA; FDrandLO_con1 <- NA; FDrandHI_con1 <- NA; 
			  FDrandLObonferroni_con1 <- NA; FDrandHIbonferroni_con1 <- NA; PDrand_con1 <- NA;
			  PDrandLO_con1 <- NA; PDrandHI_con1 <- NA; 
			  PDrandLObonferroni_con1 <- NA; PDrandHIbonferroni_con1 <- NA} else
			 
			 if(SR_con1 == SR0){FDrand_con1 <- 1; FDrandLO_con1 <- 1; FDrandHI_con1 <- 1; 
			  FDrandLObonferroni_con1 <- 1; FDrandHIbonferroni_con1 <- 1; PDrand_con1 <- 1;
			  PDrandLO_con1 <- 1; PDrandHI_con1 <- 1; 
			  PDrandLObonferroni_con1 <- 1; PDrandHIbonferroni_con1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_con1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_con1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_con1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_con1 <- mean(PDr, na.rm=T); 
				PDrandLO_con1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_con1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_con1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_con1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_con1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_con1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			

	
			#----- Spp impacts: Non-natives
			# Removing species where the threat is major (1) and causing population decline
			nonnat1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & NonNatives == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_nonnat1 <- nrow(nonnat1_spp)
			numspploss_nonnat1 <- SR0 - nrow(nonnat1_spp)
			
			if(SR_nonnat1 == 0){PD_nonnat1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(nonnat1_spp)))
				names(phymat) <- nonnat1_spp$species
				PD_nonnat1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_nonnat1 < 2){FDrand_nonnat1 <- NA; FDrandLO_nonnat1 <- NA; FDrandHI_nonnat1 <- NA; 
			  FDrandLObonferroni_nonnat1 <- NA; FDrandHIbonferroni_nonnat1 <- NA; PDrand_nonnat1 <- NA;
			  PDrandLO_nonnat1 <- NA; PDrandHI_nonnat1 <- NA; 
			  PDrandLObonferroni_nonnat1 <- NA; PDrandHIbonferroni_nonnat1 <- NA} else
			 
			 if(SR_nonnat1 == SR0){FDrand_nonnat1 <- 1; FDrandLO_nonnat1 <- 1; FDrandHI_nonnat1 <- 1; 
			  FDrandLObonferroni_nonnat1 <- 1; FDrandHIbonferroni_nonnat1 <- 1; PDrand_nonnat1 <- 1;
			  PDrandLO_nonnat1 <- 1; PDrandHI_nonnat1 <- 1; 
			  PDrandLObonferroni_nonnat1 <- 1; PDrandHIbonferroni_nonnat1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_nonnat1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_nonnat1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_nonnat1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_nonnat1 <- mean(PDr, na.rm=T); 
				PDrandLO_nonnat1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_nonnat1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_nonnat1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_nonnat1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_nonnat1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_nonnat1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			

		
			#----- Spp impacts: Pollution
			# Removing species where the threat is major (1) and causing population decline
			pol1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & pollution == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_pol1 <- nrow(pol1_spp)
			numspploss_pol1 <- SR0 - nrow(pol1_spp)
			
			if(SR_pol1 == 0){PD_pol1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(pol1_spp)))
				names(phymat) <- pol1_spp$species
				PD_pol1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_pol1 < 2){FDrand_pol1 <- NA; FDrandLO_pol1 <- NA; FDrandHI_pol1 <- NA; 
			  FDrandLObonferroni_pol1 <- NA; FDrandHIbonferroni_pol1 <- NA; PDrand_pol1 <- NA;
			  PDrandLO_pol1 <- NA; PDrandHI_pol1 <- NA; 
			  PDrandLObonferroni_pol1 <- NA; PDrandHIbonferroni_pol1 <- NA} else
			 
			 if(SR_pol1 == SR0){FDrand_pol1 <- 1; FDrandLO_pol1 <- 1; FDrandHI_pol1 <- 1; 
			  FDrandLObonferroni_pol1 <- 1; FDrandHIbonferroni_pol1 <- 1; PDrand_pol1 <- 1;
			  PDrandLO_pol1 <- 1; PDrandHI_pol1 <- 1; 
			  PDrandLObonferroni_pol1 <- 1; PDrandHIbonferroni_pol1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_pol1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_pol1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_pol1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_pol1 <- mean(PDr, na.rm=T); 
				PDrandLO_pol1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_pol1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_pol1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_pol1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_pol1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_pol1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
		
			
			#----- Spp impacts: Hybridization
			# Removing species where the threat is major (1) and causing population decline
			hyb1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & hybrid == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_hyb1 <- nrow(hyb1_spp)
			numspploss_hyb1 <- SR0 - nrow(hyb1_spp)
			
			if(SR_hyb1 == 0){PD_hyb1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(hyb1_spp)))
				names(phymat) <- hyb1_spp$species
				PD_hyb1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_hyb1 < 2){FDrand_hyb1 <- NA; FDrandLO_hyb1 <- NA; FDrandHI_hyb1 <- NA; 
			  FDrandLObonferroni_hyb1 <- NA; FDrandHIbonferroni_hyb1 <- NA; PDrand_hyb1 <- NA;
			  PDrandLO_hyb1 <- NA; PDrandHI_hyb1 <- NA; 
			  PDrandLObonferroni_hyb1 <- NA; PDrandHIbonferroni_hyb1 <- NA} else
			 
			 if(SR_hyb1 == SR0){FDrand_hyb1 <- 1; FDrandLO_hyb1 <- 1; FDrandHI_hyb1 <- 1; 
			  FDrandLObonferroni_hyb1 <- 1; FDrandHIbonferroni_hyb1 <- 1; PDrand_hyb1 <- 1;
			  PDrandLO_hyb1 <- 1; PDrandHI_hyb1 <- 1; 
			  PDrandLObonferroni_hyb1 <- 1; PDrandHIbonferroni_hyb1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_hyb1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_hyb1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_hyb1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_hyb1 <- mean(PDr, na.rm=T); 
				PDrandLO_hyb1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_hyb1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_hyb1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_hyb1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_hyb1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_hyb1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
	
			
			#----- Spp impacts: Prey depletion
			# Removing species where the threat is major (1) and causing population decline
			prey1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & prey == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_prey1 <- nrow(prey1_spp)
			numspploss_prey1 <- SR0 - nrow(prey1_spp)
			
			if(SR_prey1 == 0){PD_prey1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(prey1_spp)))
				names(phymat) <- prey1_spp$species
				PD_prey1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_prey1 < 2){FDrand_prey1 <- NA; FDrandLO_prey1 <- NA; FDrandHI_prey1 <- NA; 
			  FDrandLObonferroni_prey1 <- NA; FDrandHIbonferroni_prey1 <- NA; PDrand_prey1 <- NA;
			  PDrandLO_prey1 <- NA; PDrandHI_prey1 <- NA; 
			  PDrandLObonferroni_prey1 <- NA; PDrandHIbonferroni_prey1 <- NA} else
			 
			 if(SR_prey1 == SR0){FDrand_prey1 <- 1; FDrandLO_prey1 <- 1; FDrandHI_prey1 <- 1; 
			  FDrandLObonferroni_prey1 <- 1; FDrandHIbonferroni_prey1 <- 1; PDrand_prey1 <- 1;
			  PDrandLO_prey1 <- 1; PDrandHI_prey1 <- 1; 
			  PDrandLObonferroni_prey1 <- 1; PDrandHIbonferroni_prey1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_prey1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_prey1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_prey1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_prey1 <- mean(PDr, na.rm=T); 
				PDrandLO_prey1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_prey1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_prey1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_prey1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_prey1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_prey1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			
	
			
			#----- Spp impacts: Disease
			# Removing species where the threat is major (1) and causing population decline
			dis1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & disease == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_dis1 <- nrow(dis1_spp)
			numspploss_dis1 <- SR0 - nrow(dis1_spp)
			
			if(SR_dis1 == 0){PD_dis1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(dis1_spp)))
				names(phymat) <- dis1_spp$species
				PD_dis1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_dis1 < 2){FDrand_dis1 <- NA; FDrandLO_dis1 <- NA; FDrandHI_dis1 <- NA; 
			  FDrandLObonferroni_dis1 <- NA; FDrandHIbonferroni_dis1 <- NA; PDrand_dis1 <- NA;
			  PDrandLO_dis1 <- NA; PDrandHI_dis1 <- NA; 
			  PDrandLObonferroni_dis1 <- NA; PDrandHIbonferroni_dis1 <- NA} else
			 
			 if(SR_dis1 == SR0){FDrand_dis1 <- 1; FDrandLO_dis1 <- 1; FDrandHI_dis1 <- 1; 
			  FDrandLObonferroni_dis1 <- 1; FDrandHIbonferroni_dis1 <- 1; PDrand_dis1 <- 1;
			  PDrandLO_dis1 <- 1; PDrandHI_dis1 <- 1; 
			  PDrandLObonferroni_dis1 <- 1; PDrandHIbonferroni_dis1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_dis1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_dis1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_dis1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_dis1 <- mean(PDr, na.rm=T); 
				PDrandLO_dis1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_dis1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_dis1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_dis1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_dis1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_dis1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
			

	
			#----- Spp impacts: Inbreeding
			# Removing species where the threat is major (1) and causing population decline
			inb1_spp <- datXY %>%
			  dplyr::filter(!(trend == "Decreasing" & inbreeding == 1)) %>%
			  dplyr::select(species) %>%
			  rbind(spp_exc_df) %>% 
			  distinct()
			SR_inb1 <- nrow(inb1_spp)
			numspploss_inb1 <- SR0 - nrow(inb1_spp)
			
			if(SR_inb1 == 0){PD_inb1 <- 0} else {
				phymat <- data.frame(matrix(1, nrow=1, ncol=nrow(inb1_spp)))
				names(phymat) <- inb1_spp$species
				PD_inb1 <- picante::pd(phymat, spptree, include.root=TRUE)$PD
				}
			
			if(SR_inb1 < 2){FDrand_inb1 <- NA; FDrandLO_inb1 <- NA; FDrandHI_inb1 <- NA; 
			  FDrandLObonferroni_inb1 <- NA; FDrandHIbonferroni_inb1 <- NA; PDrand_inb1 <- NA;
			  PDrandLO_inb1 <- NA; PDrandHI_inb1 <- NA; 
			  PDrandLObonferroni_inb1 <- NA; PDrandHIbonferroni_inb1 <- NA} else
			 
			 if(SR_inb1 == SR0){FDrand_inb1 <- 1; FDrandLO_inb1 <- 1; FDrandHI_inb1 <- 1; 
			  FDrandLObonferroni_inb1 <- 1; FDrandHIbonferroni_inb1 <- 1; PDrand_inb1 <- 1;
			  PDrandLO_inb1 <- 1; PDrandHI_inb1 <- 1; 
			  PDrandLObonferroni_inb1 <- 1; PDrandHIbonferroni_inb1 <- 1} else {	
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
				FDrand <- FDrand[-1,]
				FDrand_inb1 <- mean(FDrand$FD.new, na.rm=T) 
				FDrandLO_inb1 <- quantile(FDrand$FD.new, probs=c(0.025), na.rm=T)
				FDrandHI_inb1 <- quantile(FDrand$FD.new, probs=c(0.975), na.rm=T)
				PDrand_inb1 <- mean(PDr, na.rm=T); 
				PDrandLO_inb1 <- quantile(PDr, probs=c(0.025), na.rm=T) 
				PDrandHI_inb1 <- quantile(PDr, probs=c(0.975), na.rm=T)
				FDrandLObonferroni_inb1 <- quantile(FDrand$FD.new, probs=c(bc), na.rm=T)
				FDrandHIbonferroni_inb1 <- quantile(FDrand$FD.new, probs=c(1-bc), na.rm=T)
				PDrandLObonferroni_inb1 <- quantile(PDr, probs=c(bc), na.rm=T) 
				PDrandHIbonferroni_inb1 <- quantile(PDr, probs=c(1-bc), na.rm=T) 
			  }
			rm(phymat, FDrand, rep_comm_1, PDr, community.composition, sppr, rand_comm)
		


			#----- Any communities completely gone (0 spp)?
			FD_new_0s <- c(nrow(allthreats_spp), nrow(hab1_spp),
				nrow(hunt1_spp),  
				nrow(con1_spp), 
				nrow(clim1_spp), 
				nrow(nonnat1_spp), 
				nrow(pol1_spp), 
				nrow(hyb1_spp), 
				nrow(prey1_spp), 
				nrow(dis1_spp), 
				nrow(inb1_spp))
			FD_new_0s_vec <- which(FD_new_0s == 0)
			


			#----- New functional diversity
			# Generate community composition matrix for each scenario
			comm_comp_spp <- rbind(all_spp, allthreats_spp, hab1_spp, hunt1_spp, 
				clim1_spp, con1_spp, nonnat1_spp, pol1_spp, hyb1_spp, prey1_spp, 
				dis1_spp, inb1_spp)
			comm_comp_idx <- c(rep(1, nrow(datXY)), 
				rep(2, nrow(allthreats_spp)), 
				rep(3, nrow(hab1_spp)), 
				rep(4, nrow(hunt1_spp)), 
				rep(5, nrow(clim1_spp)),  
				rep(6, nrow(con1_spp)),  
				rep(7, nrow(nonnat1_spp)), 
				rep(8, nrow(pol1_spp)), 
				rep(9, nrow(hyb1_spp)),  
				rep(10, nrow(prey1_spp)), 
				rep(11, nrow(dis1_spp)), 
				rep(12, nrow(inb1_spp))	)
			community.composition <- as.data.frame(cbind(comm_comp_idx, comm_comp_spp))
			colnames(community.composition)[1] <- "community"
			
			FD <- Getlength(xtree, community.composition)
			FD_c <- FD[2]
			FD_r_tmp1 <- t(FD_c)
			FD_r_tmp2 <- FD_r_tmp1[-1]
			
			if(any(FD_new_0s == 0) == TRUE){
				ins_vals <- rep(0, length(FD_new_0s_vec))
				FD_r <- R.utils::insert(FD_r_tmp2, ats = FD_new_0s_vec, values = ins_vals)
				} else {FD_r <- FD_r_tmp2}


	}  # end 'Current diversity' if loop
  

	# Assemble the results vector for iteration i
	res_vec <- c(SR0, SR_allthreats, SR_hab1, SR_hunt1, SR_clim1, SR_con1, SR_nonnat1, SR_pol1,  
			SR_hyb1, SR_prey1, SR_dis1, SR_inb1,  
 
       	      FD0, 
		 	FD_r,  # length 11

               	FDrand_allthreats, FDrand_hab1, FDrand_hunt1, FDrand_clim1, FDrand_con1, FDrand_nonnat1,  
			   FDrand_pol1, FDrand_hyb1, FDrand_prey1, FDrand_dis1, FDrand_inb1, 	
 
               	FDrandLO_allthreats, FDrandLO_hab1, FDrandLO_hunt1, FDrandLO_clim1, FDrandLO_con1, FDrandLO_nonnat1,  
			   FDrandLO_pol1, FDrandLO_hyb1, FDrandLO_prey1, FDrandLO_dis1, FDrandLO_inb1, 	
 
               	FDrandHI_allthreats, FDrandHI_hab1, FDrandHI_hunt1, FDrandHI_clim1, FDrandHI_con1, FDrandHI_nonnat1,  
			   FDrandHI_pol1, FDrandHI_hyb1, FDrandHI_prey1, FDrandHI_dis1, FDrandHI_inb1, 	

              	FDrandLObonferroni_allthreats, FDrandLObonferroni_hab1, FDrandLObonferroni_hunt1, FDrandLObonferroni_clim1, FDrandLObonferroni_con1, FDrandLObonferroni_nonnat1,  
			   FDrandLObonferroni_pol1, FDrandLObonferroni_hyb1, FDrandLObonferroni_prey1, FDrandLObonferroni_dis1, FDrandLObonferroni_inb1, 	
 
            	FDrandHIbonferroni_allthreats, FDrandHIbonferroni_hab1, FDrandHIbonferroni_hunt1, FDrandHIbonferroni_clim1, FDrandHIbonferroni_con1, FDrandHIbonferroni_nonnat1,  
			   FDrandHIbonferroni_pol1, FDrandHIbonferroni_hyb1, FDrandHIbonferroni_prey1, FDrandHIbonferroni_dis1, FDrandHIbonferroni_inb1, 	

               	PD0, 
			PD_allthreats/PD0, PD_hab1/PD0, PD_hunt1/PD0, PD_clim1/PD0, PD_con1/PD0,  
			   PD_nonnat1/PD0, PD_pol1/PD0, PD_hyb1/PD0, PD_prey1/PD0, PD_dis1/PD0, PD_inb1/PD0,  

               	PDrand_allthreats/PD0, PDrand_hab1/PD0, PDrand_hunt1/PD0, PDrand_clim1/PD0, PDrand_con1/PD0, PDrand_nonnat1/PD0,  
			   PDrand_pol1/PD0, PDrand_hyb1/PD0, PDrand_prey1/PD0, PDrand_dis1/PD0, PDrand_inb1/PD0, 	
 
               	PDrandLO_allthreats/PD0, PDrandLO_hab1/PD0, PDrandLO_hunt1/PD0, PDrandLO_clim1/PD0, PDrandLO_con1/PD0, PDrandLO_nonnat1/PD0,  
			   PDrandLO_pol1/PD0, PDrandLO_hyb1/PD0, PDrandLO_prey1/PD0, PDrandLO_dis1/PD0, PDrandLO_inb1/PD0, 	
 
               	PDrandHI_allthreats/PD0, PDrandHI_hab1/PD0, PDrandHI_hunt1/PD0, PDrandHI_clim1/PD0, PDrandHI_con1/PD0, PDrandHI_nonnat1/PD0,  
			   PDrandHI_pol1/PD0, PDrandHI_hyb1/PD0, PDrandHI_prey1/PD0, PDrandHI_dis1/PD0, PDrandHI_inb1/PD0, 	

              	PDrandLObonferroni_allthreats/PD0, PDrandLObonferroni_hab1/PD0, PDrandLObonferroni_hunt1/PD0, PDrandLObonferroni_clim1/PD0, PDrandLObonferroni_con1/PD0, PDrandLObonferroni_nonnat1/PD0,  
			   PDrandLObonferroni_pol1/PD0, PDrandLObonferroni_hyb1/PD0, PDrandLObonferroni_prey1/PD0, PDrandLObonferroni_dis1/PD0, PDrandLObonferroni_inb1/PD0, 	
 
            	PDrandHIbonferroni_allthreats/PD0, PDrandHIbonferroni_hab1/PD0, PDrandHIbonferroni_hunt1/PD0, PDrandHIbonferroni_clim1/PD0, PDrandHIbonferroni_con1/PD0, PDrandHIbonferroni_nonnat1/PD0,  
			   PDrandHIbonferroni_pol1/PD0, PDrandHIbonferroni_hyb1/PD0, PDrandHIbonferroni_prey1/PD0, PDrandHIbonferroni_dis1/PD0, PDrandHIbonferroni_inb1/PD0, 	
 
		      thab.01-t0.01, thab.02-t0.02, thab.03-t0.03, thab.04-t0.04, thab.05-t0.05, thab.06-t0.06,
		      thab.07-t0.07, thab.08-t0.08, thab.09-t0.09, thab.10-t0.10, thab.11-t0.11, thab.12-t0.12,
		      thunt.01-t0.01, thunt.02-t0.02, thunt.03-t0.03, thunt.04-t0.04, thunt.05-t0.05, thunt.06-t0.06,
		      thunt.07-t0.07, thunt.08-t0.08, thunt.09-t0.09, thunt.10-t0.10, thunt.11-t0.11, thunt.12-t0.12)
  
	if(res_vec[1] == 0){res_vec[] <- 0} else {res_vec <- res_vec}


out_mat[i,] <- res_vec
	
	
rm("datXY", "all_spp", "country_name", "keep_exc", "spp_exc", "spp_exc_df", "spp_list", 
	"SR0", "allthreats_spp", "hab1_spp", "hunt1_spp", "clim1_spp", "con1_spp", "nonnat1_spp", 
	"pol1_spp", "prey1_spp", "hyb1_spp", "dis1_spp", "inb1_spp",  
	 
	"numspploss_allthreats", "numspploss_hab1", "numspploss_hunt1", "numspploss_clim1", "numspploss_con1", 
	"numspploss_nonnat1", "numspploss_pol1", "numspploss_prey1", "numspploss_hyb1",  "numspploss_dis1", "numspploss_inb1", 	

	"SR_allthreats", "SR_hab1", "SR_hunt1", "SR_clim1", "SR_con1", "SR_nonnat1", "SR_pol1", "SR_hyb1", 
	"SR_prey1", "SR_dis1", "SR_inb1", 

	"i.prime", "tree", "xtree", "distances", 
	"comm_comp_idx" , "comm_comp_spp", "community.composition", "species.traits",

	"FD0", "FD_r", "FD_r_tmp1", "FD_r_tmp2", "FD_c", "FD", "FD_new_0s", "FD_new_0s_vec",

	"FDrand_allthreats", "FDrand_hab1", "FDrand_hunt1", "FDrand_clim1", "FDrand_con1", 
	"FDrand_nonnat1", "FDrand_pol1", "FDrand_hyb1", "FDrand_prey1", "FDrand_dis1", "FDrand_inb1", 
	
	"FDrandLO_allthreats", "FDrandLO_hab1", "FDrandLO_hunt1", "FDrandLO_clim1", "FDrandLO_con1", 
	"FDrandLO_nonnat1", "FDrandLO_pol1", "FDrandLO_hyb1", "FDrandLO_prey1", "FDrandLO_dis1", "FDrandLO_inb1", 
	
	"FDrandHI_allthreats", "FDrandHI_hab1", "FDrandHI_hunt1", "FDrandHI_clim1", "FDrandHI_con1", 
	"FDrandHI_nonnat1", "FDrandHI_pol1", "FDrandHI_hyb1", "FDrandHI_prey1", "FDrandHI_dis1", "FDrandHI_inb1", 
	
	"FDrandLObonferroni_allthreats", "FDrandLObonferroni_hab1", "FDrandLObonferroni_hunt1", "FDrandLObonferroni_clim1", "FDrandLObonferroni_con1", 
	"FDrandLObonferroni_nonnat1", "FDrandLObonferroni_pol1", "FDrandLObonferroni_hyb1", "FDrandLObonferroni_prey1", "FDrandLObonferroni_dis1", "FDrandLObonferroni_inb1", 
	
	"FDrandHIbonferroni_allthreats", "FDrandHIbonferroni_hab1", "FDrandHIbonferroni_hunt1", "FDrandHIbonferroni_clim1", "FDrandHIbonferroni_con1", 
	"FDrandHIbonferroni_nonnat1", "FDrandHIbonferroni_pol1", "FDrandHIbonferroni_hyb1", "FDrandHIbonferroni_prey1", "FDrandHIbonferroni_dis1", "FDrandHIbonferroni_inb1", 

	"PD0", "PD_allthreats", "PD_hab1", "PD_hunt1", "PD_clim1", "PD_con1", "PD_nonnat1", "PD_pol1", 
	"PD_hyb1", "PD_prey1", "PD_dis1", "PD_inb1", 

	"PDrand_allthreats", "PDrand_hab1", "PDrand_hunt1", "PDrand_clim1", "PDrand_con1", 
	"PDrand_nonnat1", "PDrand_pol1", "PDrand_hyb1", "PDrand_prey1", "PDrand_dis1", "PDrand_inb1", 
	
	"PDrandLO_allthreats", "PDrandLO_hab1", "PDrandLO_hunt1", "PDrandLO_clim1", "PDrandLO_con1", 
	"PDrandLO_nonnat1", "PDrandLO_pol1", "PDrandLO_hyb1", "PDrandLO_prey1", "PDrandLO_dis1", "PDrandLO_inb1", 
	
	"PDrandHI_allthreats", "PDrandHI_hab1", "PDrandHI_hunt1", "PDrandHI_clim1", "PDrandHI_con1", 
	"PDrandHI_nonnat1", "PDrandHI_pol1", "PDrandHI_hyb1", "PDrandHI_prey1", "PDrandHI_dis1", "PDrandHI_inb1", 
	
	"PDrandLObonferroni_allthreats", "PDrandLObonferroni_hab1", "PDrandLObonferroni_hunt1", "PDrandLObonferroni_clim1", "PDrandLObonferroni_con1", 
	"PDrandLObonferroni_nonnat1", "PDrandLObonferroni_pol1", "PDrandLObonferroni_hyb1", "PDrandLObonferroni_prey1", "PDrandLObonferroni_dis1", "PDrandLObonferroni_inb1", 
	
	"PDrandHIbonferroni_allthreats", "PDrandHIbonferroni_hab1", "PDrandHIbonferroni_hunt1", "PDrandHIbonferroni_clim1", "PDrandHIbonferroni_con1", 
	"PDrandHIbonferroni_nonnat1", "PDrandHIbonferroni_pol1", "PDrandHIbonferroni_hyb1", "PDrandHIbonferroni_prey1", "PDrandHIbonferroni_dis1", "PDrandHIbonferroni_inb1", 

	"t0.01", "t0.02", "t0.03", "t0.04", "t0.05", "t0.06", "t0.07", "t0.08", 
	"t0.09", "t0.10", "t0.11", "t0.12", "thab.01", "thab.02", "thab.03", "thab.04", "thab.05", "thab.06", 
	"thab.07", "thab.08", "thab.09", "thab.10", "thab.11", "thab.12", "thunt.01", "thunt.02", "thunt.03", 
	"thunt.04", "thunt.05", "thunt.06", "thunt.07", "thunt.08", "thunt.09", "thunt.10", "thunt.11", "thunt.12", 
	"hab1_tmp", "hunt1_tmp", "res_vec")
	
	
}, error=function(e){}) # end of the tryCatch function
}  # end i loop





#----------- ASSEMBLE AND SAVE THE RESULTS MATRIX -------------------------------
out_df <- as.data.frame(out_mat)
res_df <- cbind.data.frame(ctrs_area, out_df)
names(res_df) <- c("Country", "Area_sq_km",
                  "SR0", "SR.allthreats", "SR.habitat1", "SR.hunting1", "SR.climate1", "SR.conflict1", 
                  "SR.nonNatives1", "SR.pollution1", "SR.hybrid1", "SR.prey1", "SR.disease1", "SR.inbreeding1", 
         
		      "FD0", 
			"FD.allthreats", "FD.habitat1", "FD.hunting1", "FD.climate1", "FD.conflict1", 
                  "FD.nonNatives1", "FD.pollution1", "FD.hybrid1", "FD.prey1", "FD.disease1", "FD.inbreeding1", 
 
                  "FDrand.allthreats", "FDrand.habitat1", "FDrand.hunting1", "FDrand.climate1", "FDrand.conflict1", 
			"FDrand.nonNatives1", "FDrand.pollution1", "FDrand.hybrid1", "FDrand.prey1", "FDrand.disease1", "FDrand.inbreeding1", 

                  "FDrandLO.allthreats", "FDrandLO.habitat1", "FDrandLO.hunting1", "FDrandLO.climate1", "FDrandLO.conflict1", 
			"FDrandLO.nonNatives1", "FDrandLO.pollution1", "FDrandLO.hybrid1", "FDrandLO.prey1", "FDrandLO.disease1", "FDrandLO.inbreeding1", 

                  "FDrandHI.allthreats", "FDrandHI.habitat1", "FDrandHI.hunting1", "FDrandHI.climate1", "FDrandHI.conflict1", 
			"FDrandHI.nonNatives1", "FDrandHI.pollution1", "FDrandHI.hybrid1", "FDrandHI.prey1", "FDrandHI.disease1", "FDrandHI.inbreeding1", 

                  "FDrandLObonferroni.allthreats", "FDrandLObonferroni.habitat1", "FDrandLObonferroni.hunting1", "FDrandLObonferroni.climate1", "FDrandLObonferroni.conflict1", 
			"FDrandLObonferroni.nonNatives1", "FDrandLObonferroni.pollution1", "FDrandLObonferroni.hybrid1", "FDrandLObonferroni.prey1", "FDrandLObonferroni.disease1", "FDrandLObonferroni.inbreeding1", 

                  "FDrandHIbonferroni.allthreats", "FDrandHIbonferroni.habitat1", "FDrandHIbonferroni.hunting1", "FDrandHIbonferroni.climate1", "FDrandHIbonferroni.conflict1", 
			"FDrandHIbonferroni.nonNatives1", "FDrandHIbonferroni.pollution1", "FDrandHIbonferroni.hybrid1", "FDrandHIbonferroni.prey1", "FDrandHIbonferroni.disease1", "FDrandHIbonferroni.inbreeding1", 

                  "PD0", "PD.allthreats", "PD.habitat1", "PD.hunting1", "PD.climate1", "PD.conflict1", 
                  "PD.nonNatives1", "PD.pollution1", "PD.hybrid1", "PD.prey1", "PD.disease1", "PD.inbreeding1", 

                  "PDrand.allthreats", "PDrand.habitat1", "PDrand.hunting1", "PDrand.climate1", "PDrand.conflict1", 
			"PDrand.nonNatives1", "PDrand.pollution1", "PDrand.hybrid1", "PDrand.prey1", "PDrand.disease1", "PDrand.inbreeding1", 

                  "PDrandLO.allthreats", "PDrandLO.habitat1", "PDrandLO.hunting1", "PDrandLO.climate1", "PDrandLO.conflict1", 
			"PDrandLO.nonNatives1", "PDrandLO.pollution1", "PDrandLO.hybrid1", "PDrandLO.prey1", "PDrandLO.disease1", "PDrandLO.inbreeding1", 

                  "PDrandHI.allthreats", "PDrandHI.habitat1", "PDrandHI.hunting1", "PDrandHI.climate1", "PDrandHI.conflict1", 
			"PDrandHI.nonNatives1", "PDrandHI.pollution1", "PDrandHI.hybrid1", "PDrandHI.prey1", "PDrandHI.disease1", "PDrandHI.inbreeding1", 

                  "PDrandLObonferroni.allthreats", "PDrandLObonferroni.habitat1", "PDrandLObonferroni.hunting1", "PDrandLObonferroni.climate1", "PDrandLObonferroni.conflict1", 
			"PDrandLObonferroni.nonNatives1", "PDrandLObonferroni.pollution1", "PDrandLObonferroni.hybrid1", "PDrandLObonferroni.prey1", "PDrandLObonferroni.disease1", "PDrandLObonferroni.inbreeding1", 

                  "PDrandHIbonferroni.allthreats", "PDrandHIbonferroni.habitat1", "PDrandHIbonferroni.hunting1", "PDrandHIbonferroni.climate1", "PDrandHIbonferroni.conflict1", 
			"PDrandHIbonferroni.nonNatives1", "PDrandHIbonferroni.pollution1", "PDrandHIbonferroni.hybrid1", "PDrandHIbonferroni.prey1", "PDrandHIbonferroni.disease1", "PDrandHIbonferroni.inbreeding1", 

			"TraitChange.hab1.DietInv", "TraitChange.hab1.DietVert", "TraitChange.hab1.DietFish", "TraitChange.hab1.DietScav", 
			"TraitChange.hab1.DietFruit", "TraitChange.hab1.DietNect", "TraitChange.hab1.DietSeed", "TraitChange.hab1.DietHerb", 
			"TraitChange.hab1.BodyMass", "TraitChange.hab1.Ground", "TraitChange.hab1.Climbing", "TraitChange.hab1.Volant", 
			"TraitChange.hunt1.DietInv", "TraitChange.hunt1.DietVert", "TraitChange.hunt1.DietFish", "TraitChange.hunt1.DietScav", 
			"TraitChange.hunt1.DietFruit", "TraitChange.hunt1.DietNect", "TraitChange.hunt1.DietSeed", "TraitChange.hunt1.DietHerb", 
			"TraitChange.hunt1.BodyMass", "TraitChange.hunt1.Ground", "TraitChange.hunt1.Climbing", "TraitChange.hunt1.Volant")
	





# Add columns for 'bias' in FD and PD loss 
# Value explanations:
#  1: loss is significantly biased towards functionally (or phylogenetically) unique species
#  0: loss is unbiased
# -1: loss is significantly biased towards functionally (or phylogenetically) redundant species

dat0 <- res_df

#--- Functional diversity
# FD: all threats
dat0$tmp <- ifelse(dat0$FD.allthreats <= dat0$FDrandLO.allthreats, 1, 0)
dat0$tmp <- ifelse(dat0$FD.allthreats >= dat0$FDrandHI.allthreats, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.allthreats, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.allthreats <- as.factor(dat0$tmp)

# FD: habitat 1
dat0$tmp <- ifelse(dat0$FD.habitat1 <= dat0$FDrandLO.habitat1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.habitat1 >= dat0$FDrandHI.habitat1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.habitat1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.habitat1 <- as.factor(dat0$tmp)

# FD: hunting 1
dat0$tmp <- ifelse(dat0$FD.hunting1 <= dat0$FDrandLO.hunting1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.hunting1 >= dat0$FDrandHI.hunting1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hunting1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.hunting1 <- as.factor(dat0$tmp)

# FD: climate 1
dat0$tmp <- ifelse(dat0$FD.climate1 <= dat0$FDrandLO.climate1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.climate1 >= dat0$FDrandHI.climate1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.climate1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.climate1 <- as.factor(dat0$tmp)

# FD: conflict 1
dat0$tmp <- ifelse(dat0$FD.conflict1 <= dat0$FDrandLO.conflict1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.conflict1 >= dat0$FDrandHI.conflict1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.conflict1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.conflict1 <- as.factor(dat0$tmp)

# FD: non-natives 1
dat0$tmp <- ifelse(dat0$FD.nonNatives1 <= dat0$FDrandLO.nonNatives1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.nonNatives1 >= dat0$FDrandHI.nonNatives1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.nonNatives1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.nonNatives1 <- as.factor(dat0$tmp)

# FD: pollution 1
dat0$tmp <- ifelse(dat0$FD.pollution1 <= dat0$FDrandLO.pollution1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.pollution1 >= dat0$FDrandHI.pollution1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.pollution1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.pollution1 <- as.factor(dat0$tmp)

# FD: hybridization 1
dat0$tmp <- ifelse(dat0$FD.hybrid1 <= dat0$FDrandLO.hybrid1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.hybrid1 >= dat0$FDrandHI.hybrid1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hybrid1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.hybrid1 <- as.factor(dat0$tmp)

# FD: prey depletion 1
dat0$tmp <- ifelse(dat0$FD.prey1 <= dat0$FDrandLO.prey1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.prey1 >= dat0$FDrandHI.prey1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.prey1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.prey1 <- as.factor(dat0$tmp)

# FD: disease 1
dat0$tmp <- ifelse(dat0$FD.disease1 <= dat0$FDrandLO.disease1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.disease1 >= dat0$FDrandHI.disease1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.disease1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.disease1 <- as.factor(dat0$tmp)

# FD: inbreeding 1
dat0$tmp <- ifelse(dat0$FD.inbreeding1 <= dat0$FDrandLO.inbreeding1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.inbreeding1 >= dat0$FDrandHI.inbreeding1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.inbreeding1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbias.inbreeding1 <- as.factor(dat0$tmp)


#--- Phylogenetic diversity
# PD: all threats
dat0$tmp <- ifelse(dat0$PD.allthreats <= dat0$PDrandLO.allthreats, 1, 0)
dat0$tmp <- ifelse(dat0$PD.allthreats >= dat0$PDrandHI.allthreats, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.allthreats, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.allthreats <- as.factor(dat0$tmp)

# PD: habitat 1
dat0$tmp <- ifelse(dat0$PD.habitat1 <= dat0$PDrandLO.habitat1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.habitat1 >= dat0$PDrandHI.habitat1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.habitat1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.habitat1 <- as.factor(dat0$tmp)

# PD: hunting 1
dat0$tmp <- ifelse(dat0$PD.hunting1 <= dat0$PDrandLO.hunting1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.hunting1 >= dat0$PDrandHI.hunting1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hunting1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.hunting1 <- as.factor(dat0$tmp)

# PD: climate 1
dat0$tmp <- ifelse(dat0$PD.climate1 <= dat0$PDrandLO.climate1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.climate1 >= dat0$PDrandHI.climate1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.climate1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.climate1 <- as.factor(dat0$tmp)

# PD: conflict 1
dat0$tmp <- ifelse(dat0$PD.conflict1 <= dat0$PDrandLO.conflict1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.conflict1 >= dat0$PDrandHI.conflict1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.conflict1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.conflict1 <- as.factor(dat0$tmp)

# PD: non-natives 1
dat0$tmp <- ifelse(dat0$PD.nonNatives1 <= dat0$PDrandLO.nonNatives1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.nonNatives1 >= dat0$PDrandHI.nonNatives1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.nonNatives1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.nonNatives1 <- as.factor(dat0$tmp)

# PD: pollution 1
dat0$tmp <- ifelse(dat0$PD.pollution1 <= dat0$PDrandLO.pollution1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.pollution1 >= dat0$PDrandHI.pollution1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.pollution1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.pollution1 <- as.factor(dat0$tmp)

# PD: hybridization 1
dat0$tmp <- ifelse(dat0$PD.hybrid1 <= dat0$PDrandLO.hybrid1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.hybrid1 >= dat0$PDrandHI.hybrid1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hybrid1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.hybrid1 <- as.factor(dat0$tmp)

# PD: prey depletion 1
dat0$tmp <- ifelse(dat0$PD.prey1 <= dat0$PDrandLO.prey1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.prey1 >= dat0$PDrandHI.prey1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.prey1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.prey1 <- as.factor(dat0$tmp)

# PD: disease 1
dat0$tmp <- ifelse(dat0$PD.disease1 <= dat0$PDrandLO.disease1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.disease1 >= dat0$PDrandHI.disease1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.disease1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.disease1 <- as.factor(dat0$tmp)

# PD: inbreeding 1
dat0$tmp <- ifelse(dat0$PD.inbreeding1 <= dat0$PDrandLO.inbreeding1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.inbreeding1 >= dat0$PDrandHI.inbreeding1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.inbreeding1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbias.inbreeding1 <- as.factor(dat0$tmp)


#--- Functional diversity --- with Bonferroni correction
# FD: all threats
dat0$tmp <- ifelse(dat0$FD.allthreats <= dat0$FDrandLObonferroni.allthreats, 1, 0)
dat0$tmp <- ifelse(dat0$FD.allthreats >= dat0$FDrandHIbonferroni.allthreats, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.allthreats, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.allthreats <- as.factor(dat0$tmp)

# FD: habitat 1
dat0$tmp <- ifelse(dat0$FD.habitat1 <= dat0$FDrandLObonferroni.habitat1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.habitat1 >= dat0$FDrandHIbonferroni.habitat1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.habitat1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.habitat1 <- as.factor(dat0$tmp)

# FD: hunting 1
dat0$tmp <- ifelse(dat0$FD.hunting1 <= dat0$FDrandLObonferroni.hunting1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.hunting1 >= dat0$FDrandHIbonferroni.hunting1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hunting1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.hunting1 <- as.factor(dat0$tmp)

# FD: climate 1
dat0$tmp <- ifelse(dat0$FD.climate1 <= dat0$FDrandLObonferroni.climate1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.climate1 >= dat0$FDrandHIbonferroni.climate1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.climate1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.climate1 <- as.factor(dat0$tmp)

# FD: conflict 1
dat0$tmp <- ifelse(dat0$FD.conflict1 <= dat0$FDrandLObonferroni.conflict1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.conflict1 >= dat0$FDrandHIbonferroni.conflict1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.conflict1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.conflict1 <- as.factor(dat0$tmp)

# FD: non-natives 1
dat0$tmp <- ifelse(dat0$FD.nonNatives1 <= dat0$FDrandLObonferroni.nonNatives1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.nonNatives1 >= dat0$FDrandHIbonferroni.nonNatives1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.nonNatives1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.nonNatives1 <- as.factor(dat0$tmp)

# FD: pollution 1
dat0$tmp <- ifelse(dat0$FD.pollution1 <= dat0$FDrandLObonferroni.pollution1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.pollution1 >= dat0$FDrandHIbonferroni.pollution1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.pollution1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.pollution1 <- as.factor(dat0$tmp)

# FD: hybridization 1
dat0$tmp <- ifelse(dat0$FD.hybrid1 <= dat0$FDrandLObonferroni.hybrid1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.hybrid1 >= dat0$FDrandHIbonferroni.hybrid1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hybrid1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.hybrid1 <- as.factor(dat0$tmp)

# FD: prey depletion 1
dat0$tmp <- ifelse(dat0$FD.prey1 <= dat0$FDrandLObonferroni.prey1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.prey1 >= dat0$FDrandHIbonferroni.prey1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.prey1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.prey1 <- as.factor(dat0$tmp)

# FD: disease 1
dat0$tmp <- ifelse(dat0$FD.disease1 <= dat0$FDrandLObonferroni.disease1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.disease1 >= dat0$FDrandHIbonferroni.disease1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.disease1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.disease1 <- as.factor(dat0$tmp)

# FD: inbreeding 1
dat0$tmp <- ifelse(dat0$FD.inbreeding1 <= dat0$FDrandLObonferroni.inbreeding1, 1, 0)
dat0$tmp <- ifelse(dat0$FD.inbreeding1 >= dat0$FDrandHIbonferroni.inbreeding1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.inbreeding1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$FDbiasbonferroni.inbreeding1 <- as.factor(dat0$tmp)


#--- Phylogenetic diversity --- with Bonferroni correction
# PD: all threats
dat0$tmp <- ifelse(dat0$PD.allthreats <= dat0$PDrandLObonferroni.allthreats, 1, 0)
dat0$tmp <- ifelse(dat0$PD.allthreats >= dat0$PDrandHIbonferroni.allthreats, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.allthreats, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.allthreats <- as.factor(dat0$tmp)

# PD: habitat 1
dat0$tmp <- ifelse(dat0$PD.habitat1 <= dat0$PDrandLObonferroni.habitat1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.habitat1 >= dat0$PDrandHIbonferroni.habitat1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.habitat1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.habitat1 <- as.factor(dat0$tmp)

# PD: hunting 1
dat0$tmp <- ifelse(dat0$PD.hunting1 <= dat0$PDrandLObonferroni.hunting1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.hunting1 >= dat0$PDrandHIbonferroni.hunting1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hunting1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.hunting1 <- as.factor(dat0$tmp)

# PD: climate 1
dat0$tmp <- ifelse(dat0$PD.climate1 <= dat0$PDrandLObonferroni.climate1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.climate1 >= dat0$PDrandHIbonferroni.climate1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.climate1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.climate1 <- as.factor(dat0$tmp)

# PD: conflict 1
dat0$tmp <- ifelse(dat0$PD.conflict1 <= dat0$PDrandLObonferroni.conflict1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.conflict1 >= dat0$PDrandHIbonferroni.conflict1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.conflict1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.conflict1 <- as.factor(dat0$tmp)

# PD: non-natives 1
dat0$tmp <- ifelse(dat0$PD.nonNatives1 <= dat0$PDrandLObonferroni.nonNatives1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.nonNatives1 >= dat0$PDrandHIbonferroni.nonNatives1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.nonNatives1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.nonNatives1 <- as.factor(dat0$tmp)

# PD: pollution 1
dat0$tmp <- ifelse(dat0$PD.pollution1 <= dat0$PDrandLObonferroni.pollution1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.pollution1 >= dat0$PDrandHIbonferroni.pollution1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.pollution1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.pollution1 <- as.factor(dat0$tmp)

# PD: hybridization 1
dat0$tmp <- ifelse(dat0$PD.hybrid1 <= dat0$PDrandLObonferroni.hybrid1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.hybrid1 >= dat0$PDrandHIbonferroni.hybrid1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.hybrid1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.hybrid1 <- as.factor(dat0$tmp)

# PD: prey depletion 1
dat0$tmp <- ifelse(dat0$PD.prey1 <= dat0$PDrandLObonferroni.prey1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.prey1 >= dat0$PDrandHIbonferroni.prey1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.prey1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.prey1 <- as.factor(dat0$tmp)

# PD: disease 1
dat0$tmp <- ifelse(dat0$PD.disease1 <= dat0$PDrandLObonferroni.disease1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.disease1 >= dat0$PDrandHIbonferroni.disease1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.disease1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.disease1 <- as.factor(dat0$tmp)

# PD: inbreeding 1
dat0$tmp <- ifelse(dat0$PD.inbreeding1 <= dat0$PDrandLObonferroni.inbreeding1, 1, 0)
dat0$tmp <- ifelse(dat0$PD.inbreeding1 >= dat0$PDrandHIbonferroni.inbreeding1, -1, dat0$tmp)
dat0$tmp <- ifelse(dat0$SR0 == dat0$SR.inbreeding1, 0, dat0$tmp)
dat0$tmp[is.na(dat0$tmp)] <- 0
dat0$PDbiasbonferroni.inbreeding1 <- as.factor(dat0$tmp)


dat0 <- subset(dat0, select = -c(tmp))


#------------------- SAVE FILE ------------------------------------------------------
write.csv(dat0, "FuncPhylo_mammals190816.csv", row.names=F) 

end <- Sys.time()
end
elapsed <- end - start
elapsed  # 2.3 days to run on my laptop with numrand=301; 3.5 days on desktop with numrand=501

















