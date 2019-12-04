###############################################################################
#  Functional diversity analysis plots
#  June 15, 2019; last updated July 2, 2019
#  Script to read in results dataframe from diversity analysis and generate maps
#  Sara Williams
###############################################################################



# =============================================================================
#  Load packages and set working director.
# =============================================================================
library(dplyr)
library(sf)
library(tidyr)
library(raster)
library(ggplot2)
library(ggpubr)
library(smoothr)
library(cowplot)
library(extrafont)


setwd("C:/Users/saraw/Documents/FD/")
st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))


# =============================================================================
#  Load input data and set base parameters
# =============================================================================
	
	# ----------------------
	#	Projection used in analysis
	prj <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
	
	# ----------------------
	#  Countries
	ctrs <- st_read("Data/Raw/UIA_World_Countries_Boundaries.shp") %>%
		dplyr::select(Country, geometry) %>%
		st_transform(st_crs(prj))
	
	# ----------------------
	#  Continents outline
	area_thresh <- units::set_units(500, km^2)
	con_crop <- st_read("Data/Raw/p_template.shp") %>%
		st_transform(prj) %>%
		drop_crumbs(threshold = area_thresh) 
	
	# ----------------------
	#  Make continents nice for plotting
	area_thresh2 <- units::set_units(50000, km^2)
	con_out <- con_crop %>%
		st_buffer(30000) #%>%
		#fill_holes(area_thresh2)
		
	# ----------------------
	#	Read in results CSV and convert to sf
	#    Save for future use
	res_df <- read.csv("Results/......")
	#  Input latest results file here
	res_sf <- res_df %>% 
		left_join(ctrs, by = "Country") %>%
		st_as_sf(sf_column_name = "geometry", crs = prj) %>%
		mutate(ctry_sq_km = as.numeric(st_area(.))/1000000)  #%>%
		#dplyr::filter(ctry_sq_km > 999) 

	
	# ----------------------
	#  Set up vector containing plot titles
	p_titles_vec <- c("habitat_loss", "hunting", "climate_change", 
		"human_wildlife_conflict", "non_natives", "pollution", 
		"hybridization", "prey_depletion", "disease", "inbreeding")
		

# =============================================================================
#  Faceted current diversity plots
# =============================================================================
	
	# ----------------------
	#  Baseline - current taxonomic, functional and phylogenetic diversity
	curr_plot_cols <- c(3, 24, 45)
	
	# ----------------------
	#  Prep plotting data 
	curr_res_sf <- res_sf[curr_plot_cols]
	
	td_mu <- mean(curr_res_sf$SR0, na.rm = TRUE)
	td_sd <- sd(curr_res_sf$SR0, na.rm = TRUE)
	fd_mu <- mean(curr_res_sf$FD0, na.rm = TRUE)
	fd_sd <- sd(curr_res_sf$FD0, na.rm = TRUE)
	pd_mu <- mean(curr_res_sf$PD0, na.rm = TRUE)
	pd_sd <- sd(curr_res_sf$PD0, na.rm = TRUE)
	
	colnames(curr_res_sf)[1:3] <- c("Current Taxonomic Diversity", 
		"Current Functional Diversity",
		"Current Phylogenetic Diversity")
	curr_res_sf_long <- curr_res_sf %>% gather(richness, value, -geometry) 
	td_tmp <- curr_res_sf_long %>%
		dplyr::filter(richness == "Current Taxonomic Diversity") %>%
		mutate(stand_val = (value - td_mu) / td_sd)
	fd_tmp <- curr_res_sf_long %>%
		dplyr::filter(richness == "Current Functional Diversity") %>%
		mutate(stand_val = (value - fd_mu) / fd_sd)
	pd_tmp <- curr_res_sf_long %>%
		dplyr::filter(richness == "Current Phylogenetic Diversity") %>%
		mutate(stand_val = (value - pd_mu) / pd_sd)
	curr_res_sf_long_stand <- td_tmp %>%
		rbind(fd_tmp) %>%
		rbind(pd_tmp)
		
		
	
	tiff("current_div_facet.tiff", width = 6, height = 9, units = 'in', res = 300)
	# ----------------------
	#  Faceted plot
	p_curr_div <- ggplot() +
		geom_sf(data = curr_res_sf_long_stand, aes(colour = stand_val, fill = stand_val)) +
		theme(plot.background = element_blank(),
			legend.position = "right", 
			legend.direction = "vertical",
			legend.spacing.y = unit(0.25, 'cm'),
			panel.background = element_blank(), 
			panel.grid.major = element_line(colour = "transparent"),
			panel.grid.minor = element_line(colour = "transparent"),
			panel.border = element_blank(),
			strip.placement = "outside", 
			strip.background = element_rect("transparent"),
			panel.spacing = unit(0.5, "cm"), 
			plot.title = element_blank(),
			legend.title = element_text(family = "Calibri", size = 10, face="italic"),
			legend.text = element_text(family = "Calibri", size = 9),
			strip.text = element_text(family = "Calibri", size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
			legend.background = element_rect(fill="white", colour = "black")) +
		#geom_sf(data = con_out, fill = "transparent", color = "grey70", size = 0.4) +
		facet_wrap(vars(richness), ncol = 1, strip.position = "top") +
		scale_colour_gradient(high = "forestgreen", low ="white", 
			name = "Standardized \nrichness metric", na.value = "transparent") +
		scale_fill_gradient(high = "forestgreen", low ="white", 
			guide = FALSE, na.value = "black") +
		guides(colour = guide_colorbar(ticks = FALSE)) +
		coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
	p_curr_div
	dev.off()

	
	
# =============================================================================
#  Functional Diversity
# =============================================================================

	# ----------------------
	#  Raw
	func_plot_cols <- c(24, seq(26, 44, by = 2))
	fd_res_sf <- res_sf[func_plot_cols]
	colnames(fd_res_sf)[1] <- "current"
	
	# ----------------------
	#  Prep plotting data for raw value decline
	raw_fd_decline_sf <- fd_res_sf %>%
		mutate_at(vars(FD.habitat2:FD.inbreeding2),  ~(current - (current * .))) %>%
		dplyr::select(-current) 

	# ----------------------
	#  Prep plotting data for proportional decline
	prop_fd_decline_sf <- fd_res_sf %>%
		mutate_at(vars(FD.habitat2:FD.inbreeding2),  ~(ifelse(current == 0, 0, 1 - .))) %>%
		dplyr::select(-current) 

	
	for(i in 1:10){
		
		# ----------------------
		#  Set up column that is used for raw value data
		col_var_raw <- colnames(raw_fd_decline_sf[i])[1]
		
		# ----------------------
		#  Set up column that is used for proportional data
		col_var_prop <- colnames(prop_fd_decline_sf[i])[1]
		
		# ----------------------
		#  Color scale ranges
		tmp_prop <- prop_fd_decline_sf %>%
			dplyr::select(col_var_prop)
		high_prop_fd_decline_sf <- tmp_prop %>%
			dplyr::filter_at(vars(-geometry), any_vars(. > 0.69999))
		rest_prop_fd_decline_sf <- tmp_prop %>%
			dplyr::filter_at(vars(-geometry), any_vars(. < 0.699991))
		
		
		# ----------------------
		#  Generate plots
		raw_FD <- ggplot() +
			geom_sf(data = raw_fd_decline_sf, 
				aes(colour = get(noquote(col_var_prop)),  
					fill = get(noquote(col_var_prop))), 
				size = 0.4) +
			scale_colour_gradient(high = "red", 
				low ="white", 
				name = "Raw              \ndecline", 
				na.value = "transparent", 
				limits=c(0, 70), 
				guide = FALSE) +
			scale_fill_gradient(high = "red", 
				low ="white", 
				name = "Raw              \ndecline", 
				na.value = "transparent", 
				limits=c(0, 70),
				breaks = c(10, 30, 50, 70),
				labels = c("10", "30", "50", "70")) +
			guides(fill = guide_colorbar(ticks = FALSE)) +
			geom_sf(data = con_out, fill = "transparent", color = "grey70", size = 0.3) +
			coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) +
			theme(plot.background = element_blank(),
				legend.title = element_text(family = "Calibri", size = 10, face="italic"),
				legend.text = element_text(family = "Calibri", size = 9),
				legend.background = element_rect(fill="white", colour = "black"),
				legend.position = "right", 
				legend.direction = "vertical",
				legend.spacing.x = unit(0.5, 'cm'),
				legend.spacing.y = unit(0.2, 'cm'),
				panel.background = element_blank(), 
				panel.grid.major = element_line(colour = "transparent"),
				panel.grid.minor = element_line(colour = "transparent"),
				panel.border = element_blank(), 
				plot.margin = unit(c(3,3,3,3), "lines"))
		
		
		prop_FD <- ggplot() +
			geom_sf(data = rest_prop_fd_decline_sf, 
				aes(colour = get(noquote(col_var_prop)),  
					fill = get(noquote(col_var_prop))), 
				size = 0.4) +
			geom_sf(data = high_prop_fd_decline_sf, 
				colour = "red", 
				fill = "red") +
			scale_colour_gradient(high = "red", 
				low ="white", 
				name = "Proportional \ndecline", 
				na.value = "transparent", 
				limits=c(0, 0.7), 
				guide = FALSE) +
			scale_fill_gradient(high = "red", 
				low ="white", 
				name = "Proportional \ndecline", 
				na.value = "transparent", 
				limits=c(0, 0.7),
				breaks = c(0.1, 0.3, 0.5, 0.7),
				labels = c("0.1", "0.3", "0.5", "0.7+")) +
			guides(fill = guide_colorbar(ticks = FALSE)) +
			geom_sf(data = con_out, fill = "transparent", color = "grey60", size = 0.3) +
			coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) +
			theme(plot.background = element_blank(),
				legend.title = element_text(family = "Calibri", size = 10, face="italic"),
				legend.text = element_text(family = "Calibri", size = 9),
				legend.background = element_rect(fill="white", colour = "black"),
				legend.position = "right", 
				legend.direction = "vertical",
				legend.spacing.x = unit(0.5, 'cm'),
				legend.spacing.y = unit(0.2, 'cm'),
				panel.background = element_blank(), 
				panel.grid.major = element_line(colour = "transparent"),
				panel.grid.minor = element_line(colour = "transparent"),
				panel.border = element_blank(),
				plot.margin = unit(c(3,3,3,3), "lines"))
			
		
		threat <- strsplit(col_var_raw, "_")[[1]][1]
		save_nam <- paste(threat, "_prop_raw", sep = "")


		# ----------------------
		#  Draw plots
		p <- plot_grid(raw_FD, prop_FD, labels = "AUTO", ncol = 1, nrow = 2, label_size = 11)
		ggsave(paste(save_nam, ".tiff", sep = ""), width = 9, height = 9, units = 'in', dpi = 300)
		
		}



	
# =============================================================================
#  Phylogenetic Diversity
# =============================================================================
	
	# ----------------------
	#  Raw
	phylo_plot_cols <- c(45, seq(47, 65, by = 2))
	pd_res_sf <- res_sf[phylo_plot_cols]
	colnames(pd_res_sf)[1] <- "current"
	
	# ----------------------
	#  Prep plotting data for raw value decline
	raw_pd_decline_sf <- pd_res_sf %>%
		mutate_at(vars(PD.habitat2:PD.inbreeding2),  ~(current - (current * .))) %>%
		dplyr::select(-current) 

	# ----------------------
	#  Prep plotting data for proportional decline
	prop_pd_decline_sf <- pd_res_sf %>%
		mutate_at(vars(PD.habitat2:PD.inbreeding2),  ~(ifelse(current == 0, 0, 1 - .))) %>%
		dplyr::select(-current) 
	
	
	for(i in 1:10){
		
		# ----------------------
		#  Set up column that is used for raw value data
		col_var_raw <- colnames(raw_pd_decline_sf[i])[1]
		
		# ----------------------
		#  Set up column that is used for proportional data
		col_var_prop <- colnames(prop_pd_decline_sf[i])[1]
		
		# ----------------------
		#  Color scale ranges
		tmp_prop <- prop_pd_decline_sf %>%
			dplyr::select(col_var_prop)
		high_prop_pd_decline_sf <- tmp_prop %>%
			dplyr::filter_at(vars(-geometry), any_vars(. > 0.69999))
		rest_prop_pd_decline_sf <- tmp_prop %>%
			dplyr::filter_at(vars(-geometry), any_vars(. < 0.699991))
		
		
		# ----------------------
		#  Generate plots
		raw_PD <- ggplot() +
			geom_sf(data = raw_pd_decline_sf, 
				aes(colour = get(noquote(col_var_prop)),  
					fill = get(noquote(col_var_prop))), 
				size = 0.4) +
			scale_colour_gradient(high = "red", 
				low ="white", 
				name = "Raw              \ndecline", 
				na.value = "transparent", 
				limits=c(0, 10), 
				guide = FALSE) +
			scale_fill_gradient(high = "red", 
				low ="white", 
				name = "Raw              \ndecline", 
				na.value = "transparent", 
				limits=c(0, 10),
				breaks = c(2, 4, 6, 8, 10),
				labels = c("2", "4", "6", "8", "10")) +
			guides(fill = guide_colorbar(ticks = FALSE)) +
			geom_sf(data = con_out, fill = "transparent", color = "grey70", size = 0.3) +
			coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) +
			theme(plot.background = element_blank(),
				legend.title = element_text(family = "Calibri", size = 10, face="italic"),
				legend.text = element_text(family = "Calibri", size = 9),
				legend.background = element_rect(fill="white", colour = "black"),
				legend.position = "right", 
				legend.direction = "vertical",
				legend.spacing.x = unit(0.5, 'cm'),
				legend.spacing.y = unit(0.2, 'cm'),
				panel.background = element_blank(), 
				panel.grid.major = element_line(colour = "transparent"),
				panel.grid.minor = element_line(colour = "transparent"),
				panel.border = element_blank(), 
				plot.margin = unit(c(3,3,3,3), "lines"))
		
		
		prop_PD <- ggplot() +
			geom_sf(data = rest_prop_pd_decline_sf, 
				aes(colour = get(noquote(col_var_prop)),  
					fill = get(noquote(col_var_prop))), 
				size = 0.4) +
			geom_sf(data = high_prop_pd_decline_sf, 
				colour = "red", 
				fill = "red") +
			scale_colour_gradient(high = "red", 
				low ="white", 
				name = "Proportional \ndecline", 
				na.value = "transparent", 
				limits=c(0, 0.7), 
				guide = FALSE) +
			scale_fill_gradient(high = "red", 
				low ="white", 
				name = "Proportional \ndecline", 
				na.value = "transparent", 
				limits=c(0, 0.7),
				breaks = c(0.1, 0.3, 0.5, 0.7),
				labels = c("0.1", "0.3", "0.5", "0.7+")) +
			guides(fill = guide_colorbar(ticks = FALSE)) +
			geom_sf(data = con_out, fill = "transparent", color = "grey60", size = 0.3) +
			coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) +
			theme(plot.background = element_blank(),
				legend.title = element_text(family = "Calibri", size = 10, face="italic"),
				legend.text = element_text(family = "Calibri", size = 9),
				legend.background = element_rect(fill="white", colour = "black"),
				legend.position = "right", 
				legend.direction = "vertical",
				legend.spacing.x = unit(0.5, 'cm'),
				legend.spacing.y = unit(0.2, 'cm'),
				panel.background = element_blank(), 
				panel.grid.major = element_line(colour = "transparent"),
				panel.grid.minor = element_line(colour = "transparent"),
				panel.border = element_blank(),
				plot.margin = unit(c(3,3,3,3), "lines"))
			
		
		threat <- strsplit(col_var_raw, "_")[[1]][1]
		save_nam <- paste(threat, "_prop_raw", sep = "")


		# ----------------------
		#  Draw plots
		p <- plot_grid(raw_PD, prop_PD, labels = "AUTO", ncol = 1, nrow = 2, label_size = 11)
		ggsave(paste(save_nam, ".tiff", sep = ""), width = 9, height = 9, units = 'in', dpi = 300)
		
		}

	
	
# =============================================================================
#  Taxonomic Diversity
# =============================================================================
	
	# ----------------------
	#  Raw
	taxo_plot_cols <- c(3, seq(5, 23, by = 2))
	td_res_sf <- res_sf[taxo_plot_cols]
	colnames(td_res_sf)[1] <- "current"
	
	# ----------------------
	#  Prep plotting data for raw value decline
	raw_td_decline_sf <- td_res_sf %>%
		mutate_at(vars(SR.habitat2:SR.inbreeding2),  ~(current - .)) %>%
		dplyr::select(-current) 

	# ----------------------
	#  Prep plotting data for proportional decline
	prop_td_decline_sf <- td_res_sf %>%
		mutate_at(vars(SR.habitat2:SR.inbreeding2),  ~(ifelse(current == 0, 0, (current - .) / current))) %>%
		dplyr::select(-current) 

	
	for(i in 1:10){
		
		# ----------------------
		#  Set up column that is used for raw value data
		col_var_raw <- colnames(raw_td_decline_sf[i])[1]
		
		# ----------------------
		#  Set up column that is used for proportional data
		col_var_prop <- colnames(prop_td_decline_sf[i])[1]
		
		# ----------------------
		#  Color scale ranges
		tmp_prop <- prop_td_decline_sf %>%
			dplyr::select(col_var_prop)
		high_prop_td_decline_sf <- tmp_prop %>%
			dplyr::filter_at(vars(-geometry), any_vars(. > 0.69999))
		rest_prop_td_decline_sf <- tmp_prop %>%
			dplyr::filter_at(vars(-geometry), any_vars(. < 0.699991))
		
		# ----------------------
		#  Generate plots
		raw_TD <- ggplot() +
			geom_sf(data = raw_td_decline_sf, 
				aes(colour = get(noquote(col_var_prop)),  
					fill = get(noquote(col_var_prop))), 
				size = 0.4) +
			scale_colour_gradient(high = "red", 
				low ="white", 
				name = "Raw              \ndecline", 
				na.value = "transparent", 
				limits=c(0, 240), 
				guide = FALSE) +
			scale_fill_gradient(high = "red", 
				low ="white", 
				name = "Raw              \ndecline", 
				na.value = "transparent", 
				limits=c(0, 240),
				breaks = c(60, 120, 180, 240),
				labels = c("60", "120", "180", "240")) +
			guides(fill = guide_colorbar(ticks = FALSE)) +
			geom_sf(data = con_out, fill = "transparent", color = "grey70", size = 0.3) +
			coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) +
			theme(plot.background = element_blank(),
				legend.title = element_text(family = "Calibri", size = 10, face="italic"),
				legend.text = element_text(family = "Calibri", size = 9),
				legend.background = element_rect(fill="white", colour = "black"),
				legend.position = "right", 
				legend.direction = "vertical",
				legend.spacing.x = unit(0.5, 'cm'),
				legend.spacing.y = unit(0.2, 'cm'),
				panel.background = element_blank(), 
				panel.grid.major = element_line(colour = "transparent"),
				panel.grid.minor = element_line(colour = "transparent"),
				panel.border = element_blank(), 
				plot.margin = unit(c(3,3,3,3), "lines"))
		
		
		prop_TD <- ggplot() +
			geom_sf(data = rest_prop_td_decline_sf, 
				aes(colour = get(noquote(col_var_prop)),  
					fill = get(noquote(col_var_prop))), 
				size = 0.4) +
			geom_sf(data = high_prop_td_decline_sf, 
				colour = "red", 
				fill = "red") +
			scale_colour_gradient(high = "red", 
				low ="white", 
				name = "Proportional \ndecline", 
				na.value = "transparent", 
				limits=c(0, 0.7), 
				guide = FALSE) +
			scale_fill_gradient(high = "red", 
				low ="white", 
				name = "Proportional \ndecline", 
				na.value = "transparent", 
				limits=c(0, 0.7),
				breaks = c(0.1, 0.3, 0.5, 0.7),
				labels = c("0.1", "0.3", "0.5", "0.7+")) +
			guides(fill = guide_colorbar(ticks = FALSE)) +
			geom_sf(data = con_out, fill = "transparent", color = "grey60", size = 0.3) +
			coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) +
			theme(plot.background = element_blank(),
				legend.title = element_text(family = "Calibri", size = 10, face="italic"),
				legend.text = element_text(family = "Calibri", size = 9),
				legend.background = element_rect(fill="white", colour = "black"),
				legend.position = "right", 
				legend.direction = "vertical",
				legend.spacing.x = unit(0.5, 'cm'),
				legend.spacing.y = unit(0.2, 'cm'),
				panel.background = element_blank(), 
				panel.grid.major = element_line(colour = "transparent"),
				panel.grid.minor = element_line(colour = "transparent"),
				panel.border = element_blank(),
				plot.margin = unit(c(3,3,3,3), "lines"))
			
		
		threat <- strsplit(col_var_raw, "_")[[1]][1]
		save_nam <- paste(threat, "_prop_raw", sep = "")


		# ----------------------
		#  Draw plots
		p <- plot_grid(raw_TD, prop_TD, labels = "AUTO", ncol = 1, nrow = 2, label_size = 11)
		ggsave(paste(save_nam, ".tiff", sep = ""), width = 9, height = 9, units = 'in', dpi = 300)
		
		}


	
	
	
	
	##### BELOW HAS HAD NO UPDATES YET FOR COUNTRY-BASED RESULTS #######
	
# =============================================================================
#  Bias
# =============================================================================
	
	# ----------------------
	#  Pull in bias data set from "2_Analysis_graphs.R" script.
	#   Object called res_sf_bias?
	
	# ----------------------
	#  Remove points outside continents, on tiny islands and with initial species 
	#   richness of 0.	
	bias_dat_tmp <- st_intersection(res_sf_bias, con_crop) %>%
		dplyr::filter(SR0 > 0)
	bias_dat <- st_erase(bias_dat_tmp, casp_sea)
	bias_dat <- st_erase(bias_dat, black_sea)
	bias_dat <- st_erase(bias_dat, bangla)
	bias_dat <- st_erase(bias_dat, st_buffer(aral, 1))
	
	# ----------------------
	#  Generate plots
	cols <- c("-2" = "blue", "-1" = "light sky blue", "0" = "white", "1" = "pink", "2" = "red")
	
	bias_FD <- ggplot() +
		geom_sf(data = bias_dat, aes(colour = FDbias.habitat1), shape = 16, size = 0.8) +
		theme(plot.background = element_blank(),
			legend.position = "right", 
			legend.direction = "vertical",
			legend.spacing.y = unit(0.25, 'cm'),
			legend.title = element_text(size = 9),
			legend.text = element_text(size = 8),
			panel.background = element_blank(), 
			panel.grid.major = element_line(colour = "transparent"),
			panel.grid.minor = element_line(colour = "transparent"),
			panel.border = element_blank()) +
		#scale_colour_gradientn(colours = c("red", "white", "blue"), name = "Bias               ", 
			#na.value = "transparent", limits = c(-0.5, 0.5)) +
		scale_colour_manual(values = cols, name = "Bias               ",  na.value = "transparent") +
		guides(colour = guide_colorbar(ticks = FALSE)) +
		geom_sf(data = con_out, fill = "transparent", color = "grey70") +
		coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) 	
	
	# ----------------------
	#  Generate plots
	bias_PD <- ggplot() +
		geom_sf(data = bias_dat, aes(colour = PDbias.habitat1), shape = 16, size = 0.8) +
		theme(plot.background = element_blank(),
			legend.position = "right", 
			legend.direction = "vertical",
			legend.spacing.y = unit(0.25, 'cm'),
			legend.title = element_text(size = 9),
			legend.text = element_text(size = 8),
			panel.background = element_blank(), 
			panel.grid.major = element_line(colour = "transparent"),
			panel.grid.minor = element_line(colour = "transparent"),
			panel.border = element_blank()) +
		#scale_colour_gradientn(colours = c("red", "white", "blue"), name = "Bias               ", 
			#na.value = "transparent", limits = c(-0.5, 0.5)) +
		scale_colour_manual(values = cols, name = "Bias               ",  na.value = "transparent") +
		guides(colour = guide_colorbar(ticks = FALSE)) +
		geom_sf(data = con_out, fill = "transparent", color = "grey70") +
		coord_sf(crs = st_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")) 	
	
	# ----------------------
	#  Draw plots
	plot_grid(bias_FD, bias_PD, labels = "AUTO", ncol = 1, nrow = 2, label_size = 11)
