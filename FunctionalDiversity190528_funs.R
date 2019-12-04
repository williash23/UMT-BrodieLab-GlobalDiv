#  Functions used in Functional Diversity analysis


#  Function to prepare mammal range polygons (from IUCN shapefile) for use within 
#    FD analysis.
#   Must supply file name of mammal or bird range file 

prep_range_polys <- function(fn = NULL, ctrs = NULL){
	poly_tmp <- st_read(fn) %>%
		dplyr::filter(legend != "Extinct") %>% ## Now removes extinct species
		dplyr::select(binomial) %>%
		tidyr::separate(binomial, c("genus", "spp"), " ") %>%
		dplyr::mutate(species = paste(genus, spp, sep = "_")) %>%
		dplyr::select(species)
	poly <- st_transform(poly_tmp, st_crs(ctrs)) # Match countries sf object
	return(poly)
	}
	

#  Example function calls:
#mam <- prep_range_polys(fn = "TERRESTRIAL_MAMMALS.shp", pts = global_pts)
#spp_list <- find_int(range_polys = mam, pts = global_pts, p = i) 



#  Function to calculate FD per community. From: Owen Petchey
#   https://github.com/opetchey/ttl-resources/blob/master/functional_diversity/FD.example.1.r


Getlength <- function(xtree, comp=NA){
  if(!is.data.frame(comp))
      result <- Getlength.inner(xtree)
  if(is.data.frame(comp)){
      S <- tapply(comp[,2], comp[,1], function(x) length(x))
      FD <- tapply(comp[,2], comp[,1], function(x) Getlength.inner(list(xtree[[1]], xtree[[2]][!is.na(match(dimnames(xtree[[2]])[[1]], x)),]))) 
      FD <- FD/Getlength.inner(xtree)
      result <- data.frame(S=S, FD.new=FD)
  }
  result
}
Getlength.inner <- function(xtree){   
    if(!is.matrix(xtree[[2]]))
        result <- 0
    if(is.matrix(xtree[[2]]))
        result = sum(xtree[[1]][colSums(xtree[[2]]) != 0 & colSums(xtree[[2]]) < length(xtree[[2]][,1])])
    result
}


## 17/1/03. Written by Jens Schumacher. Please acknowledge as appropriate.

Xtree <- function(h)
    ## evaluate species branch matrix (sensu Petchey&Gaston) from a dendrogram
    ## tested for results of hclust and agnes
    ## hclust - hierarchical clustering 
    ## agnes - agglomerative clustering
    
    ## used components:
    ## merge - history of cluster merging
    ## height - actual heights at merging
    ## order - permutation to achieve nice output (needed only for agnes)
{

    species.names <- h$labels
    
    
    H1 <- matrix(0, length(h$order), 2 * length(h$order) - 2)
    l <- vector("numeric", 2 * length(h$order) - 2)
    for(i in 1:(length(h$order) - 1)) {
                                        # evaluate branch lengths
                                        #
        if(h$merge[i, 1] < 0) {
            l[2 * i - 1] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 1], 2 * i - 1] <- 1
        }
        else {
            l[2 * i - 1] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 1]]]
            H1[, 2 * i - 1] <- H1[, 2 * h$merge[i, 1] - 1] + H1[
                                                                , 2 * h$merge[i, 1]]
        }
        if(h$merge[i, 2] < 0) {
            l[2 * i] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 2], 2 * i] <- 1
        }
        else {
            l[2 * i] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 2]]]
            H1[, 2 * i] <- H1[, 2 * h$merge[i, 2] - 1] + H1[, 2 *
                                                            h$merge[i, 2]]
        }
    }
    dimnames(H1) <- list(species.names,NULL)  
    list(h2.prime=l, H1=H1)
    ## l contains the length of all the tiny branches
    ## H1: each row represents one species, each column represents one branch
    ##     1 indicates that a branch is part of the pathway from species to top of the dendrogram
    ##     0 otherwise
}