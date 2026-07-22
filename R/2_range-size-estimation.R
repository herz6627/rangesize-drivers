# Sarah Herzog
# Fall 2024


# library -----------------------------------------------------------------
library(CoordinateCleaner)
library(beepr) # convenient beep when done running
library(raster)
library(ridigbio)
library(spocc)
library(taxize)
library(sf)
library(tidyverse)

# Load files --------------------------------------------------------------
konzaSp <- read_csv("./data/PPS011.csv") # konza species list
head(konzaSp) # (first few rows are pretty empty)

source("./R/func_mcp-sh.R") # modified version of adehabitatHR::mcp

# Begin code --------------------------------------------------------------
# rename columns to not include spaces
names(konzaSp)<-gsub("\\s","_",names(konzaSp))

# Add column of species name joined together
konzaSp <- konzaSp %>%
  rename(epithet = species) %>%
  mutate(species = paste(genus, epithet),
         species = str_to_sentence(species))
# make a character vector of unique names and remove all species names that are not actually species names
konzaSpList <- unique(konzaSp$species)
konzaSpList <- konzaSpList[!konzaSpList %in% 
  c(". .", "Annual forb", "Cyperus spp.", "Aster spp.", "Carex spp.","Allium spp.","Erigeron spp.", " Euphorbia spp.")]

### cycle through each species and pull occurrence data.
# a list of species can be searched at once (takes a while), but since
# I want to know what species was being searched in an efficient way I am looping through the sp
downloadkeys <- as.character(rep(NA, length(konzaSpList)))
citations <- rep(NA, length(konzaSpList))
# NOTE: This takes about 2 weeks to run 
for (i in 1:length(konzaSpList)) {
 
  # find species
  taxonKey_i <- rgbif::name_backbone(konzaSpList[i])  # get taxon key for species 
  species_search <- rgbif::occ_download_wait( # wait for an occurrence download to be done, because they must be done in sequence -- this command gets rid of the warnings about wait time that Sarah was dealing with
    rgbif::occ_download(# this command gets unlimited records AND generates a citation
      # and, all of these downloads will show up on AML's GBIF account online. NB
      rgbif::pred("hasGeospatialIssue", FALSE),user= "amlouthan", pwd= "mYwhu3-wukfaz-qahpaf", email= "amlouthan@ksu.edu", 
      rgbif::pred("hasCoordinate", TRUE), 
      rgbif::pred("occurrenceStatus","PRESENT"),
      rgbif::pred_not(rgbif::pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))), 
      rgbif::pred("taxonKey", taxonKey_i$usageKey)  ,
      rgbif::pred_or(rgbif::pred_lt("coordinateUncertaintyInMeters",1000),rgbif::pred_isnull(
        "coordinateUncertaintyInMeters")), #coordinateUncertaintyInMeters is less 1000 meter or is left blank.
      format = "SIMPLE_CSV"
    )
  ) 
  downloadkeys[i] <- species_search$key
  citations[i] <-  rgbif::gbif_citation(species_search$key)$download
  write.csv(downloadkeys, file= "./output/2_GBIF_download_keys.csv") # putting these write commands inside the loop because this loop fails all the time
  write.csv(citations, file= "./output/2_GBIF_citations_for_pub.csv")
}
write.csv(konzaSpList, file= "./output/2_GBIF_species_names.csv")


# This section is much faster---
# make vector of packages needed
packages_needed <- c("rgbif", "dplyr", "CoordinateCleaner", "sp", "tidyverse", "tidyr", 
                     "raster", "ggplot2", "maps","stringr", "adehabitatHR")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}



# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}


# load data---- 
downloadkeys <-   read.csv("./output/2_GBIF_download_keys.csv")[,"x"] # putting these write commands inside the loop because this loop fails all the time

sn_cap <- read.csv("./output/2_GBIF_species_names.csv")[,2]


# cleaning occurence data & calculating convex hull area---- 

good_species <- rep(NA, length(sn_cap))
# Trimming to the desired area
land <- rnaturalearth::ne_countries(returnclass = "sf", continent = c("North America", "South America")) %>%
  st_union()
all_areas <- as.data.frame(matrix(NA, nrow=length(sn_cap), ncol= 2))
for (i in 1:length(sn_cap)){
    # code that runs in both parallel & not 
    good_species_i <- TRUE
    
    # download GBIF data from the GBIF website
    d <- rgbif::occ_download_get(key = downloadkeys[i], path= ".", overwrite=FALSE) %>%
      rgbif::occ_download_import()
    file.remove(paste(downloadkeys[i], ".zip", sep= "")) # then remove from working directory; files are huge
    
    # remove rows with occurrence issues 
    # reference: https://data-blog.gbif.org/post/gbif-filtering-guide/
    
    d <- d %>%
      filter(coordinatePrecision < 0.01 | is.na(coordinatePrecision)) %>% 
      filter(!coordinateUncertaintyInMeters %in% c(301,3036,999,9999)) %>% # remove any records taht were assigned the default value for coordinate uncertainty-- often erroneous
      filter(!decimalLatitude == 0 | !decimalLongitude == 0) %>%
      CoordinateCleaner::cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
      CoordinateCleaner::cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
      CoordinateCleaner::cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
      CoordinateCleaner::cc_sea() %>% # remove from ocean 
      distinct(decimalLongitude,decimalLatitude,speciesKey,datasetKey, .keep_all = TRUE) # removes potentially duplicated recrods
    
    if (dim(d)[1]< 5 ) {good_species_i <- FALSE} else {
      # select relevant columns 
      d <- d[ , c("species","scientificName", "decimalLongitude", 
                  "decimalLatitude")]
      
      xy_data_temp <- data.frame(x = d$decimalLongitude, y = d$decimalLatitude, species = sn_cap[i])
      xy_data_temp <- sf::st_as_sf(xy_data_temp,
                                   coords = c("x", "y"),
                                   crs = st_crs(land))
      # Extract the points intersecting with land mass
      land_points_i <- xy_data_temp %>%
        filter(st_intersects(., land, sparse = FALSE)[,1])
      if (dim(d)[1]< 5 ) {good_species_i <- FALSE} else {
        output_file <- file.path("./output", paste0("2_",gsub(" ", "_", sn_cap[i]), ".csv"))
        # write.csv(d, file = output_file, row.names = FALSE) # you do not actually need to write the coordinates to the file, and this uses memory and time
        
        land_points_i <- land_points_i %>%
          st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")  # just in case

        all_areas[i,] <- mcp(as(land_points_i, "Spatial"), percent = 99, unin = "m", unout = "km2")@data
      }

  good_species[i] <- good_species_i
}

}
write.csv(good_species, file= "./output/2_good_species.csv")
write.csv(all_areas, file= "./output/2_range-size-df.csv")


# Find the quantiles of range size
quant <- quantile(all_areas$area, probs = c(0.05, .25, .5, .75, 0.95))
# bin range sizes into the quantiles
areas <- areas %>%
  mutate(range_size_class = case_when(area < quant[1] ~ 1,
                                      between(area, quant[1], quant[2]) ~ 2,
                                      between(area, quant[2], quant[3]) ~ 3,
                                      between(area, quant[3], quant[4]) ~ 4,
                                      between(area, quant[4], quant[5]) ~ 5,
                                      area > quant[5] ~ 6))


####### END