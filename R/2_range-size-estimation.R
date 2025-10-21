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
# make a character vector of unique names and remove first name since it is just '.."
konzaSpList <- unique(konzaSp$species)[-1]

### cycle through each species and pull occurrence data.
# a list of species can be searched at once (takes a while), but since
# I want to know what species was being searched in an efficient way I am looping through the sp
spocc_df <- data.frame()
# NOTE: This takes ~2 hours to run
for (species in 1:length(konzaSpList)) {
  # find species
  spocc_temp <- occ(query = konzaSpList[species], from = c('gbif','idigbio'), has_coords = T, limit = 5000) # NOTE: this function does not make a DOI for the download, which makes it harder to cite in manuscripts
  # convert to df
  spocc_df_temp <- occ2df(spocc_temp)
  # keep only relevant columns
  spocc_df_temp <- spocc_df_temp[,colnames(spocc_df_temp) %in% c("name", "longitude", "latitude", "key", "prov", "search_sp")]
  # add species being searched
  spocc_df_temp$search_sp = konzaSpList[species]
  # add to table
  spocc_df <- rbind(spocc_df, spocc_df_temp)
  # give an update
  print(paste("Done with species", species))
}
beepr::beep()
write.csv(spocc_df, "./output/2_spocc_df.csv", row.names = FALSE)


# Occurrence data cleaning -------------------------------------------------
spocc_df <- read_csv("./output/2_spocc_df.csv")
occurrences <- spocc_df

head(occurrences)

### make sure we have all the species:
occurrences %>%
  dplyr::select(search_sp) %>%
  distinct() %>%
  count()
# so we are missing a couple (3)
# find which species we are missing
x <- occurrences %>%
  dplyr::select(search_sp) %>%
  distinct() %>%
  pull(search_sp)
#from occurrence list
setdiff(konzaSpList, x) #finds which sp are in Konza list but not occurrence list and it looks like it is just genera so we are fine
#from konza
setdiff(x, konzaSpList) #which is 0, meaning we are good to go

### How many occurrences per species are included in this file?
species_count <- occurrences %>%
  group_by(search_sp) %>%
  tally() %>% 
  arrange(n)
species_count

### get rid of NAs and species with only a few occurrences
smallSpCount <- species_count %>% 
  filter(n < 5) %>% 
  pull(search_sp)
# then remove those species and remove extra columns
occurrences_reduced <- occurrences %>% 
  filter(!search_sp %in% smallSpCount) %>% 
  rename("species" = "search_sp") %>% 
  dplyr::select(species, longitude, latitude, key)

## save csv
write.csv(occurrences_reduced, "./output/2_occurrences.csv", row.names = FALSE)

# Clean locality information  ---------------------------------------------
# Removing impossible points
occurrences_reduced <- occurrences_reduced %>%
  filter(latitude != 0, longitude != 0)

# Remove points that are botanical gardens/other biodiversity institutions (CoordinateCleaner::)
occurrences_cleaned <- cc_inst(occurrences_reduced, lon = "longitude", lat = "latitude", species = species)

# Remove duplicates
occurrences_cleaned <- occurrences_cleaned %>% 
  distinct()

# Save
write.csv(occurrences_cleaned, "./output/2_occurrences_cleaned.csv", row.names = FALSE)

# Trimming to the desired area
land <- rnaturalearth::ne_countries(returnclass = "sf", continent = c("North America", "South America")) %>%
  st_union()

# Next, convert merged dataset into a spatial file.
spList <- occurrences_cleaned %>%
  dplyr::select(species) %>%
  distinct() %>%
  pull(species)

land_points <- data.frame()
for (i in 1:length(spList)) { # (need to work with a subset so we dont use up all our storage, hence the loop) # takes ~20 min to run
  #select species
  occurrences_cleaned_temp <- occurrences_cleaned %>%
    filter(species %in% spList[i])
  # convert to spatial format
  xy_data_temp <- data.frame(x = occurrences_cleaned_temp$longitude, y = occurrences_cleaned_temp$latitude, key = occurrences_cleaned_temp$key)
  xy_data_temp <- sf::st_as_sf(xy_data_temp,
                               coords = c("x", "y"),
                               crs = st_crs(land))
  # Extract the points intersecting with land mass
  land_points_temp <- xy_data_temp %>%
    filter(st_intersects(., land, sparse = FALSE)[,1]) %>%
    mutate(species = spList[i])
  land_points <- rbind(land_points, land_points_temp)

  print(paste("done with species", i))
}
beep()

# need to remove any species that occur only a few times
land_points <- land_points %>% #
  group_by(species) %>%
  filter(n() > 5)

# did we lose any species?
setdiff(spList, land_points$species) #nope, which is good

# save as shape file
st_write(land_points, "./output/2_land_points.shp")
 
# Get convex hull range size ----------------------------------------------
if(!exists("land_points")) { # read in file if not in environment
  land_points <- st_read("./output/2_land_points.shp")
}
land_points <- land_points %>%
  st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84") %>%  # just in case
  select(-key) # =remove key column

### calculate MCP
mcp_out <- mcp(land_points, percent = 99, unin = "m", unout = "km2") # copycat code from  adehabitatHR::, just added clipping to terrestrial areas
head(mcp_out[[2]])
names(mcp_out[[1]]) <- mcp_out[[2]]$sp # name polygons for easier access

# convert to df
mcp_df <- as.data.frame(mcp_out[[2]])

 
##### get some stats on range sizes ####
areas <- mcp_out[[2]]

# Find the quantiles of range size
quant <- quantile(areas$area, probs = c(0.05, .25, .5, .75, 0.95))
# bin range sizes into the quantiles
areas <- areas %>%
  mutate(range_size_class = case_when(area < quant[1] ~ 1,
                                      between(area, quant[1], quant[2]) ~ 2,
                                      between(area, quant[2], quant[3]) ~ 3,
                                      between(area, quant[3], quant[4]) ~ 4,
                                      between(area, quant[4], quant[5]) ~ 5,
                                      area > quant[5] ~ 6))


### export
write.csv(areas, "./output/2_range-size-df.csv", row.names = FALSE)

# save occurrence info for publication
land_points <- st_read("./output/2_land_points.shp")

temp <- land_points %>%
  # clean up lat long
  mutate(geometry = st_as_text(geometry)) %>% 
  extract(geometry, into = c("lon", "lat"),
          regex = "POINT \\((-?[0-9.]+) (-?[0-9.]+)\\)",
          convert = TRUE)
temp <- temp %>% 
  mutate(
    source = case_when(
      grepl("^[0-9]+$", key) ~ "gbif",
      grepl("^[a-f0-9\\-]{36}$", key, ignore.case = TRUE) ~ "idigbio", # UUID (Universally Unique Identifier)
      TRUE ~ "unknown"
    )
  )

# save
write_csv(temp, "./output/2_occurrence-records-for-pub.csv")

####### END