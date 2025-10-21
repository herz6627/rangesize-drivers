# Sarah Herzog
# Fall 2024


# library -----------------------------------------------------------------
library(MuMIn)
library(lme4)
library(tidyverse)

# Load files --------------------------------------------------------------
# vegetation data
dat <- read_csv("./output/1_dat-all.csv") %>% 
  mutate(across(c(year, watershed, transect, plot, species, burned, not_burned, grazing, cattle, bison, fri, family, growthcode, funcgroup, origin, census_year), as.factor))  # make sure everything is in the right format

# grasshopper data:
# http://lter.konza.ksu.edu/content/cgr02-sweep-sampling-grasshoppers-konza-prairie-lter-watersheds
grasshopper_dat2 <- read_csv("./data/CGR023.csv") %>%  # bison. CGR023:  This dataset is the age and sex of the listed species.  All 10 bags are combined for this dataset.
  rename_with(tolower) %>%  # standardize col name formats
  mutate(watershed = tolower(watershed), # make sure all watersheds are lowercase
         species = str_to_sentence(species) # fix capitalization differences
  ) %>% 
  filter(!str_detect(species, "Tettigoniidae")) # remove katydids

# Welti 2018 data
grasshopper_dat_welti <- read_csv("./output/3_range-gh-dat.csv")

# Campbell 1974 grasshopper diet data
grasshopper_dat_campbell <-
  read_csv("./data/grasshopperDiets/Campbell ghopdiets.csv") %>%
  dplyr::select(1) %>%
  rename(species = ...1) %>%
  separate_wider_delim(species, delim = " ", names = c("genus", "epithet", NA), too_few = "align_start") %>%  # split columns, removing variety names
  mutate(epithet = case_when(epithet == "sp." ~ "spp.",
                             is.na(epithet) ~ "spp.", # fill in blank epithets with spp.
                             TRUE ~ epithet)) %>%
  unite(species, genus, epithet, remove = F, na.rm = T, sep = " ") # combine genus and epithet, now cleaned
# North Dakota Agricultural Experiment Station (Fargo) & Mulkern, G. B. (1969). Food habits and preferences of grassland grasshoppers of the north central Great Plains. Agricultural Experiment Station, North Dakota State University.
grasshopper_dat_mulkern <-
  read_csv("./data/grasshopperDiets/mulkern ghop diets.csv") %>%
  dplyr::select(1) %>% # just want species names
  rename(species = ...1) %>%
  separate_wider_delim(species, delim = " ", names = c("genus", "epithet", NA), too_few = "align_start") %>%  # split columns, removing variety names
  mutate(epithet = case_when(epithet == "sp." ~ "spp.", # clean up epithets
                             TRUE ~ epithet)) %>%
  unite(species, genus, epithet, remove = F, na.rm = T, sep = " ") # combine genus and epithet, now cleaned


#functions
source("./R/4_functions.R") # big analysis loops

# Begin code --------------------------------------------------------------
# parameters
years_included = c(2003, 2020)
watersheds_included = Reduce(intersect, list(dat$watershed,grasshopper_dat2$watershed)) # get watersheds that occur in all datasets
soils_included = c("tu") # just want uplands (tully). "fl" (florence) is the other option. slope type was already removed in 1_
dats <- c("all", "gh_welti", "gh_camp", "gh_mulk") # data available to filter species by
min_cc_change <- 3 # minimum number of coverclass transistions for a species to be included in the dataset
min_dat <- 40 # minimum number of total observations for a species to be included in the dataset
model_selection = F # should model selection be run on each species' model


# Data cleaning and formatting --------------------------------------------
# trim dat based on parameters above
dat <- 
  dat %>%
  filter(year %in% years_included[1]:years_included[2], # years of interest
         watershed %in% watersheds_included,  # watersheds of interest
         soiltype %in% soils_included, # soil types of interest
         !grepl("spp.", .$species) # remove genus level identifications
  )

# get general species data
sp_info <- dat %>%
  select(species, family, growthcode, funcgroup, origin) %>%
  distinct()

# find any rows with NAs
temp <- dat[rowSums(is.na(dat)) > 0,]
print(paste(c("The following species are missing data:", unique(temp$species)))) # comes out to hybrids or genus level ID, so we are fine to drop them

dat <- dat %>%
  mutate(unq_tran = paste0(watershed, transect)) %>%
  group_by(species) %>%
  filter(
         !(n_distinct(midpoint) == 1), # remove species where all observations are 0 (due to filtering out years and watersheds above)
         n_distinct(midpoint) > min_cc_change, # drop species with fewer transitions than specified above
         n() > min_dat, # remove species with less than 20 total observations
         ) %>% 
  ungroup() %>%
  drop_na() %>%     # drop observations with missing data
  mutate(
    midpoint.1 = ifelse(midpoint != "0", 1, 0), # midpoint values to binary present/absent
    midpoint_next.1 = ifelse(midpoint_next != "0", 1, 0),
    colonization_event = ifelse(midpoint.1 == 0 & midpoint_next.1 == 1, 1, 0), # indicate extinction or colonization events
    extinction_event = ifelse(midpoint.1 == 1 & midpoint_next.1 == 0, 1, 0),
    fri = as.numeric(fri), # factoral leveling causes weird results
    n_grasshop = ((n_grasshop - mean(n_grasshop))/sd(n_grasshop)), # standardize values
    n_grasshop_previous = ((n_grasshop_previous - mean(n_grasshop_previous))/sd(n_grasshop_previous)),
    spei = ((spei - mean(spei))/sd(spei))
  )

sp <- unique(dat$species) # get available species

### initiate -----------------------------------------------------------
# using all data
big_func(dat, file_name = "all", model_selection = model_selection)

# trimming to just species eaten by grasshoppers
# campbell data
dat_gh <- dat %>% 
  filter(species %in% grasshopper_dat_campbell$species)
big_func(dat = dat_gh, file_name = "gh-campbell", model_selection = model_selection)
# mulkern data
dat_gh <- dat %>% 
  filter(species %in% grasshopper_dat_mulkern$species)
big_func(dat = dat_gh, file_name = "gh-mulkern", model_selection = model_selection)
# welti data
dat_gh <- dat %>% 
  filter(species %in% grasshopper_dat_welti$species)
big_func(dat = dat_gh, file_name = "gh-welti", model_selection = model_selection)
# all gh consumed species
dat_gh <- dat %>% 
  filter(species %in% unique(c(grasshopper_dat_welti$species, grasshopper_dat_mulkern$species, grasshopper_dat_campbell$species)))
big_func(dat = dat_gh, file_name = "gh", model_selection = model_selection)


# Adding bison:gh ---------------------------------------------------------
# using all data
big_func(dat, file_name = "all-with-bison-GH-interaction", param = c("bison + spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri + bison:n_grasshop"), model_selection = model_selection)
# grasshopper subsets
# campbell data
dat_gh <- dat %>% 
  filter(species %in% grasshopper_dat_campbell$species)
big_func(dat = dat_gh, file_name = "gh-campbell-with-bison-GH-interaction", model_selection = model_selection, param = c("bison + spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri + bison:n_grasshop"))
# mulkern data
dat_gh <- dat %>% 
  filter(species %in% grasshopper_dat_mulkern$species)
big_func(dat = dat_gh, file_name = "gh-mulkern-with-bison-GH-interaction", model_selection = model_selection, param = c("bison + spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri + bison:n_grasshop"))
# welti data
dat_gh <- dat %>% 
  filter(species %in% grasshopper_dat_welti$species)
big_func(dat = dat_gh, file_name = "gh-welti-with-bison-GH-interaction", model_selection = model_selection, param = c("bison + spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri + bison:n_grasshop"))
# all gh consumed species
dat_gh <- dat %>% 
  filter(species %in% unique(c(grasshopper_dat_welti$species, grasshopper_dat_mulkern$species, grasshopper_dat_campbell$species)))
big_func(dat = dat_gh, file_name = "gh-with-bison-GH-interaction", model_selection = model_selection, param = c("bison + spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri + bison:n_grasshop"))



# Running bison and GH datasets separately --------------------------------

# run a model with just bison (no GH) coeff
big_func(dat, file_name = "all-just-bison-coeff-no-GH", param = c("bison + spei + years_since_last_burn + fri"), model_selection = model_selection)

# grasshoppers without bison present
dat_temp <- dat %>% 
  filter(bison == 0)
big_func(dat_temp, file_name = "all-just-GH-coeff-no-bison", param = c("spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri"), model_selection = model_selection)
# grasshoppers with bison present
dat_temp <- dat %>% 
  filter(bison == 1)
big_func(dat_temp, file_name = "all-just-GH-coeff-yes-bison", param = c("spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri"), model_selection = model_selection)


# odds and ends -----------------------------------------------------------
# # number of non 0-0 transitions
# temp <- dat %>% filter(coverclass != 0 | coverclass_next != 0) 
# nrow(temp)
# # average per species
# temp2 <- temp %>% 
#   group_by(species) %>% 
#   summarise(n = n()) %>%   # number of observations per species
#   filter(n > min_dat) # these species get dropped when making the models
# temp2 %>% 
#   ungroup() %>% 
#   summarise(med = median(n))
# 
# # gh abundance
# dat %>% 
#   select(year, watershed, n_grasshop) %>% 
#   distinct() %>% 
#   summarise(med = median(n_grasshop),
#             min = min(n_grasshop),
#             max = max(n_grasshop))
# # SPEI
# dat %>% 
#   select(year, spei) %>% 
#   distinct() %>% 
#   summarise(med = median(spei),
#             min = min(spei),
#             max = max(spei))
# 
# # watersheds included
# dat %>% 
#   select(watershed) %>% 
#   distinct() %>% 
#   pull()


# END
