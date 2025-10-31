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
  # left_join( # add range size info
  #   select(
  #     rename(
  #       read_csv("./output/2_range-size-df.csv"), # load range size file
  #       species = sp, # match formatting with pvc dataset
  #       range_est = area 
  #       ),
  #     !range_size_class # remove extra column
  #     )
  #   ) 

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

# # bison diet data from Littleford-Colquhoun et al. 2024
# bison_dat <- 
#   read_csv("./data/rsos240136_si_004.csv")
# names(bison_dat) <- c(str_replace_all(names(bison_dat)[1:7], " ", "_"), paste(sep = "_", bison_dat[1,-(1:7)], bison_dat[2,-(1:7)], bison_dat[3,-(1:7)])) # clean up the messy column names
# bison_dat <- 
#   bison_dat[-(1:3),] %>%  # drop redundant rows
#   pivot_longer(cols = 8:last_col(), names_to = "var", values_to = "sequence_counts") %>%  # move column names to rows
#   select(Species, Assigned_by_local_or_global_reference_library, var, sequence_counts) %>% 
#   filter(!is.na(Species), # remove missing species (those identified to family or genus)
#          sequence_counts > 0,# dont care about no reads
#          str_detect(var, "Bison")) %>%   # only want bison data
#   mutate(sequence_counts = as.numeric(sequence_counts)) %>% 
#   group_by(Species, Assigned_by_local_or_global_reference_library) %>% 
#   summarise(n = sum(sequence_counts))
# # species found in bison diet using local database
# bison_dat_local <- bison_dat %>% 
#   filter(Assigned_by_local_or_global_reference_library == "Local") %>% 
#   pull(Species)
# # species found in bison diet using global database
# bison_dat_global <- bison_dat %>% 
#   filter(Assigned_by_local_or_global_reference_library == "Global") %>% 
#   pull(Species)

#functions
# source("./R/4_func-cor-test.R") # function for driver correlation tests
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
  #mutate(midpoint_change_to_next = midpoint_change_to_next+0.000001) %>%  # prevent any zeros interfering with log()
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


# # Gleditisa
# dat %>% 
#   filter(species == "Gleditsia triacanthos",
#          coverclass != 0 | coverclass_next != 0) %>% 
#   ggplot(aes(x = coverclass, y = coverclass_next, color = grazing)) +
#   geom_jitter(size = 3, alpha = 0.5)
# # heatmap
# dat %>% 
#     filter(species == "Gleditsia triacanthos",
#            coverclass != 0 | coverclass_next != 0) %>% 
#     group_by(grazing, coverclass, coverclass_next) %>% 
#     summarise(n = n()) %>% 
#   ggplot(aes(x = coverclass, y = coverclass_next, fill = n)) +
#   geom_raster() +
#   facet_grid(~grazing) +
#   paletteer::scale_fill_paletteer_c("grDevices::Purple-Yellow", direction = -1)
# # Density structured model ------------------------------------------------
# 
# 
# ## Bison -------------------------------------------------------------------
# # to get colonization and extinction in the same model as midpoint change
# cc_talley <- dat %>%
#   group_by(species) %>% 
#   filter(sum(bison == 1) >= min_dat & sum(bison == 0) >= min_dat) %>% # filter species missing observations in one of the bison treatments
#   mutate(coverclass = paste0("CC", coverclass),
#          coverclass_next = paste0("CC", coverclass_next)) %>% 
#   group_by(species, bison) %>% 
#   count(coverclass, coverclass_next) %>%
#   pivot_wider(names_from = coverclass, values_from = n,names_prefix = "from_", values_fill = 0) %>%  # pivot wider so columns = coverclass. replace missing values with 0
#   ungroup() %>% 
#   complete(species, bison, coverclass_next) %>%  # add missing cover classes to coverclass_next
#   replace(is.na(.), 0) # replace newly made NAs from 'complete' with 0 as we want 0 observations for transitions that dont exist
# if(!'from_CC7' %in% names(cc_talley)) cc_talley <- cc_talley %>% add_column(from_CC7 = 0) # add any missing columns, as would be the case if all the species do not occur in all the examined coverclasses (e.g. rare species that dont get to CC7). All CC should be accounted for in the rows should be included if they were included in the factor levels
# 
# # overwrite previous 'sp' to match current filtered data
# sp = unique(cc_talley$species)
# 
# # make matrix for each species
# trans.array <- array(dim = c(length(unique(cc_talley$coverclass_next)), length(unique(cc_talley$coverclass_next)), length(sp), 2), dimnames = list(c(0:(length(unique(cc_talley$coverclass_next))-1)), c(0:0:(length(unique(cc_talley$coverclass_next))-1)), c(sp), c("ungrazed", "bison"))) # make blank array for transition matrices to go in
# prob.array <-array(dim = c(length(unique(cc_talley$coverclass_next)), length(unique(cc_talley$coverclass_next)), length(sp), 2), dimnames = list(c(0:(length(unique(cc_talley$coverclass_next))-1)), c(0:(length(unique(cc_talley$coverclass_next))-1)), c(sp), c("ungrazed", "bison"))) # make blank array for prob of transition matrices to go in
# for (s in 1:length(sp)) { # by each species make a transition matrix
#   for (t in 0:1) { # bison treatment
#     temp <- cc_talley %>% 
#       filter(bison == t) 
#     my.mat <- temp %>%
#       filter(species == sp[s]) %>%   # filter to just 1 species
#       select(4:ncol(cc_talley)) %>%
#       as.matrix()
#     
#     # add to transition array
#     trans.array[,,sp[s], c("ungrazed", "bison")[t+1]] <- my.mat
#     
#     # divide rows by col totals to get prob. of coverclass transitioning
#     my.prob.mat <- apply(my.mat, 2, function(x){x/sum(x)})
#     my.prob.mat[is.na(my.prob.mat)] <- 0 # replace NaN from dividing cols with no transitions
#     # add to prob array
#     prob.array[,,sp[s], c("ungrazed", "bison")[t+1]] <- my.prob.mat
#     
#   }
# }
# 
# mean_CC_func <- function(x) {
#   stableCC <-  popbio::eigen.analysis(x)$stable.stage
#   out <- 
#     mean(c(
#       rep(0, round(stableCC[1]*100)),
#       rep(0.5, round(stableCC[2]*100)),
#       rep(3.5, round(stableCC[3]*100)), 
#       rep(15, round(stableCC[4]*100)), 
#       rep(37.5, round(stableCC[5]*100)), 
#       rep(62.5, round(stableCC[6]*100)), 
#       rep(85, round(stableCC[7]*100)), 
#       rep(97.5, round(stableCC[8]*100))
#     )) # mean percent cover
#   return(out)
# }
# # get stable stages
# meanCC <- as.data.frame(apply(prob.array, 3:4, function(x) mean_CC_func(x)))
# # effect size of bison on percent cover
# meanCC$effect <- log(meanCC[,"bison"]/meanCC[,"ungrazed"])
# # write
# rownames_to_column(meanCC, var = "species") %>% 
#   write_csv("./output/4_density-model-output-bison.csv")
# 
# ## Grasshopper -------------------------------------------------------------------
# # to get colonization and extinction in the same model as midpoint change
# cc_talley <- dat %>%
#   group_by(species, watershed) %>% 
#   mutate(gh_intensity = case_when(n_grasshop >= quantile(n_grasshop, 0.9) ~ "high",
#                                   n_grasshop <= quantile(n_grasshop, 0.1) ~ "low")) %>%  # calculate precentile cutoffs 
#   filter(
#     !is.na(gh_intensity) # filter to most extreme values based on percentile
#          ) %>% 
#   mutate(coverclass = paste0("CC", coverclass),
#          coverclass_next = paste0("CC", coverclass_next)) %>% 
#   group_by(species, gh_intensity) %>% 
#   count(coverclass, coverclass_next) 
# temp <- # get species with enough non-0 observations
#   cc_talley %>%
#   group_by(species, gh_intensity) %>% 
#   filter(coverclass != "CC0") %>% 
#   summarise(n = sum(n)) %>% 
#   filter(n > 50) %>%  # filter out species with fewer than this observations
#   group_by(species) %>% 
#   summarise(n = n()) %>% # count how many rows to see if we have both high and low GH levels
#   filter(n == 2) %>%  # remove species with too few GH levels
#   pull(species) # get just sp names
# 
# cc_talley <- cc_talley %>% 
#   filter(species %in% temp) %>% 
#   pivot_wider(names_from = coverclass, values_from = n,names_prefix = "from_", values_fill = 0) %>%  # pivot wider so columns = coverclass. replace missing values with 0
#   ungroup() %>% 
#   complete(species, gh_intensity, coverclass_next) %>%  # add missing cover classes to coverclass_next
#   replace(is.na(.), 0) # replace newly made NAs from 'complete' with 0 as we want 0 observations for transitions that dont exist
# if(!'from_CC7' %in% names(cc_talley)) cc_talley <- cc_talley %>% add_column(from_CC7 = 0) # add any missing columns, as would be the case if all the species do not occur in all the examined coverclasses (e.g. rare species that dont get to CC7). All CC should be accounted for in the rows should be included if they were included in the factor levels
# 
# # overwrite previous 'sp' to match current filtered data
# sp = unique(cc_talley$species)
# 
# # make matrix for each species
# trans.array <- array(dim = c(length(unique(cc_talley$coverclass_next)), length(unique(cc_talley$coverclass_next)), length(sp), 2), dimnames = list(c(0:(length(unique(cc_talley$coverclass_next))-1)), c(0:0:(length(unique(cc_talley$coverclass_next))-1)), c(sp), c("low", "high"))) # make blank array for transition matrices to go in
# prob.array <-array(dim = c(length(unique(cc_talley$coverclass_next)), length(unique(cc_talley$coverclass_next)), length(sp), 2), dimnames = list(c(0:(length(unique(cc_talley$coverclass_next))-1)), c(0:(length(unique(cc_talley$coverclass_next))-1)), c(sp), c("low", "high"))) # make blank array for prob of transition matrices to go in
# for (s in 1:length(sp)) { # by each species make a transition matrix
#   for (t in c("low", "high")) { # GH level
#     temp <- cc_talley %>% 
#       filter(gh_intensity == t) 
#     my.mat <- temp %>%
#       filter(species == sp[s]) %>%   # filter to just 1 species
#       select(4:ncol(cc_talley)) %>%
#       as.matrix()
#     
#     # add to transition array
#     trans.array[,,sp[s], t] <- my.mat
#     
#     # divide rows by col totals to get prob. of coverclass transitioning
#     my.prob.mat <- apply(my.mat, 2, function(x){x/sum(x)})
#     my.prob.mat[is.na(my.prob.mat)] <- 0 # replace NaN from dividing cols with no transitions
#     # add to prob array
#     prob.array[,,sp[s], t] <- my.prob.mat
#     
#   }
# }
# 
# # get stable stages
# meanCC <- as.data.frame(apply(prob.array, 3:4, function(x) mean_CC_func(x)))
# # effect size of bison on percent cover
# meanCC$effect <- log(meanCC[,"high"]/meanCC[,"low"])
# # write
# rownames_to_column(meanCC, var = "species") %>% 
#   write_csv("./output/4_density-model-output-gh.csv")
# 

# END
