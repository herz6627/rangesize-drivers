
# library -----------------------------------------------------------------
library(forecast)
library(SPEI)
library(tidyverse)

# data --------------------------------------------------------------------
# 10.6073/pasta/b768b10f9b17bafc68194a4aaa8e53c2 (link is external)
pvc <- read_csv("./data/PVC021.csv", col_types = cols(SoilType = col_character()))

# 10.6073/pasta/b1e152cc621a32c7aa623bafc016ce6c (link is external)
pbg <- read_csv("./data/PBG011.csv", col_types = cols(SoilType = col_character()))

# Konza plant list
sp_list <- read_csv("./data/PPS011.csv") %>%  # Konza plant species list
  rename(SpeCode = code, GrowthCode = "growthform", FuncGroup = "lifespan") %>% # rename to match pvc
  dplyr::rename(epithet = species) %>% 
  mutate(genus = str_to_sentence(genus)) %>% # capitalize genus
  unite(species, c(genus, epithet), sep = " ", remove = F)

# watershed burn history 
burn_dat <- 
  read_csv("./data/kfh011.csv") %>% 
  rename_with(tolower) %>% 
  mutate(burned = 1,
         not_burned = 0,
         watershed = tolower(watershed)) %>% 
  complete(watershed, year, fill = list(burned = 0)) %>% 
  complete(watershed, year, fill = list(not_burned = 1)) %>% # rows with NA data are the rows being added here
  group_by(watershed) %>% 
  group_by(grp = cumsum(not_burned == 0)) %>% # group based on when the least time it was burned. needed for next row
  mutate(years_since_last_burn = cumsum(not_burned)) %>% 
  ungroup() %>% 
  select(!grp) %>%  # remove extra column
  mutate(watershed = case_when(watershed == "1d" ~ "001d", # rename watersheds to match pvc
                               watershed == "2c" ~ "002c",
                               watershed == "2d"~ "002d",
                               watershed == "4a" ~ "004a", 
                               watershed == "4b" ~ "004b",
                               watershed == "4f" ~ "004f",
                               watershed == "fa" ~ "00fa",
                               watershed == "fb" ~ "00fb",
                               watershed == "wa" ~ "00wa",
                               watershed == "wb" ~ "00wb",
                               watershed == "20b" ~ "020b",
                               watershed == "spa" ~ "0spa",
                               watershed == "spb" ~ "0spb",
                               watershed == "sua" ~ "0sua",
                               watershed == "sub" ~ "0sub",
                               watershed == "n1a" ~ "n01a",
                               watershed == "n1b" ~ "n01b",
                               watershed == "n4a" ~ "n04a",
                               watershed == "n4d" ~ "n04d",
                               watershed == "n20a" ~ "n20a",
                               watershed == "n20b" ~ "n20b" ,
                               watershed == "r1a"  ~ "r01a" ,
                               watershed == "r1b" ~ "r01b" ,
                               watershed == "r20a" ~ "r20a" ,
                               watershed == "r20b" ~ "r20b" ,
                               watershed == "n2a" ~ "n02a" ,
                               watershed == "n2b" ~ "n02b",
                               watershed == "c1a" ~ "c01a",
                               watershed ==  "c1sb" ~ "c1sb",
                               watershed ==  "c3a" ~ "c03a",
                               watershed ==  "c3b" ~ "c03b",
                               watershed ==  "c3c" ~ "c03c",
                               watershed ==  "c3sa" ~ "c3sa",
                               watershed ==  "c3sb" ~ "c3sb",
                               watershed ==  "c3sc" ~ "c3sc",
                               TRUE ~ watershed)) %>% 
  select(watershed, year, burned, not_burned, years_since_last_burn) %>% 
  distinct() # remove duplicated years from when there were multiple burns in a watershed. I checked this F2023, and there were no issues with discrepencies in if the watershed was burned that year. 


# start --------------------------------------------------------------
# matching up useful meta-data
colnames(pbg)
colnames(pvc)
pbg <- pbg %>% 
  rename(Cover = CoverClass,
         AB_genus = Ab_genus,
         AB_species = Ab_species,
         DataCode = Datacode,
         WaterShed = Watershed,
         RecType = Rectype) %>% 
  mutate(Transect = str_to_upper(Transect),
         WaterShed = str_to_lower(WaterShed))
head(pvc)
pvc = rbind(pvc, pbg) # join both datasets

pvc <- pvc %>% 
  rename(Year = RecYear, Month = RecMonth, Day = RecDay, Watershed = WaterShed, CoverClass = Cover) %>%  # Rename columns
  filter(Transect!="E", # remove depreciated transects
         SoilType != "s" # eliminate all slope soil type
  ) %>%  
  mutate(Midpoint = case_when(CoverClass == 7 ~97.5, # create col with midpoint for each cover-class
                              CoverClass == 6 ~85,
                              CoverClass == 5 ~62.5,
                              CoverClass == 4 ~32.5,
                              CoverClass == 3 ~15,
                              CoverClass == 2 ~3.5,
                              CoverClass == 1 ~0.5))

# Replacing watershed designtations for reversal watersheds
pvc$Watershed[pvc$Watershed == "020a"] = "r01a"
pvc$Watershed[pvc$Watershed == "001a"] = "r20a"
pvc$Watershed[pvc$Watershed == "020d"] = "r01b"
pvc$Watershed[pvc$Watershed == "001c"] = "r20b"

# adding additional columns
pvc <- pvc %>% 
  mutate(
    grazing = case_when( # add grazing 
      Watershed %in% c("r20b", "001d", "004b", "020b", "002c", "002d", "004a", "004f", "r01b", "00wb", "0spa", "0spb", "00fa", "00fb", "00wa", "0sua", "0sub", "r20a", "r01a") ~ "ungrazed",
      Watershed %in% c("c01a", "c03a", "c03b", "c03c", "c3sa", "c3sb", "c3sc", "c1sb") ~ "cattle",
      Watershed %in% c("n01b", "n04d", "n20b", "n02a", "n02b", "n04a", "n20a", "n01a") ~ "bison"),
    
    fri = case_when( # add burning frequency
      Watershed %in% c("001d", "r01b", "00wb", "0spa", "0spb", "00fa", "00fb", "00wa", "0sua", "0sub", "r01a", "n01b", "n01a", "c01a", "c1sb") ~ 1,
      Watershed %in% c("002c", "002d", "n02a", "n02b") ~ 2,
      Watershed %in% c("c03a", "c03b", "c03c", "c3sa", "c3sb", "c3sc") ~ 3,
      Watershed %in% c("004b", "004a", "004f", "n04d", "n04a") ~ 4,
      Watershed %in% c("r20b", "020b", "r20a", "n20b", "n20a") ~ 20)) %>% 
  unite(unq_plot_yr, c(Watershed, SoilType, Transect, Plot, Year), remove = F) %>%  # unique plot ID
  unite(unq_plot, c(Watershed, SoilType, Transect, Plot), remove = F) %>% 
  inner_join(sp_list, by = "SpeCode") %>%  # add species meta-data
  left_join(select(burn_dat, watershed, year, burned, not_burned, years_since_last_burn), by = c("Watershed" = "watershed", "Year" = "year")) %>%  # add burn history
  group_by(unq_plot, Year, species) %>% 
  dplyr::slice_max(order_by = CoverClass, n = 1, with_ties = F) %>%  # !!! important: if there are multiple observations in a year, pick the higher CC value
  ungroup()
### Rename transect E to D and delete old transect D
# First, eliminate the old D transect from 1983 to 1987
pvc <- pvc[!(pvc$Watershed == "n20b" & pvc$SoilType == "f" & pvc$Transect == "D" & pvc$Year < 1987),]
# Rename transect E to D for years before 1987
pvc <- pvc %>% mutate(Transect = recode(Transect,'E' = "D"))

# get some basic info for methods -----------------------------------------
# number of observations
pvc %>% 
  group_by(species) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 

# cover
pvc %>%
  mutate(midpoint = case_when(CoverClass == 0 ~ 0, # convert to midpoint from coverclass
                              CoverClass == 1 ~ 0.005,
                              CoverClass == 2 ~ 0.035,
                              CoverClass == 3 ~ 0.15,
                              CoverClass == 4 ~ 0.375,
                              CoverClass == 5 ~ 0.625,
                              CoverClass == 6 ~ 0.85,
                              CoverClass == 7 ~ 0.975)) %>% 
  group_by(species) %>% 
  summarise(ave_cover = mean(midpoint)) %>% 
  arrange(desc(ave_cover)) 

# grass v forbs
pvc %>% 
  group_by(GrowthCode) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 


# remove unwanted treatments/issues  --------------------------------------
pvc <- pvc %>% 
  dplyr::filter(origin == "n") %>%  # remove non-natives (keep only natives)
  dplyr::filter(!str_detect(Watershed, '^c')) # ***** remove cattle watersheds

# Add 0 CC and reformat  ----------------------------------------------------------------

pvc_long <- pvc %>% 
  # add blank rows for each species for each plot
  select(unq_plot, Year, Watershed, SoilType, Transect, Plot, species, CoverClass) %>%  # trim down data
  complete(Year, Watershed, SoilType, Transect, Plot, species, fill = list(CoverClass = 0, Midpoint = 0)) %>%  # fill new observation data with 0
  unite(unq_plot, c(Watershed, SoilType, Transect, Plot), remove = F) %>%   # add back unique plot ID
  group_by(Watershed, species) %>%  # remove observations when plant has not been observed in watershed
  filter(any(CoverClass != 0)) %>%
  mutate(midpoint = case_when(CoverClass == 0 ~ 0, # convert to midpoint from coverclass
                              CoverClass == 1 ~ 0.005,
                              CoverClass == 2 ~ 0.035,
                              CoverClass == 3 ~ 0.15,
                              CoverClass == 4 ~ 0.375,
                              CoverClass == 5 ~ 0.625,
                              CoverClass == 6 ~ 0.85,
                              CoverClass == 7 ~ 0.975),  # convert CC to midpoint, as scale of change is larger (and biologically relevant) with percentage cover
         CoverClass_next = lead(CoverClass, order_by = Year), # need to add next years data (t+1) 
         #CoverClass_previous = lag(CoverClass, order_by = Year),
         midpoint_next = lead(midpoint, order_by = Year),
         midpoint_previous = lag(midpoint, order_by = Year)
         ) %>% 
  # add back meta data
  left_join(select(burn_dat, watershed, year, burned, not_burned, years_since_last_burn), by = c("Watershed" = "watershed", "Year" = "year")) %>%  # add burn history
  mutate(
    grazing = case_when( # add grazing 
      Watershed %in% c("r20b", "001d", "004b", "020b", "002c", "002d", "004a", "004f", "r01b", "00wb", "0spa", "0spb", "00fa", "00fb", "00wa", "0sua", "0sub", "r20a", "r01a") ~ "ungrazed",
      Watershed %in% c("c01a", "c03a", "c03b", "c03c", "c3sa", "c3sb", "c3sc", "c1sb") ~ "cattle",
      Watershed %in% c("n01b", "n04d", "n20b", "n02a", "n02b", "n04a", "n20a", "n01a") ~ "bison"),
    bison = case_when(Watershed %in% c("n01b", "n04d", "n20b", "n02a", "n02b", "n04a", "n20a", "n01a") ~ 1, # have a bison specific column
                      .default = 0),
    cattle = case_when(Watershed %in% c("c01a", "c03a", "c03b", "c03c", "c3sa", "c3sb", "c3sc", "c1sb") ~ 1, # have a cattle specific column (which has been filtered out above)
                      .default = 0),
    fri = case_when( # add burning frequency
      Watershed %in% c("001d", "r01b", "00wb", "0spa", "0spb", "00fa", "00fb", "00wa", "0sua", "0sub", "r01a", "n01b", "n01a", "c01a", "c1sb") ~ 1,
      Watershed %in% c("002c", "002d", "n02a", "n02b") ~ 2,
      Watershed %in% c("c03a", "c03b", "c03c", "c3sa", "c3sb", "c3sc") ~ 3,
      Watershed %in% c("004b", "004a", "004f", "n04d", "n04a") ~ 4,
      Watershed %in% c("r20b", "020b", "r20a", "n20b", "n20a") ~ 20))  %>% 
  left_join(select(sp_list, species, family, GrowthCode, FuncGroup, origin)) %>% # add species meta-data
  rename_with(tolower) %>% 
  mutate(soiltype = case_when(soiltype == "t" ~ "tu", # fix r reading soil types as logicals
                       soiltype == "f" ~ "fl"),
         census_year = paste0(year-1, "_", year) # indicate census year 
         )
beepr::beep() # convenient beep
write_csv(pvc_long, "./output/1_PVC-long.csv") 




 
# Add climate data --------------------------------------------------------
pvc_long <- read_csv("./output/1_PVC-long.csv")

# Konza climate data
# http://lter.konza.ksu.edu/content/awe01-meteorological-data-konza-prairie-headquarters-weather-station
clim_dat_raw <- read_csv("./data/AWE012.csv", col_types = cols(.default = "?",
                                                               TMAX = "n",
                                                               TMIN = "n",
                                                               TAVE = "n",
                                                               DHUMID = "n",
                                                               DSRAD = "n",
                                                               DPPT = "n",
                                                               SMAX = "n",
                                                               SMIN = "n",
                                                               S_AVE = "n",
                                                               WAVE = "n"), na=c("",".", "NA")) %>%  #daily
  rename_with(tolower) # standardize col name formats

# Konza HQ precipitation. More accurate than the AWE01 tipping bucket method (per the AWE01 metadata)
# http://lter.konza.ksu.edu/content/awe01-meteorological-data-konza-prairie-headquarters-weather-station
precip_dat_raw <- read_csv("./data/APT011.csv", col_types = c(.default = "?",
                                                              ppt = "n"), na=c("",".", "NA")) # missing values = ".". These should be replaced as NA when converting to numeric

# average precip and temp data - for methods
print(paste0("Mean annual temp for konza from 1983 to 2023 is ",
             (clim_dat_raw %>% 
              filter(recyear > 1983 & recyear < 2023) %>% 
              group_by(recyear) %>% 
              summarise(at = mean(tave, na.rm = T)) %>% # ave daily precip totals to get mean annual temp
              summarise(mean_t = mean(at, na.rm = T)) %>% # mean anual temp across years
                pull() %>% 
                round(., 2)), 
              " C"))  

print(paste0("Mean annual precip for konza from 1983 to 2023 is ",
             (precip_dat_raw %>% 
               mutate(recyear = year(mdy(RecDate))) %>% 
               filter(recyear > 1983 & recyear < 2023) %>% 
               group_by(recyear) %>% 
               summarise(appt = sum(ppt, na.rm = T)) %>% 
               summarise(mean_ppt = mean(appt, na.rm = T)) %>% 
               pull() %>% 
               round(., 2)),
             " mm"))  

# average climate data over census year
clim_dat <- clim_dat_raw %>% 
  mutate(census_year = case_when(recmonth < 6 ~ paste0(recyear-1, "_", recyear),
                                 recmonth >= 6 ~ paste0(recyear, "_", recyear+1))) %>% 
  group_by(census_year) %>%
  summarise(temp = mean(tave, na.rm = T),
            precip = mean(dppt, na.rm = T))

# get drought index
# average climate data by month
clim_by_month <- clim_dat_raw %>% 
  group_by(recyear, recmonth) %>% 
  summarise(tmean_mean = mean(tave, na.rm = T)) # temp means
clim_by_month <- 
  precip_dat_raw %>%  # precip means
  mutate(RecDate = lubridate::mdy(RecDate), # convert to date format
         recyear = year(RecDate),
         recmonth = month(RecDate)) %>% 
  group_by(recyear, recmonth) %>% 
  summarise(pr_mean = mean(ppt, na.rm = T)) %>% 
  ungroup() %>% 
  add_row(recyear = c(1986,1989, 2009, 2016), recmonth = c(1, 11, 1, 1), pr_mean = 0) %>%   # adding some dummy data for months when precip was too low to register
  full_join(clim_by_month) %>% # add precip and temp dat together
  filter(recyear %in% 1983:2022) %>%  # trim data
  arrange(recyear) #cosmetic

# calculate potential evapotranspiration (PET)
clim_by_month$pet <- SPEI::thornthwaite(Tave = clim_by_month$tmean_mean, lat = 39.09306) # konza latitude
# calculate climatic water balance (BAL)
clim_by_month$bal <- clim_by_month$pr_mean - clim_by_month$pet
# convert to time series (ts) for convenience
clim_by_month <- ts(clim_by_month[, -c(1, 2)], end = c(2022, 12), frequency = 12) #end = c(2020, 12)
plot(clim_by_month)

# One and twelve-months SPEI
spei1 <- SPEI::spei(clim_by_month[, "bal"], 1, na.rm = T)
spei12 <- SPEI::spei(clim_by_month[, "bal"], 12, na.rm = T)

# join
spei_res <- full_join(data.frame(spei1.z=as.matrix(spei1$fitted), date=zoo::as.Date(time(spei1$fitted))),
                      data.frame(spei12.z=as.matrix(spei12$fitted), date=zoo::as.Date(time(spei12$fitted))))
# # SPEI values between -1 and 1 are considered near normal for a given area, whereas values below -1 signify drought and values above 1 signify unusually moist conditions. Because drought conditions fluctuate naturally, it is helpful to look at average conditions over several years to explore how drought is connected to long-term climate change.5

# average by year 
spei.annual <- 
  spei_res %>%
  mutate(recyear = year(date), # fix col types
         recmonth = month(date),
         census_year = case_when(recmonth < 6 ~ paste0(recyear-1, "_", recyear), # get "census year" - AKA time between plant monitoring surveys
                                 recmonth >= 6 ~ paste0(recyear, "_", recyear+1))) %>% 
  group_by(census_year) %>%
  summarise(spei1.mean = mean(spei1.z), #average
            spei1.sd = sd(spei1.z), # sd
            spei12.mean = mean(spei12.z),
            spei12.sd = sd(spei12.z))


# Add grasshopper data ----------------------------------------------------
# Konza grasshopper data: 
# http://lter.konza.ksu.edu/content/cgr02-sweep-sampling-grasshoppers-konza-prairie-lter-watersheds
grasshopper_dat2 <- read_csv("./data/CGR023.csv") %>%  # bison. CGR023:  This dataset is the age and sex of the listed species.  All 10 bags are combined for this dataset.
  rename_with(tolower) %>%  # standardize col name formats
  mutate(watershed = tolower(watershed), # make sure all watersheds are lowercase
         species = str_to_sentence(species), # fix capitalization differences
         species = case_when( # species names are a hot mess. NOTE: dont use 'specode' to identify grasshopper species, as there are many errors.
           species == "Ageneotett deorum" ~ "Ageneotettix deorum",
           species ==  "Arphia species"  ~  "Arphia spp.",
           species == "Brachystol magna" ~ "Brachystola magna",
           species == "Schistocer lineata" ~ "Schistocerca lineata",
           species == "Paratylotr brunneri" ~ "Paratylotropidia brunneri",
           species == "Paratylota brunneri" ~ "Paratylotropidia brunneri", 
           species == "Campylacan olivacea" ~ "Campylacantha olivacea",
           species == "Hesperotet speciosus" ~ "Hesperotettix speciosus",
           species == "Hesperotet viridis" ~ "Hesperotettix viridis",
           species == "Hesperotet spp." ~ "Hesperotettix spp.",
           species == "Hesperotet species" ~ "Hesperotettix spp.",
           species == "Phoetaliot nebrascen" ~ "Phoetaliotes nebrascensis",
           species == "Melanoplus sanguinip" ~ "Melanoplus sanguinipes",
           species == "Melanoplus femurrubr" ~ "Melanoplus femurrubrum",
           species == "Melanoplus different" ~ "Melanoplus differentialis",
           species == "Melanoplus bivittatu" ~ "Melanoplus bivittatus",
           species == "Melanoplus species" ~ "Melanoplus spp.",
           species == "Syrbula admirabil" ~ "Syrbula admirabilis",
           species == "Mermiria bivitatta" ~ "Mermiria bivittata",
           species == "Pseudopoma brachypte" ~ "Pseuodopomala brachyptera",
           species == "Boopedon auriventr" ~ "Boopedon auriventris",
           species == "Mermiria species" ~ "Mermiria spp.",
           species == "Chortophag viridifas" ~ "Chortophaga viridifasciata",
           species == "Arphia xanthopte" ~ "Arphia xanthoptera",
           species == "Hadrotetti trifascia" ~ "Hadrotettix trifasciatus",
           species == "Hippiscus rugosus" ~ "Hippiscus ocelote",
           species == "Pardalopho haldemani" ~ "Pardalophora haldemani",
           species == "Arphia species" ~ "Arphia spp.",
           species == "Schistocer obscura" ~ "Schistocerca obscura",
           species == "Encoptolop sordidus" ~ "Encoptolophus sordidus",
           species == "Melanoplus angustipe" ~ "Melanoplus angustipennis",
           species == "Chortophag viridifas" ~ "Chortophaga viridifasciata",
           species == "Xanthippus corallipe" ~ "Xanthippus corallipes",
           species == "Encoptolop subgracilis" ~ "Encoptolophus subgracilis",
           species == "Pardalopho spp." ~ "Pardalophora spp.",
           species == "Pardalphor spp." ~ "Pardalophora spp.",
           species == "Encoptolp spp." ~ "Encoptolphus spp.",
           species == "Oedipodinae" ~ "Oedipodinae spp.",
           TRUE ~ species)   # fix spelling errors
  ) 

# remove non-acridid species, as these were accounted for starting in 2014
non_acridid_genera <- c("Oecanthus", "Tettigoniidae", "Gryllidae",
                        "Scudderia", "Orchelimum", "Neoconocephalus",
                        "Arethaea", "Amblycorypha", "Pediodectes", "Tettigidea")
grasshopper_dat2 <- grasshopper_dat2 %>% 
  filter(!str_detect(species, paste0("^(", paste(non_acridid_genera, collapse = "|"), ")")))


# add up grasshoppers over species and year
grasshopper_dat <- 
  grasshopper_dat2 %>% 
  group_by(recyear, watershed) %>% 
  summarise(n_grasshop = sum(total, na.rm = T)/n_distinct(repsite)) %>%  # divide number of grasshoppers by sampling effort in that watershed. Only a couple years was there uneven sampling in a watershed. 
  group_by(watershed) %>% 
  mutate(n_grasshop_previous = lag(n_grasshop, order_by = recyear)) %>% 
  rename(year = recyear) # match with pvc

# merge datasets ----------------------------------------------------------
dat_all <- pvc_long %>% 
  left_join(clim_dat) %>% 
  left_join(select(
    rename(spei.annual, spei = spei12.mean), census_year, spei)) %>% 
  left_join(grasshopper_dat) %>% 
  mutate(across(c(year, watershed, transect, plot, species, midpoint, midpoint_previous, burned, not_burned, grazing, cattle, bison, fri, family, growthcode, funcgroup, origin, census_year), as.factor)) # make sure everything is in the right format

write_csv(dat_all, "./output/1_dat-all.csv")
# END

