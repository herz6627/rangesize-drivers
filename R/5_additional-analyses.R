# Sarah Herzog
# Fall 2024


# library -----------------------------------------------------------------
library(MuMIn)
library(lme4)
library(tidyverse)

# Load files --------------------------------------------------------------
# range size information
ranges <- read_csv("./output/2_range-size-df.csv") %>% 
  select(!range_size_class) %>% 
  rename(species = sp,
         range_size = area) %>% 
  mutate(range_size = range_size/10^7)
# plant info
# http://lter.konza.ksu.edu/content/pps01-konza-prairie-plant-species-list
sp_info <- read_csv("./data/PPS011.csv") %>% 
  rename(epithet = species) %>% 
  unite(species, genus, epithet, remove = F, sep = " ") %>% 
  mutate(species = str_to_sentence(species))



#  functions --------------------------------------------------------------
source("./R/5_functions.R")


model_selection = F # needs to match same variable in 4_*
include_interaction = T # do we want to include results with bison:GH interaction effects?
include_zeros = F
transform_coeffs = T # if FALSE: still takes absolute values, but no log()
filterXL = T # should GH analyses/coeffs contain plant species with XL range sizes. if T, removes largest (>2) range sizes
var_names = c('(Intercept)', "bison1", "fri","n_grasshop", "spei", "years_since_last_burn", 'bison1:n_grasshop', "bison1_mod", "n_grasshop_mod")


# set theme for figures -------------------------------------------------------------------
theme_pretty <- theme_set(
  theme_light()
)
theme_pretty <- theme_update(
  strip.background = element_rect(
    color="black", fill="gray80", size=1.5, linetype="blank"),
  strip.text = element_text(colour = 'black')
)


# load and format files --------------------------------------------------------------
if (include_interaction == T) {
  # interaction effect ------------------------------------------------------
  # colonization models
  col_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # midpoint change models
  cc_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # extinction models
  ext_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
}

if (include_interaction == F) {
  # colonization models
  col_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # midpoint change models
  cc_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # extinction models
  ext_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
}

# colonization
col_dat <- mod_list_to_df(col_mods) %>% 
  left_join(ranges) %>%  # add range size for each species
  select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
  {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
  mutate(
    temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
    bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
    n_grasshop_mod = n_grasshop+temp) %>% 
  select(-temp) %>% 
  # if zeros are included, we need to make them non-zero before loging
  {if (include_zeros == T & transform_coeffs == T) mutate(., 
                                                          across(any_of(var_names), abs), 
                                                          across(all_of(var_names), ~.+0.001), 
                                                          across(all_of(var_names), log)) else if(include_zeros == F & transform_coeffs == T) mutate(., 
                                                                                                                                                     across(all_of(var_names), abs), 
                                                                                                                                                     across(all_of(var_names), log)) else mutate(., # if no transformation, just absolute values
                                                                                                                                                                                                 across(all_of(var_names), abs)) } # no zeros for log. response variables are now transformed to log(abs())

# cc change
cc_dat <- mod_list_to_df(cc_mods) %>% 
  left_join(ranges) %>%  # add range size for each species
  select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
  {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
  mutate(
    temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
    bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
    n_grasshop_mod = n_grasshop+temp) %>% 
  select(-temp) %>% 
  # if zeros are included, we need to make them non-zero before loging
  {if (include_zeros == T & transform_coeffs == T) mutate(., 
                                                          across(any_of(var_names), abs), 
                                                          across(all_of(var_names), ~.+0.001), 
                                                          across(all_of(var_names), log)) else if(include_zeros == F & transform_coeffs == T) mutate(., 
                                                                                                                                                     across(all_of(var_names), abs), 
                                                                                                                                                     across(all_of(var_names), log)) else mutate(., # if no transformation, just absolute values
                                                                                                                                                                                                 across(all_of(var_names), abs)) } # no zeros for log. response variables are now transformed to log(abs())
# extinction
ext_dat <- mod_list_to_df(ext_mods) %>% 
  left_join(ranges) %>%  # add range size for each species
  select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
  {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
  mutate(
    temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
    bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
    n_grasshop_mod = n_grasshop+temp) %>% 
  select(-temp) %>% 
  # if zeros are included, we need to make them non-zero before loging
  {if (include_zeros == T & transform_coeffs == T) mutate(., 
                                                          across(any_of(var_names), abs), 
                                                          across(all_of(var_names), ~.+0.001), 
                                                          across(all_of(var_names), log)) else if(include_zeros == F & transform_coeffs == T) mutate(., 
                                                                                                                                                     across(all_of(var_names), abs), 
                                                                                                                                                     across(all_of(var_names), log)) else mutate(., # if no transformation, just absolute values
                                                                                                                                                                                                 across(all_of(var_names), abs)) } # no zeros for log. response variables are now transformed to log(abs())
# combine data and format -------------------------------------------------
# !@!!! MAKE SURE THIS IS CORRECT WHEN CHANGING THE DATASETS
dat <- col_dat %>% 
  mutate(mod = "colonization") %>% 
  full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
  full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
  left_join(sp_info) %>% # add species info
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%  # organize columns
  pivot_longer(cols = all_of(c('(Intercept)', var_names)), names_to = "var", values_to = "val") %>% 
  mutate(var = as.factor(var),
         # val = val,
         # type = case_when(
         #   var %in% c("bison1", "log(n_grasshop)", "log(n_grasshop_previous)") ~ "biotic",
         #   var %in% c("spei", 'fri2', "fri4", 'fri20', "years_since_last_burn", "log(years_since_last_burn + 1)", "fri") ~ "abiotic")
  ) %>% 
  filter(var != "(Intercept)") 



# stats for manuscript ----------------------------------------------------

# # how many coefficients are used for each model
# dat %>% 
#   group_by(mod, var) %>% 
#   summarise(n = sum(!is.na(val))) %>% 
#   pivot_wider(names_from = var, values_from = n) %>% 
#   write_csv("./output/5_summary-table-n-coeff.csv")
# 
# # range sizes
# dat %>% 
#   summarise(
#     med = median(range_size),
#     sd = sd(range_size),
#     min = min(range_size),
#     max = max(range_size)
#   )
# 
# # effect sizes
# dat %>%
#   # mutate(val = log(abs(val))) %>% 
#   group_by(var, mod) %>%
#   summarise(
#     med = median(val, na.rm = T),
#     mean = mean(val, na.rm = T),
#     sd = sd(val, na.rm = T),
#     min = min(val, na.rm = T),
#     max = max(val, na.rm = T)
#   ) %>%
#   arrange(desc(mean)) %>% 
#   View()
# 
# # primary drivers
# dat %>%
#   filter(var %in% c("bison1", "fri", "n_grasshop", "n_grasshop_previous", "years_since_last_burn", "spei"),
#          range_size < 2) %>%
#   group_by(var, mod) %>%
#   summarise(
#     med = median(val, na.rm = T),
#     mean = mean(val, na.rm = T),
#     sd = sd(val, na.rm = T),
#     min = min(val, na.rm = T),
#     max = max(val, na.rm = T)
#   ) %>%
#   arrange(desc(mean)) %>%
#   View()



# coefficient values graph ------------------------------------------------

# # all
# plot_coefficients(dat, filterXL = TRUE, coeff_type = "all")
# ggsave("./figs/5_coeff-values-log.png", height = 3, width = 6)

# primary coeffs
plot_coefficients(dat, filterXL = TRUE, coeff_type = "primary")
ggsave("./figs/5_coeff-values-log-primary.png", height = 3, width = 4)

# secondary coeffs
plot_coefficients(dat, filterXL = TRUE, coeff_type = "secondary")
ggsave("./figs/5_coeff-values-log-secondary.png", height = 3, width = 6)


# models and model testing using absolute values of coefficients and log transformation ----------------------------------
pdf("./figs/5_col-resids.pdf")
col_mod_test_abs <- col_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T)
dev.off()
pdf("./figs/5_cc-resids.pdf")
cc_mod_test_abs <- cc_dat %>%
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%
  # filter(species != "Gleditsia triacanthos") %>%
  mod_test(include_weights = T)
dev.off()
pdf("./figs/5_ext-resids.pdf")
ext_mod_test_abs <- ext_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T)
dev.off()
# which models have range size as a predictor?
res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
res

# plot
# get linear model predictions
pred_col <- bind_rows(col_mod_test_abs$predicted) %>% 
  mutate(mod = "colonization")
pred_cc <- bind_rows(cc_mod_test_abs$predicted) %>% 
  mutate(mod = "midpoint_change")
pred_ext <- bind_rows(ext_mod_test_abs$predicted) %>% 
  mutate(mod = "extinction")
pred_dat <- bind_rows(pred_col, pred_cc, pred_ext)


make_fig(dat_df = dat, pred_df = pred_dat, res, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop"), height = 5.5, width = 4)
# # plot with only drivers of interest
# make_fig(dat_df = dat, pred_df = pred_dat, sig_coeffs_list = res, height = 4, width = 4, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop"), fld = "./figs/5_mod-res-subset-primary.png")
# # plot with non-drivers of interest
# make_fig(dat_df = dat, pred_df = pred_dat, sig_coeffs_list = res, height = 6.5, width = 4, filter.var.vec = c("fri", "bison1", "spei", "n_grasshop", "years_since_last_burn",  "n_grasshop_previous"), fld = "./figs/5_mod-res-subset-secondary.png")

# save
res %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  rename("range.coeff" = 2) %>% 
  mutate(mod = case_when(mod == "cc.bison1.n_grasshop" ~ "cc.bison1:n_grasshop", # change character for interactions
                         mod == "col.bison1.n_grasshop" ~ "col.bison1:n_grasshop",
                         mod == "ext.bison1.n_grasshop" ~ "ext.bison1:n_grasshop",
                         .default = mod)) %>% 
  separate_wider_delim(1, delim = ".", names = c("mod", "driver")) %>% 
  write_csv("./output/5_summary-table-range-coeff-vals.csv")


# removing outliers (extra large range size) --------------------------------------------------------
# what sp to remove
XL_sp <- dat %>% 
  filter(range_size < 2) %>% 
  pull(species) %>% 
  unique()

pdf("./figs/5_col-resids-noXL.pdf")
col_mod_test_abs <- col_dat %>% 
  filter(species %in% XL_sp) %>% 
  select('(Intercept)', bison1, fri,n_grasshop,  spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T)
dev.off()
pdf("./figs/5_cc-resids-noXL.pdf")
cc_mod_test_abs <- cc_dat %>% 
  filter(species %in% XL_sp) %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T)
dev.off()
pdf("./figs/5_ext-resids-noXL.pdf")
ext_mod_test_abs <- ext_dat %>% 
  filter(species %in% XL_sp) %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T)
dev.off()


# which models have range size as a predictor?
res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
res



# plot
# get linear model predictions
pred_col <- bind_rows(col_mod_test_abs$predicted) %>% 
  mutate(mod = "colonization")
pred_cc <- bind_rows(cc_mod_test_abs$predicted) %>% 
  mutate(mod = "midpoint_change")
pred_ext <- bind_rows(ext_mod_test_abs$predicted) %>% 
  mutate(mod = "extinction")
pred_dat <- bind_rows(pred_col, pred_cc, pred_ext)

# make_fig(dat_df = dat, pred_df = pred_dat, sig_coeffs_list = res, height = 9, width = 4, fld = "./figs/5_mod-res-noXLrange.png", filter.sp = T, filter.sp.vec = XL_sp)
# # plot with only drivers of interest
make_fig(dat_df = dat, pred_df = pred_dat,  sig_coeffs_list = res, height = 5.5, width = 4, 
         # filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop", "fri", "n_grasshop_previous", "years_since_last_burn"), 
         filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop"),
         fld = "./figs/5_mod-res-subset-noXL-primary.png", filter.sp = T,  filter.sp.vec = XL_sp)
# # plot with non-drivers of interest
# make_fig(dat_df = dat,pred_df = pred_dat,  sig_coeffs_list = res, height = 6.5, width = 4, 
#          filter.var.vec = c("bison1", "spei", "n_grasshop"), 
#          fld = "./figs/5_mod-res-subset-noXL-secondary.png", filter.sp = T,  filter.sp.vec = XL_sp)


# save
res %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  rename("range.coeff" = 2) %>% 
  mutate(mod = case_when(mod == "cc.bison1.n_grasshop" ~ "cc.bison1:n_grasshop", # change character for interactions
                         mod == "col.bison1.n_grasshop" ~ "col.bison1:n_grasshop",
                         mod == "ext.bison1.n_grasshop" ~ "ext.bison1:n_grasshop",
                         .default = mod)) %>% 
  separate_wider_delim(1, delim = ".", names = c("mod", "driver")) %>% 
  write_csv("./output/5_summary-table-range-coeff-vals-noXLrange.csv")



# control for phylogeny ------------------------------------
library("V.PhyloMaker2")

# species list
# http://lter.konza.ksu.edu/content/pps01-konza-prairie-plant-species-list
sp_info <- read_csv("./data/PPS011.csv") %>%
  rename(epithet = species) %>%
  unite(species, genus, epithet, remove = F, sep = " ") %>%
  mutate(species = str_to_sentence(species))

temp <- sp_info %>% # *** NOTE: using all species found at konza to make our tree. this should give better age estimates
  select(species, genus, family) %>% # reduce
  mutate(old_sp = species) %>% 
  mutate(
    across(c(species, genus, family), str_to_sentence), # capitalize
    family = case_when( # fix old family names
      genus == "Asclepias" ~ "Apocynaceae",
      genus == "Ulmus" ~ "Ulmaceae",
      .default = family
    )
  ) %>%
  distinct() # remove duplicates



## colonization ---------------------------------
# make tree. Species lists are different between the demographic rates. Just going to make each tree separately since the names are different than are in the dat dfs...
tree_out <- temp %>% 
  filter(species %in% col_dat$species,
         # species %in% XL_sp
         ) %>%  # also removing XL range sp
  phylo.maker(., scenarios=c("S1"))
tree <- tree_out$scenario.1
tree <- multi2di(tree) # root tree
# test model
col_mod_test_abs <- col_dat %>% 
  filter(species %in% XL_sp) %>% 
  select('(Intercept)', bison1, fri,n_grasshop,  spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, phylo = tree)

## change in cover ---------------------------------
# make tree. Species lists are different between the demographic rates. Just going to make each tree separately since the names are different than are in the dat dfs...
tree_out <- temp %>% 
  filter(species %in% cc_dat$species,
         # species %in% XL_sp
         ) %>%  # also removing XL range sp
  phylo.maker(., scenarios=c("S1"))
tree <- tree_out$scenario.1
tree <- multi2di(tree) # root tree
# test model
cc_mod_test_abs <- cc_dat %>% 
  filter(
    # species %in% XL_sp
         ) %>% 
  select('(Intercept)', bison1, fri,n_grasshop,  spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, phylo = tree)

## extinction ---------------------------------
# make tree. Species lists are different between the demographic rates. Just going to make each tree separately since the names are different than are in the dat dfs...
tree_out <- temp %>% 
  filter(species %in% ext_dat$species,
         # species %in% XL_sp
         ) %>%  # also removing XL range sp
  phylo.maker(., scenarios=c("S1"))
tree <- tree_out$scenario.1
tree <- multi2di(tree) # root tree
# test model
ext_mod_test_abs <- ext_dat %>% 
  filter(
    # species %in% XL_sp
         ) %>% 
  select('(Intercept)', bison1, fri,n_grasshop,  spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, phylo = tree)



# which models have range size as a predictor?
res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
res

# cant really make a figure since our raw data doesn't match the scale of the model (weights and scaling for phylo)
# so just make a table
var_lookup <- c(
  bison1 = "Bison",
  fri = "Burn interval",
  spei = "Climate",
  n_grasshop = "Grasshopper",
  years_since_last_burn = "Time since fire",
  `bison1.n_grasshop` = "Bison:GH",
  `bison1:n_grasshop` = "Bison:GH",
  bison1_mod = "Bison_mod",
  n_grasshop_mod = "GH_mod"
)

mod_lookup <- c(
  col = "Colonization",
  cc = "Change in cover",
  ext = "Extirpation"
)


res %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  rename("range.coeff" = 2) %>% 
  mutate(mod = case_when(mod == "cc.bison1.n_grasshop" ~ "cc.bison1:n_grasshop", # change character for interactions
                         mod == "col.bison1.n_grasshop" ~ "col.bison1:n_grasshop",
                         mod == "ext.bison1.n_grasshop" ~ "ext.bison1:n_grasshop",
                         .default = mod)) %>% 
  separate_wider_delim(1, delim = ".", names = c("mod", "driver")) %>% 
  filter(driver %in% c('fri', 'spei', 'years_since_last_burn', 'n_grasshop', 'bison1')) %>% # only care about these drivers for manuscript
  mutate(
    mod = recode(mod, !!!mod_lookup), # !!! unpacks the named vector so each named element becomes its own argument for recode()
    driver = recode(driver, !!!var_lookup),
  ) %>% 
  write_csv("./output/5_summary-table-range-coeff-vals-noXLrange-PHYLO.csv")


# Just using grass species ------------------------------------------------
grass_sp <- dat %>% 
  filter(growthform == "g",
         range_size < 1.4 # removing XL range sizes! we have a big outlier. NOTE that this is smaller than the above range size limit
  ) %>% 
  select(species) %>% 
  distinct() %>% 
  pull(species)

col_mod_test_grass <- col_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  filter(species %in% grass_sp) %>%
  mod_test(include_lifespan = F, include_weights = T) # not enough data to use lifespan

cc_mod_test_grass <- cc_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  filter(species %in% grass_sp) %>%
  mod_test(include_lifespan = F, include_weights = T)

ext_mod_test_grass <- ext_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  filter(species %in% grass_sp) %>%
  mod_test(include_lifespan = F, include_weights = T)

# which models have range size as a predictor?
res <- get_sig_coeffs(col_mod_test_grass, cc_mod_test_grass, ext_mod_test_grass)
res

# plot
# get linear model predictions
pred_col <- bind_rows(col_mod_test_grass$predicted) %>% 
  mutate(mod = "colonization")
pred_cc <- bind_rows(cc_mod_test_grass$predicted) %>% 
  mutate(mod = "midpoint_change")
pred_ext <- bind_rows(ext_mod_test_grass$predicted) %>% 
  mutate(mod = "extinction")
pred_dat <- bind_rows(pred_col, pred_cc, pred_ext)


make_fig(dat_df = dat, pred_df = pred_dat, res, height = 5.5, width = 4, fld = "./figs/5_mod-results-just-grass-sp.png", filter.sp = T, filter.sp.vec = grass_sp, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop"))

# save
res %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "mod") %>% 
  rename("range.coeff" = 2) %>% 
  mutate(mod = case_when(mod == "cc.bison1.n_grasshop" ~ "cc.bison1:n_grasshop", # change character for interactions
                         mod == "col.bison1.n_grasshop" ~ "col.bison1:n_grasshop",
                         mod == "ext.bison1.n_grasshop" ~ "ext.bison1:n_grasshop",
                         .default = mod)) %>% 
  separate_wider_delim(1, delim = ".", names = c("mod", "driver")) %>% 
  write_csv("./output/5_summary-table-range-coeff-vals-just-grass-sp.csv")


# GH sp ----------------------------------------------------------
dats <- c("gh-campbell", "gh-mulkern", "gh-welti", "gh")
for (i in 1:length(dats)) {
  if (include_interaction == F) {
    
    
    # colonization models
    col_mods <-  list.files(path = paste0("./output/4_sp-mod-out/", dats[i],"/sp-colonization/."), pattern = ".RDS", full.names = T) %>%
      lapply(., readRDS) %>% 
      setNames( # name list objects by species' name
        str_remove( # remove characters after species name
          str_remove(list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"/sp-colonization/."), pattern = ".RDS", full.names = F), ".*_"), # remove prefix
          "-.*")
      )
    # midpoint change models
    cc_mods <-  list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"/sp-cc-change/."), pattern = ".RDS", full.names = T) %>%
      lapply(., readRDS) %>% 
      setNames( # name list objects by species' name
        str_remove( # remove characters after species name
          str_remove(list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"/sp-cc-change/."), pattern = ".RDS", full.names = F), ".*_"), # remove prefix
          "-.*")
      )
    # extinction models
    ext_mods <-  list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"/sp-extinction/."), pattern = ".RDS", full.names = T) %>%
      lapply(., readRDS) %>% 
      setNames( # name list objects by species' name
        str_remove( # remove characters after species name
          str_remove(list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"/sp-extinction/."), pattern = ".RDS", full.names = F), ".*_"), # remove prefix
          "-.*")
      )
    
    # colonization
    col_dat <- mod_list_to_df(col_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      select('(Intercept)', bison1, fri,n_grasshop,  spei, years_since_last_burn, n, everything()) %>% # organize columns
      mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
    # cc change
    cc_dat <- mod_list_to_df(cc_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, n, everything()) %>% # organize columns
      mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>% # if values are missing, make them zero
      mutate(bison1_mod = bison1+`bison1:n_grasshop`, # modify the bison and GH coef with the interaction value
             n_grasshop_mod = n_grasshop+`bison1:n_grasshop`)
    # extinction
    ext_dat <- mod_list_to_df(ext_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, n, everything()) %>% # organize columns
      mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
    
    
    
    sp <- unique(c(col_dat$species, cc_dat$species, ext_dat$species)) # species included in this dataset
    
    # model testing
    col_mod_test_abs <- col_dat %>% 
      mutate(
        # across(1:7, ~.+1), # no zeros for log
        across(1:7, abs),
        across(1:7, ~.+0.001), # no zeros for log
        across(1:7, log), # response variables are now transformed to log(abs())
      ) %>% 
      mod_test(., plot.title = paste0(dats[i], " log, col"))
    cc_mod_test_abs <- cc_dat %>% 
      mutate(
        # across(1:7, ~.+1), # no zeros for log
        across(1:7, abs),
        across(1:7, ~.+0.001), # no zeros for log
        across(1:7, log), # response variables are now transformed to log(abs())
      ) %>% 
      mod_test(., plot.title = paste0(dats[i], " log cc"))
    ext_mod_test_abs <- ext_dat %>% 
      mutate(
        # across(1:7, ~.+1), # no zeros for log
        across(1:7, abs),
        across(1:7, ~.+0.001), # no zeros for log
        across(1:7, log), # response variables are now transformed to log(abs())
      ) %>% 
      mod_test(.,  plot.title = paste0(dats[i], " log ext"))
    
    
    dat <- col_dat %>% 
      mutate(mod = "colonization") %>% 
      full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
      full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
      left_join(sp_info) %>% # add species info
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, n, everything()) %>%  # organize columns
      pivot_longer(cols = all_of(c('(Intercept)', var_names)), names_to = "var", values_to = "val") %>% 
      mutate(var = as.factor(var),
             val = val,
             # type = case_when(
             #   var %in% c("bison1", "log(n_grasshop)", "log(n_grasshop_previous)") ~ "biotic",
             #   var %in% c("spei", 'fri2', "fri4", 'fri20', "years_since_last_burn", "log(years_since_last_burn + 1)", "fri") ~ "abiotic")
      ) %>% 
      filter(var != "(Intercept)")
  }
  
  if (include_interaction == T) {
    
    # colonization models
    col_mods <-  list.files(path = paste0("./output/4_sp-mod-out/", dats[i],"-with-bison-GH-interaction","/sp-colonization/."), pattern = ".RDS", full.names = T) %>%
      lapply(., readRDS) %>% 
      setNames( # name list objects by species' name
        str_remove( # remove characters after species name
          str_remove(list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"-with-bison-GH-interaction","/sp-colonization/."), pattern = ".RDS", full.names = F), ".*_"), # remove prefix
          "-.*")
      )
    # midpoint change models
    cc_mods <-  list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"-with-bison-GH-interaction","/sp-cc-change/."), pattern = ".RDS", full.names = T) %>%
      lapply(., readRDS) %>% 
      setNames( # name list objects by species' name
        str_remove( # remove characters after species name
          str_remove(list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"-with-bison-GH-interaction","/sp-cc-change/."), pattern = ".RDS", full.names = F), ".*_"), # remove prefix
          "-.*")
      )
    # extinction models
    ext_mods <-  list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"-with-bison-GH-interaction","/sp-extinction/."), pattern = ".RDS", full.names = T) %>%
      lapply(., readRDS) %>% 
      setNames( # name list objects by species' name
        str_remove( # remove characters after species name
          str_remove(list.files(path =  paste0("./output/4_sp-mod-out/", dats[i],"-with-bison-GH-interaction","/sp-extinction/."), pattern = ".RDS", full.names = F), ".*_"), # remove prefix
          "-.*")
      )
    
    # colonization
    col_dat <- mod_list_to_df(col_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
      mutate(
        temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
        bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
        n_grasshop_mod = n_grasshop+temp) %>% 
      select(-temp) %>% 
      {if (filterXL == T) {filter(., species %in% XL_sp)} else .} %>% # remove XL range sizes if needed
      # if zeros are included, we need to make them non-zero before loging
      {if (include_zeros == T) mutate(., 
                                      across(any_of(var_names), abs), 
                                      across(all_of(var_names), ~.+0.001), 
                                      across(all_of(var_names), log)) else mutate(., 
                                                                                  across(all_of(var_names), abs), 
                                                                                  across(all_of(var_names), log))} # no zeros for log. response variables are now transformed to log(abs())
    
    # cc change
    cc_dat <- mod_list_to_df(cc_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
      mutate(
        temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
        bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
        n_grasshop_mod = n_grasshop+temp) %>% 
      select(-temp) %>% 
      {if (filterXL == T) {filter(., species %in% XL_sp)} else .} %>% 
      # if zeros are included, we need to make them non-zero before loging
      {if (include_zeros == T) mutate(., 
                                      across(any_of(var_names), abs), 
                                      across(all_of(var_names), ~.+0.001), 
                                      across(all_of(var_names), log)) else mutate(., 
                                                                                  across(all_of(var_names), abs), 
                                                                                  across(all_of(var_names), log))} # no zeros for log. response variables are now transformed to log(abs())
    
    # extinction
    ext_dat <- mod_list_to_df(ext_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
      mutate(
        temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
        bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
        n_grasshop_mod = n_grasshop+temp) %>% 
      select(-temp) %>% 
      {if (filterXL == T) {filter(., species %in% XL_sp)} else .} %>% # remove XL range sizes if needed
      # if zeros are included, we need to make them non-zero before loging
      {if (include_zeros == T) mutate(., 
                                      across(any_of(var_names), abs), 
                                      across(all_of(var_names), ~.+0.001), 
                                      across(all_of(var_names), log)) else mutate(., 
                                                                                  across(all_of(var_names), abs), 
                                                                                  across(all_of(var_names), log))} # no zeros for log. response variables are now transformed to log(abs())
    
    
    sp <- unique(c(col_dat$species, cc_dat$species, ext_dat$species)) # species included in this dataset
    
    # model testing
    col_mod_test_abs <- col_dat %>% 
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
      mod_test(., plot.title = paste0(dats[i], " log, col"))
    cc_mod_test_abs <- cc_dat %>% 
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
      mod_test(., plot.title = paste0(dats[i], " log cc"))
    ext_mod_test_abs <- ext_dat %>% 
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
      mod_test(.,  plot.title = paste0(dats[i], " log ext"))
    
    ## combine data and format -------------------------------------------------
    dat <- col_dat %>% 
      mutate(mod = "colonization") %>% 
      full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
      full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
      left_join(sp_info) %>% # add species info
      select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%  # organize columns
      pivot_longer(cols = all_of(c('(Intercept)', var_names)), names_to = "var", values_to = "val") %>% 
      mutate(var = as.factor(var),
      ) %>% 
      filter(var != "(Intercept)") 
  }
  
  # which models have range size as a predictor?
  sig_coeffs <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
  # range_size coefficient values for the given driver models
  print(paste0("results for log transformed:", dats[i]))
  print(sig_coeffs)
  
  # save
  sig_coeffs %>% 
    unlist() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "mod") %>% 
    rename("range.coeff" = 2) %>% 
    mutate(mod = case_when(mod == "cc.bison1.n_grasshop" ~ "cc.bison1:n_grasshop", # change character for interactions
                           mod == "col.bison1.n_grasshop" ~ "col.bison1:n_grasshop",
                           mod == "ext.bison1.n_grasshop" ~ "ext.bison1:n_grasshop",
                           .default = mod)) %>% 
    separate_wider_delim(1, delim = ".", names = c("mod", "driver")) %>% 
    write_csv("./output/5_summary-table-range-coeff-vals.csv")
  
  # figure
  # get linear model predictions
  pred_col <- bind_rows(col_mod_test_abs$predicted) %>%
    mutate(mod = "colonization")
  pred_cc <- bind_rows(cc_mod_test_abs$predicted) %>%
    mutate(mod = "midpoint_change")
  pred_ext <- bind_rows(ext_mod_test_abs$predicted) %>%
    mutate(mod = "extinction")
  pred_dat <- bind_rows(pred_col, pred_cc, pred_ext)
  
  
  make_fig(dat_df = dat, pred_df = pred_dat, sig_coeffs, height = 5.5, width = 4, fld = paste0("./figs/5_mod-results-",dats[i], ".png"), filter.sp = T, filter.sp.vec = sp, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop"))
  
  # if(dats[i] == "gh") {
  #   # plot with only drivers of interest
  #   make_fig(sig_coeffs_list = sig_coeffs, height = 4, width = 4, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1:n_grasshop", "bison1.n_grasshop", "fri", "years_since_last_burn"), fld = paste0("./figs/5_mod-results-",dats[i], "-primary.png"), filter.sp = T, filter.sp.vec = sp)
  #   # plot with non-drivers of interest
  #   make_fig(sig_coeffs_list = sig_coeffs, height = 6.5, width = 4, filter.var.vec = c("bison1", "spei","n_grasshop"), fld = paste0("./figs/5_mod-results-",dats[i], "-secondary.png"), filter.sp = T, filter.sp.vec = sp)
  #   
  # }
}




# oddball info for manuscript ---------------------------------------------
if (include_interaction == T) {
  ## interaction effect ------------------------------------------------------
  # colonization models
  col_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # midpoint change models
  cc_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # extinction models
  ext_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
}

if (include_interaction == F) {
  # colonization models
  col_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # midpoint change models
  cc_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
  # extinction models
  ext_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
    lapply(., readRDS) %>% 
    setNames( # name list objects by species' name
      str_remove( # remove characters after species name
        str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
        "-.*")
    )
}

# colonization
col_dat <- mod_list_to_df(col_mods) %>% 
  left_join(ranges) %>%  # add range size for each species
  select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
  {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
  mutate(
    temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
    bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
    n_grasshop_mod = n_grasshop+temp) %>% 
  select(-temp)
# cc change
cc_dat <- mod_list_to_df(cc_mods) %>% 
  left_join(ranges) %>%  # add range size for each species
  select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
  {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
  mutate(
    temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
    bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
    n_grasshop_mod = n_grasshop+temp) %>% 
  select(-temp)
# extinction
ext_dat <- mod_list_to_df(ext_mods) %>% 
  left_join(ranges) %>%  # add range size for each species
  select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
  {if (include_zeros == T) mutate(., across(where(is.numeric), ~replace(., is.na(.), 0))) else .} %>%  # if values are missing, make them zero
  mutate(
    temp = replace_na(`bison1:n_grasshop`, 0), # want to get bison or GH values, even if there is no interaction value
    bison1_mod = bison1+temp, # modify the bison and GH coef with the interaction value
    n_grasshop_mod = n_grasshop+temp) %>% 
  select(-temp) 
# NOTE coeffs are not transformed like they are above

# !@!!! MAKE SURE THIS IS CORRECT WHEN CHANGING THE DATASETS
dat <- col_dat %>% 
  mutate(mod = "colonization") %>% 
  full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
  full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
  left_join(sp_info) %>% # add species info
  select('(Intercept)', bison1, fri,n_grasshop, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%  # organize columns
  pivot_longer(cols = all_of(c('(Intercept)', var_names)), names_to = "var", values_to = "val") %>% 
  mutate(var = as.factor(var),
         # val = val,
         # type = case_when(
         #   var %in% c("bison1", "log(n_grasshop)", "log(n_grasshop_previous)") ~ "biotic",
         #   var %in% c("spei", 'fri2', "fri4", 'fri20', "years_since_last_burn", "log(years_since_last_burn + 1)", "fri") ~ "abiotic")
  ) %>% 
  filter(var != "(Intercept)")

if (filterXL == T) {
  dat <- dat %>% 
    filter(range_size < 2) 
}





# how many coefficients are used for each model
dat %>% 
  group_by(mod, var) %>% 
  summarise(n = sum(!is.na(val))) %>% 
  pivot_wider(names_from = var, values_from = n) %>% 
  write_csv("./output/5_summary-table-n-coeff.csv")

# range sizes
dat %>% 
  summarise(
    med = median(range_size),
    sd = sd(range_size),
    min = min(range_size),
    max = max(range_size)
  )

# effect sizes
dat %>%
  # mutate(val = log(abs(val))) %>% 
  group_by(var, mod) %>%
  summarise(
    med = median(val, na.rm = T),
    mean = mean(val, na.rm = T),
    sd = sd(val, na.rm = T),
    min = min(val, na.rm = T),
    max = max(val, na.rm = T)
  ) %>%
  arrange(desc(mean)) %>% 
  View()
# log transformed
# dat %>%
#   mutate(val = log(abs(val))) %>%
#   group_by(var) %>%
#   summarise(
#     med = median(val, na.rm = T),
#     sd = sd(val, na.rm = T),
#     min = min(val, na.rm = T),
#     max = max(val, na.rm = T)
#   ) %>%
#   arrange(desc(sd))

# number of bison at konza
read_csv("./data/CBH011.csv", show_col_types = FALSE) %>% 
  group_by(RecYear) %>% 
  summarise(n = sum(NumofFemale + NumofMale)) %>%
  # arrange(desc(RecYear))
  filter(RecYear >= 2013 & RecYear <= 2023) %>%
  summarise(mean_bison = mean(n))


#### END
