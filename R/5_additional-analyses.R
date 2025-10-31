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
var_names = c('(Intercept)', "bison1", "fri","n_grasshop", "n_grasshop_previous", "spei", "years_since_last_burn", 'bison1:n_grasshop', "bison1_mod", "n_grasshop_mod")


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
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%  # organize columns
  pivot_longer(cols = 1:10, names_to = "var", values_to = "val") %>% 
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
# primary drivers
dat %>%
  filter(var %in% c("bison1", "fri", "n_grasshop", "n_grasshop_previous", "years_since_last_burn", "spei")) %>%
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



# coefficient values graph ------------------------------------------------

# all
plot_coefficients(dat, filterXL = TRUE, coeff_type = "all")
ggsave("./figs/5_coeff-values-log.png", height = 3, width = 6)

# primary coeffs
plot_coefficients(dat, filterXL = TRUE, coeff_type = "primary")
ggsave("./figs/5_coeff-values-log-primary.png", height = 3, width = 4)

# secondary coeffs
plot_coefficients(dat, filterXL = TRUE, coeff_type = "secondary")
ggsave("./figs/5_coeff-values-log-secondary.png", height = 3, width = 6)




# 
# temp <- dat %>%
#   mutate(var = case_when(var == 'bison1' ~ "Bison",
#                          var == 'fri' ~ "FRI",
#                          # var == 'n_grasshop' ~ "GH[t]",
#                          var == 'n_grasshop' ~ "GH",
#                          var == 'n_grasshop_previous' ~ "GH[t-1]",
#                          var == 'years_since_last_burn' ~ "TSF",
#                          var == 'spei' ~ "SPEI",
#                          var == "bison1:n_grasshop" ~ "Bison:GH",
#                          var == 'bison1_mod' ~ "Bison_mod",
#                          var == 'n_grasshop_mod' ~ "GH_mod"
#   ),
#   ) %>% 
#   {if (filterXL == T) {filter(., range_size < 2)} else .}
# dat_text <- data.frame( # add panel labels
#   label = c("A", "B", "C"),
#   mod   = c("colonization", "midpoint_change", "extinction")
# )
# 
# 
# # histogram to check patterns
# temp %>% 
#   ggplot(aes(x = val)) +
#   geom_histogram() +
#   facet_wrap(~var)
# 
# # all coeffs
# temp %>% 
#   # ggplot( aes(x=var, y=val)) +
#   ggplot( aes(x=var, y=val)) +
#   geom_hline(yintercept = 0, color = "gray70") +
#   geom_jitter(color="gray", size=0.4, alpha=0.9) +
#   geom_violin(alpha = 0.6, fill="gray") +
#   # paletteer::scale_fill_paletteer_d("ggthemes::excel_Median") +
#   labs(x = "Driver", y = "Effect size") +
#   scale_x_discrete("Driver", labels = parse(text = sort(unique(temp$var)))) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, lineheight = 0.5)) +
#   facet_wrap(~factor(mod, levels = c("colonization", "midpoint_change", "extinction"), labels = c("Colonization", "Change in cover", "Extirpation"))) + 
#   geom_text(
#     data    = dat_text,
#     mapping = aes(x = Inf, y = Inf, label = label),
#     hjust   = 2,
#     vjust   = 2,
#     label.size = 0.5,
#     fontface="bold"
#   )
# ggsave("./figs/5_coeff-values-log.png", height = 3, width = 6)
# # primary coeffs
# temp %>% 
#   filter(var %in% c("FRI","GH[t-1]","SPEI","TSF","Bison","GH" )) %>% # remove transformed coeffs
#   # filter(var %in% c("Bison","FRI","GH","GH[t-1]","SPEI","TSF","Bison:GH","Bison_mod","GH_mod" )) %>% 
#   ggplot( aes(x=var, y=val)) +
#   geom_hline(yintercept = 0, color = "gray70") +
#   geom_jitter(color="gray", size=0.4, alpha=0.9) +
#   geom_violin(alpha = 0.6, fill="gray") +
#   # paletteer::scale_fill_paletteer_d("ggthemes::excel_Median") +
#   labs(x = "Driver", y = "Effect size") +
#   scale_x_discrete("Driver", labels = parse(text = sort(unique(filter(temp,var %in% c("FRI","GH[t-1]","SPEI","TSF","Bison","GH" ))$var)))) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, lineheight = 0.5)) +
#   facet_wrap(~factor(mod, levels = c("colonization", "midpoint_change", "extinction"), labels = c("Colonization", "Change in cover", "Extirpation"))) + 
#   geom_text(
#     data    = dat_text,
#     mapping = aes(x = Inf, y = Inf, label = label),
#     hjust   = 2,
#     vjust   = 2,
#     label.size = 0.5,
#     fontface="bold"
#   )
# ggsave("./figs/5_coeff-values-log-primary.png", height = 3, width = 6)
# # secondary coeffs
# temp %>% 
#   filter(var %in% c("Bison_mod","GH_mod", "Bison:GH")) %>% # remove untransformed coeffs
#   # filter(var %in% c("Bison","FRI","GH","GH[t-1]","SPEI","TSF","Bison:GH","Bison_mod","GH_mod" )) %>% 
#   ggplot( aes(x=var, y=val)) +
#   geom_hline(yintercept = 0, color = "gray70") +
#   geom_jitter(color="gray", size=0.4, alpha=0.9) +
#   geom_violin(alpha = 0.6, fill="gray") +
#   # paletteer::scale_fill_paletteer_d("ggthemes::excel_Median") +
#   labs(x = "Driver", y = "Effect size") +
#   scale_x_discrete("Driver", labels = parse(text = sort(unique(filter(temp,var %in% c("Bison_mod","GH_mod", "Bison:GH"))$var)))) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, lineheight = 0.5)) +
#   facet_wrap(~factor(mod, levels = c("colonization", "midpoint_change", "extinction"), labels = c("Colonization", "Change in cover", "Extirpation"))) + 
#   geom_text(
#     data    = dat_text,
#     mapping = aes(x = Inf, y = Inf, label = label),
#     hjust   = 2,
#     vjust   = 2,
#     label.size = 0.5,
#     fontface="bold"
#   )
# ggsave("./figs/5_coeff-values-log-secondary.png", height = 3, width = 6)



# models and model testing using absolute values of coefficients and log transformation ----------------------------------
pdf("./figs/5_col-resids.pdf")
col_mod_test_abs <- col_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, transform_coeffs = transform_coeffs)
dev.off()
pdf("./figs/5_cc-resids.pdf")
cc_mod_test_abs <- cc_dat %>%
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%
  # filter(species != "Gleditsia triacanthos") %>%
  mod_test(include_weights = T, transform_coeffs = transform_coeffs)
dev.off()
pdf("./figs/5_ext-resids.pdf")
ext_mod_test_abs <- ext_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, transform_coeffs = transform_coeffs)
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


make_fig(dat_df = dat, pred_df = pred_dat, res, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1.n_grasshop", "bison1:n_grasshop"), height = 9, width = 4)
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
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, transform_coeffs = transform_coeffs)
dev.off()
pdf("./figs/5_cc-resids-noXL.pdf")
cc_mod_test_abs <- cc_dat %>% 
  filter(species %in% XL_sp) %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, transform_coeffs = transform_coeffs)
dev.off()
pdf("./figs/5_ext-resids-noXL.pdf")
ext_mod_test_abs <- ext_dat %>% 
  filter(species %in% XL_sp) %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  mod_test(include_weights = T, transform_coeffs = transform_coeffs)
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
make_fig(dat_df = dat, pred_df = pred_dat,  sig_coeffs_list = res, height = 4, width = 4, 
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


# Just using grass species ------------------------------------------------
grass_sp <- dat %>% 
  filter(growthform == "g") %>% 
  select(species) %>% 
  distinct() %>% 
  pull(species)

col_mod_test_grass <- col_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  filter(species %in% grass_sp) %>%
  mod_test(include_lifespan = F, include_weights = T, transform_coeffs = transform_coeffs) # not enough data to use lifespan

cc_mod_test_grass <- cc_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  filter(species %in% grass_sp) %>%
  mod_test(include_lifespan = F, include_weights = T, transform_coeffs = transform_coeffs)

ext_mod_test_grass <- ext_dat %>% 
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
  filter(species %in% grass_sp) %>%
  mod_test(include_lifespan = F, include_weights = T, transform_coeffs = transform_coeffs)

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


make_fig(dat_df = dat, pred_df = pred_dat, res, height = 9, width = 4, fld = "./figs/5_mod-results-just-grass-sp.png", filter.sp = T, filter.sp.vec = grass_sp)

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

# # no absolute values or log --------------------------------------------------------
# col_mod_test <- col_dat %>% 
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%  # if values are missing, make them zero
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
#   # filter(
#   #   species != "Gleditsia triacanthos"
#   #   ) %>% 
#   # mutate(across(c(bison1:fri), ~replace(., . > 0, NA))) %>%  # just want negative values. NA values will get removed in the mod_test func
#   mod_test(include_weights = T)
# cc_mod_test <- cc_dat %>% 
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%  # if values are missing, make them zero
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
#   # filter(
#   #   species != "Gleditsia triacanthos"
#   # ) %>% 
#   # mutate(across(c(bison1:fri), ~replace(., . > 0, NA))) %>%  # just want negative values. NA values will get removed in the mod_test func
#   mod_test(include_weights = T)
# ext_mod_test <- ext_dat %>% 
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%  # if values are missing, make them zero
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
#   # filter(
#   #   species != "Gleditsia triacanthos"
#   # ) %>% 
#   # mutate(across(c(bison1:fri), ~replace(., . > 0, NA))) %>%  # just want negative values. NA values will get removed in the mod_test func
#   mod_test(include_weights = T)
# 
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test, cc_mod_test, ext_mod_test)
# res
# # plot
# make_fig(res, height = 7, width = 4, fld = "./figs/5_mod-res-noTransform.png", filter.sp = T, filter.sp.vec = unique(dat$species)[unique(dat$species) != "Gleditsia triacanthos"], val.transform = F)
# 
# # save
# res %>% 
#   unlist() %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var = "mod") %>% 
#   rename("range.coeff" = 2) %>% 
#   separate_wider_delim(1, delim = ".", names = c("mod", "driver")) %>% 
#   write_csv("./output/5_summary-table-range-coeff-vals-noTransform.csv")
# 
# 
# # just negative effects --------------------------------------------------------
# col_mod_test <- col_dat %>% 
#   filter(
#     species != "Gleditsia triacanthos"
#     ) %>%
#   mutate(across(c(bison1:fri), ~replace(., . > 0, NA))) %>%  # just want negative values. NA values will get removed in the mod_test func
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
#   mutate(across(1:7, abs),
#          across(1:7, ~.+1), # no zeros for log
#          across(1:7, log)) %>%
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# cc_mod_test <- cc_dat %>% 
#   filter(
#   species != "Gleditsia triacanthos" # big outlier
#   ) %>%
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
#   mutate(across(c(bison1:fri), ~replace(., . > 0, NA))) %>%  # just want negative values. NA values will get removed in the mod_test func
#   mutate(across(1:7, abs),
#          across(1:7, ~.+1), # no zeros for log
#          across(1:7, log)) %>%
#   filter(range_size < 2) %>%
#   mod_test(include_lifespan = F,include_weights = T)
# ext_mod_test <- ext_dat %>% 
#   filter(
#     species != "Gleditsia triacanthos"
#   ) %>%
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
#   mutate(across(c(bison1:fri), ~replace(., . > 0, NA))) %>%  # just want negative values. NA values will get removed in the mod_test func
#   mutate(across(1:7, abs),
#          across(1:7, ~.+1), # no zeros for log
#          across(1:7, log)) %>%
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# 
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test, cc_mod_test, ext_mod_test)
# res
# 
# 
# # interaction effect ------------------------------------------------------
# # colonization models
# col_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # midpoint change models
# cc_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # extinction models
# ext_mods <-  list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-with-bison-GH-interaction/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )

# # load and format files --
# # colonization
# col_dat <- mod_list_to_df(col_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', n, everything()) %>%  # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# # cc change
# cc_dat <- mod_list_to_df(cc_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', n, everything()) %>%  # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# 
# # extinction
# ext_dat <- mod_list_to_df(ext_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# 
# # models and model testing using absolute values of coefficients and log transformation --
# col_mod_test_abs <- col_dat %>% 
#   mutate(across(1:8, abs),
#          across(1:8, ~.+1), # no zeros for log
#          across(1:8, log), # response variables are now transformed to log(abs())
#          # range_size = log(range_size) 
#   ) %>% 
#   # filter(species != "Gleditsia triacanthos") %>%
#   mod_test(include_weights = T)
# cc_mod_test_abs <- cc_dat %>% 
#   mutate(across(1:8, abs),
#          across(1:8, ~.+1), # no zeros for log
#          across(1:8, log), 
#          # range_size = log(range_size)
#   ) %>% 
#   # filter(species != "Gleditsia triacanthos") %>%
#   mod_test(include_weights = T)
# ext_mod_test_abs <- ext_dat %>% 
#   mutate(across(1:8, abs),
#          across(1:8, ~.+1), # no zeros for log
#          across(1:8, log),
#          # range_size = log(range_size)
#   ) %>% 
#   # filter(species != "Gleditsia triacanthos") %>%
#   mod_test(include_weights = T)
# 
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
# res


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
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
      mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
    # cc change
    cc_dat <- mod_list_to_df(cc_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
      mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>% # if values are missing, make them zero
      mutate(bison1_mod = bison1+`bison1:n_grasshop`, # modify the bison and GH coef with the interaction value
             n_grasshop_mod = n_grasshop+`bison1:n_grasshop`)
    # extinction
    ext_dat <- mod_list_to_df(ext_mods) %>% 
      left_join(ranges) %>%  # add range size for each species
      select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
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
      mod_test(., plot.title = paste0(dats[i], " log, col"), transform_coeffs = transform_coeffs)
    cc_mod_test_abs <- cc_dat %>% 
      mutate(
        # across(1:7, ~.+1), # no zeros for log
        across(1:7, abs),
        across(1:7, ~.+0.001), # no zeros for log
        across(1:7, log), # response variables are now transformed to log(abs())
      ) %>% 
      mod_test(., plot.title = paste0(dats[i], " log cc"), transform_coeffs = transform_coeffs)
    ext_mod_test_abs <- ext_dat %>% 
      mutate(
        # across(1:7, ~.+1), # no zeros for log
        across(1:7, abs),
        across(1:7, ~.+0.001), # no zeros for log
        across(1:7, log), # response variables are now transformed to log(abs())
      ) %>% 
      mod_test(.,  plot.title = paste0(dats[i], " log ext"), transform_coeffs = transform_coeffs)
    
    
    dat <- col_dat %>% 
      mutate(mod = "colonization") %>% 
      full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
      full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
      left_join(sp_info) %>% # add species info
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>%  # organize columns
      pivot_longer(cols = 1:7, names_to = "var", values_to = "val") %>% 
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
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
      mod_test(., plot.title = paste0(dats[i], " log, col"), transform_coeffs = transform_coeffs)
    cc_mod_test_abs <- cc_dat %>% 
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
      mod_test(., plot.title = paste0(dats[i], " log cc"), transform_coeffs = transform_coeffs)
    ext_mod_test_abs <- ext_dat %>% 
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>% 
      mod_test(.,  plot.title = paste0(dats[i], " log ext"), transform_coeffs = transform_coeffs)
    
    # combine data and format -------------------------------------------------
    # !@!!! MAKE SURE THIS IS CORRECT WHEN CHANGING THE DATASETS
    dat <- col_dat %>% 
      mutate(mod = "colonization") %>% 
      full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
      full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
      left_join(sp_info) %>% # add species info
      select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%  # organize columns
      pivot_longer(cols = 1:10, names_to = "var", values_to = "val") %>% 
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
  
  
  make_fig(dat_df = dat, pred_df = pred_dat, sig_coeffs, height = 9, width = 4, fld = paste0("./figs/5_mod-results-",dats[i], ".png"), filter.sp = T, filter.sp.vec = sp)
  
  if(dats[i] == "gh") {
    # plot with only drivers of interest
    make_fig(sig_coeffs_list = sig_coeffs, height = 4, width = 4, filter.var.vec = c("bison1_mod", "n_grasshop_mod", "bison1:n_grasshop", "bison1.n_grasshop", "fri", "n_grasshop_previous", "years_since_last_burn"), fld = paste0("./figs/5_mod-results-",dats[i], "-primary.png"), filter.sp = T, filter.sp.vec = sp)
    # plot with non-drivers of interest
    make_fig(sig_coeffs_list = sig_coeffs, height = 6.5, width = 4, filter.var.vec = c("bison1", "spei","n_grasshop"), fld = paste0("./figs/5_mod-results-",dats[i], "-secondary.png"), filter.sp = T, filter.sp.vec = sp)
    
  }
}

# # bison effects only (no GH coeff) ------------------------------------------------------
# # colonization models
# col_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-bison-coeff-no-GH/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-bison-coeff-no-GH/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # midpoint change models
# cc_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-bison-coeff-no-GH/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-bison-coeff-no-GH/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # extinction models
# ext_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-bison-coeff-no-GH/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-bison-coeff-no-GH/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# 
# # load and format files --
# # colonization
# col_dat <- mod_list_to_df(col_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', bison1, fri, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# 
# # cc change
# cc_dat <- mod_list_to_df(cc_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', bison1, fri, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# # extinction
# ext_dat <- mod_list_to_df(ext_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', bison1, fri, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# 
# # models and model testing using absolute values of coefficients and log transformation --
# col_mod_test_abs <- col_dat %>% 
#   mutate(across(1:5, abs),
#          across(1:5, ~.+1), # no zeros for log
#          across(1:5, log), # response variables are now transformed to log(abs())
#          # range_size = log(range_size) 
#   ) %>% 
#   # filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# cc_mod_test_abs <- cc_dat %>% 
#   mutate(across(1:5, abs),
#          across(1:5, ~.+1), # no zeros for log
#          across(1:5, log), 
#          # range_size = log(range_size)
#   ) %>% 
#   # filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# ext_mod_test_abs <- ext_dat %>% 
#   mutate(across(1:5, abs),
#          across(1:5, ~.+1), # no zeros for log
#          across(1:5, log),
#          # range_size = log(range_size)
#   ) %>% 
#   # filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
# res
# # plot
# make_fig(res, height = 7, width = 4, 
#          fld = "./figs/5_mod-res-Bison-noGH.png", 
#          filter.sp = T, filter.sp.vec = XL_sp
# )
# 
# # GH effects only (no bison grazing) ------------------------------------------------------
# # colonization models
# col_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-no-bison/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-no-bison/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # midpoint change models
# cc_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-no-bison/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-no-bison/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # extinction models
# ext_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-no-bison/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-no-bison/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# 
# # format files --
# # colonization
# col_dat <- mod_list_to_df(col_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', fri, n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# # cc change
# cc_dat <- mod_list_to_df(cc_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# # extinction
# ext_dat <- mod_list_to_df(ext_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# 
# # combine data and format --
# dat <- col_dat %>% 
#   mutate(mod = "colonization") %>% 
#   full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
#   full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
#   left_join(sp_info) %>% # add species info
#   pivot_longer(cols = 1:6, names_to = "var", values_to = "val") %>% 
#   mutate(var = as.factor(var),
#          val = val,
#          type = case_when(
#            var %in% c("bison1", "log(n_grasshop)", "log(n_grasshop_previous)") ~ "biotic",
#            var %in% c("spei", 'fri2', "fri4", 'fri20', "years_since_last_burn", "log(years_since_last_burn + 1)", "fri") ~ "abiotic"
#          )) %>% 
#   filter(var != "(Intercept)")
# 
# # models and model testing using absolute values of coefficients and log transformation --
# col_mod_test_abs <- col_dat %>% 
#   mutate(across(1:6, abs),
#          across(1:6, ~.+1), # no zeros for log
#          across(1:6, log), # response variables are now transformed to log(abs())
#          # range_size = log(range_size) 
#   ) %>% 
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# cc_mod_test_abs <- cc_dat %>% 
#   mutate(across(1:6, abs),
#          across(1:6, ~.+1), # no zeros for log
#          across(1:6, log), 
#          # range_size = log(range_size)
#   ) %>% 
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# ext_mod_test_abs <- ext_dat %>% 
#   mutate(across(1:6, abs),
#          across(1:6, ~.+1), # no zeros for log
#          across(1:6, log),
#          # range_size = log(range_size)
#   ) %>% 
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
# res
# # plot
# make_fig(res, height = 7, width = 4, 
#          fld = "./figs/5_mod-res-GH-noBisonPresent.png", 
#          filter.sp = T, filter.sp.vec = XL_sp
#          )
# 
# 
# # GH effects only (yes bison grazing) ------------------------------------------------------
# # colonization models
# col_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-yes-bison/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-yes-bison/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # midpoint change models
# cc_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-yes-bison/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-yes-bison/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # extinction models
# ext_mods <-  list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-yes-bison/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all-just-GH-coeff-yes-bison/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# 
# # format files --
# # colonization
# col_dat <- mod_list_to_df(col_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# # cc change
# cc_dat <- mod_list_to_df(cc_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# # extinction
# ext_dat <- mod_list_to_df(ext_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   select('(Intercept)', fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) %>% # organize columns
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) # if values are missing, make them zero
# 
# 
# # combine data and format --
# dat <- col_dat %>% 
#   mutate(mod = "colonization") %>% 
#   full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
#   full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
#   left_join(sp_info) %>% # add species info
#   pivot_longer(cols = 1:6, names_to = "var", values_to = "val") %>% 
#   mutate(var = as.factor(var),
#          val = val,
#          type = case_when(
#            var %in% c("bison1", "log(n_grasshop)", "log(n_grasshop_previous)") ~ "biotic",
#            var %in% c("spei", 'fri2', "fri4", 'fri20', "years_since_last_burn", "log(years_since_last_burn + 1)", "fri") ~ "abiotic"
#          )) %>% 
#   filter(var != "(Intercept)")
# 
# # models and model testing using absolute values of coefficients and log transformation --
# col_mod_test_abs <- col_dat %>% 
#   mutate(across(1:6, abs),
#          across(1:6, ~.+1), # no zeros for log
#          across(1:6, log), # response variables are now transformed to log(abs())
#          # range_size = log(range_size) 
#   ) %>% 
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# cc_mod_test_abs <- cc_dat %>% 
#   mutate(across(1:6, abs),
#          across(1:6, ~.+1), # no zeros for log
#          across(1:6, log), 
#          # range_size = log(range_size)
#   ) %>% 
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# ext_mod_test_abs <- ext_dat %>% 
#   mutate(across(1:6, abs),
#          across(1:6, ~.+1), # no zeros for log
#          across(1:6, log),
#          # range_size = log(range_size)
#   ) %>% 
#   filter(range_size < 2) %>%
#   mod_test(include_weights = T)
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
# res
# # plot
# make_fig(res, height = 7, width = 4, 
#          fld = "./figs/5_mod-res-GH-yesBisonPresent.png",
#          filter.sp = T, filter.sp.vec = XL_sp
#          )
# 
# 
# # instead of range size, use taxon age ------------------------------------
# library("V.PhyloMaker2")
# library("ggtree")
# 
# # reload files from start
# # plant info
# # http://lter.konza.ksu.edu/content/pps01-konza-prairie-plant-species-list
# sp_info <- read_csv("./data/PPS011.csv") %>% 
#   rename(epithet = species) %>% 
#   unite(species, genus, epithet, remove = F, sep = " ") %>% 
#   mutate(species = str_to_sentence(species))
# # colonization models
# col_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-colonization/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-colonization/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # midpoint change models
# cc_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-cc-change/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-cc-change/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# # extinction models
# ext_mods <-  list.files(path = "./output/4_sp-mod-out/all/sp-extinction/.", pattern = ".RDS", full.names = T) %>%
#   lapply(., readRDS) %>% 
#   setNames( # name list objects by species' name
#     str_remove( # remove characters after species name
#       str_remove(list.files(path = "./output/4_sp-mod-out/all/sp-extinction/.", pattern = ".RDS", full.names = F), ".*_"), # remove prefix
#       "-.*")
#   )
# 
# 
# # load and format files --
# # colonization
# col_dat <- mod_list_to_df(col_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%  # if values are missing, make them zero
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) # organize columns
# 
# # cc change
# cc_dat <- mod_list_to_df(cc_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%  # if values are missing, make them zero
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) # organize columns
# # extinction
# ext_dat <- mod_list_to_df(ext_mods) %>% 
#   left_join(ranges) %>%  # add range size for each species
#   select(!starts_with("watershed")) %>%  # remove watershed coeficients if needed
#   mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%  # if values are missing, make them zero
#   select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, n, everything()) # organize columns
# # format species names
# temp <- sp_info %>% # *** NOTE: using all species found at konza to make our tree. this should give better age estimates
#   select(species, genus, family) %>% # reduce
#   mutate(
#     across(everything(), str_to_sentence), # capitalize
#     family = case_when( # fix old family names
#       genus == "Asclepias" ~ "Apocynaceae",
#       genus == "Ulmus" ~ "Ulmaceae",
#       .default = family
#     )
#   ) %>%  
#   distinct() # remove duplicates
# 
# # make tree
# result <- phylo.maker(temp, scenarios=c("S1"))
# 
# tree.summary <- result$species.list # info about how species were matched to the tree
# tree <- result$scenario.1
# 
# ### plot the phylogenies with node ages displayed.
# ggtree(tree) +
#   theme_tree2() + # add scale
#   geom_tiplab(align=TRUE, size = 1) # tip labels
# 
# # get genus age. Note that these values are time since root (so smaller values are older lineages)
# # using our tree *which is limited by how many species we have
# temp <- temp %>% 
#   mutate(species = str_replace(species, " ", "_")) # reformat sp to match tree
# age <-build.nodes.2(tree, temp)[,c("genus", "bn.bl")] # build nodes .1 and .2 give the same numbers
# age[age == ""] <- NA # blanks need to be NA
# age <- drop_na(age)
# head(arrange(age, bn.bl))
# 
# # information from the mega-tree that our smaller tree is based on
# age.bigtree <- nodes.info.1.TPL %>%  # mega tree data summary
#   filter(genus %in% unique(temp$genus)) %>% 
#   select(genus, bn.bl)
# head(arrange(age.bigtree, bn.bl))
# 
# # compare ages between the methods
# age %>% 
#   rename(local_tree_age = bn.bl) %>% 
#   full_join(age.bigtree) %>% 
#   rename(big_tree_age = bn.bl) %>% 
#   ggplot(aes(x = local_tree_age, y = big_tree_age)) +
#   geom_abline(intercept = 0, slope = 1, color = "lightblue") +
#   geom_point() 
# 
# # # show lengths
# # ggtree(tree) +
# #   theme_tree2() + # add scale
# #   geom_tiplab(align=TRUE, linesize=.5) + # tip labels
# #   geom_text2(aes(label = round(branch.length, 2)), hjust = 1, vjust = 1.2, size = 2) +
# #   ggplot2::xlim(0, 175)
# 
# 
# # models and model testing using absolute values of coefficients and log transformation 
# col_mod_test_abs <- col_dat %>% 
#   mutate(genus = word(species, 1)) %>%  # get genus name to match genus age with
#   left_join(age.bigtree) %>% 
#   select(!c(range_size, genus)) %>% 
#   rename(range_size = bn.bl) %>%  # !!! NOTE: I am just renaming the genus age as range_size because I am lazy and dont want to fix my function
#   mutate(across(1:7, abs),
#          across(1:7, ~.+1), # no zeros for log
#          across(1:7, log), # response variables are now transformed to log(abs())
#          # range_size = log(range_size) 
#   ) %>% 
#   mod_test(include_weights = T)
# cc_mod_test_abs <- cc_dat %>% 
#   mutate(genus = word(species, 1)) %>%  # get genus name to match genus age with
#   left_join(age.bigtree) %>% 
#   select(!c(range_size, genus)) %>% 
#   rename(range_size = bn.bl) %>%  # !!! NOTE: I am just renaming the genus age as range_size because I am lazy and dont want to fix my function
#   mutate(across(1:7, abs),
#          across(1:7, ~.+1), # no zeros for log
#          across(1:7, log), 
#          # range_size = log(range_size)
#   ) %>% 
#   # filter(species != "Gleditsia triacanthos") %>% 
#   mod_test(include_weights = T)
# ext_mod_test_abs <- ext_dat %>%   
#   mutate(genus = word(species, 1)) %>%  # get genus name to match genus age with
#   left_join(age.bigtree) %>% 
#   select(!c(range_size, genus)) %>% 
#   rename(range_size = bn.bl) %>%  # !!! NOTE: I am just renaming the genus age as range_size because I am lazy and dont want to fix my function
#   mutate(across(1:7, abs),
#          across(1:7, ~.+1), # no zeros for log
#          across(1:7, log),
#          # range_size = log(range_size)
#   ) %>% 
#   mod_test(include_weights = T)
# 
# # which models have range size as a predictor?
# res <- get_sig_coeffs(col_mod_test_abs, cc_mod_test_abs, ext_mod_test_abs)
# res
# 
# # plot
# # combine data and format 
# dat <- col_dat %>% 
#   mutate(mod = "colonization") %>% 
#   full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
#   full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
#   left_join(sp_info) %>% # add species info
#   pivot_longer(cols = 1:7, names_to = "var", values_to = "val") %>% 
#   mutate(var = as.factor(var),
#          val = val,
#          type = case_when(
#            var %in% c("bison1", "log(n_grasshop)", "log(n_grasshop_previous)") ~ "biotic",
#            var %in% c("spei", 'fri2', "fri4", 'fri20', "years_since_last_burn", "log(years_since_last_burn + 1)", "fri") ~ "abiotic"
#          )) %>% 
#   filter(var != "(Intercept)") %>% 
#   mutate(genus = str_to_sentence(genus)) %>% 
#   left_join(age.bigtree) %>% 
#   select(!range_size) %>% # NOTE !! I AM STILL USING THE TITLE RANGE SIZE DESPITE IT NOW BEING GENUS AGE
#   rename(range_size = bn.bl)
# 
# make_fig(res, height = 7, width = 4, fld = "./figs/5_mod-res-genus-age-NOTRANGESIZE.png", filter.sp = T, filter.sp.vec = XL_sp)
# 
# # range size versus genus age
# ranges %>% 
#   mutate(genus = word(species, 1)) %>% 
#   left_join(age.bigtree) %>% 
#   ggplot(aes(y = range_size, x = bn.bl)) +
#   geom_point() +
#   labs( y = "Range size", x = "Genus age")

# oddball info for manuscript ---------------------------------------------
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
  select('(Intercept)', bison1, fri,n_grasshop, n_grasshop_previous, spei, years_since_last_burn, 'bison1:n_grasshop', bison1_mod, n_grasshop_mod, n, everything()) %>%  # organize columns
  pivot_longer(cols = 1:10, names_to = "var", values_to = "val") %>% 
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


# OLD ---------------------------------------------------------------------
# mod_test_abs <- function(dat, plot.title = NULL){ # helper function
#   out1 <- list() # list to fill with full dredge output
#   out2 <- list() # list to fill with just coefficients from dredge output
#   params <- (ncol(dat)-2) # number of parameters
#   message(paste("detected ", params, " parameters"))
#   for (col in 2:params) { # iterate over columns (drivers)
#     temp <- dat %>% 
#       select(all_of(col), range_size) %>% 
#       drop_na()
#     message(paste0("\n",colnames(temp)[1], " has n = ", nrow(temp) ," observations for the model")) # helpful info
#     colnames(temp)[1] <-  sub("log\\(", "", colnames(temp)[1]) # rename columns to remove any transformation notation (e.g. log())
#     colnames(temp)[1] <-  sub(" \\+ 1\\)", "", colnames(temp)[1]) 
#     colnames(temp)[1] <-  sub("\\)", "", colnames(temp)[1]) 
#     
#     mod <- lm(as.formula(paste0("abs(", names(temp)[1], ")", " ~ log(range_size)")), data = temp, na.action = "na.fail") # make model
#     # residuals
#     simulationOutput <- DHARMa::simulateResiduals(fittedModel = mod, plot = F)
#     plot(simulationOutput, title = paste0(plot.title, " ", colnames(temp)[1]))
#     
#     out1[[col]] <- MuMIn::dredge(mod) # test with AICc
#     names(out1)[col] <- names(temp)[1] # name object in list as driver
#     
#     out2[[col]] <- MuMIn::get.models(out1[[col]], subset = 1)[[1]]$coefficients # get best model
#     names(out2)[col] <- names(temp)[1] # name object in list as driver
#     
#     #print(paste0("done with parameter ", col))
#   }
#   out <- list(out1, out2)
#   names(out) <- c("full_dredge_out", "coef_dredge")
#   return(out)
# }
## # ## Correlation -------------------------------------------------------------
# cor_test <- function(dat){ # helper function
#   params <- (ncol(dat)-2) # number of parameters
#   message(paste("detected ", params, " parameters"))
#   out <- data.frame(ci_low = numeric(),
#                     ci_high = numeric(),
#                     R = numeric(),
#                     t = numeric(),
#                     driver = character(),
#                     sig = character(),
#                     stringsAsFactors=FALSE) # list to fill with output
#   for (col in 2:params) { # iterate over columns (drivers)
#     temp <- dat %>%
#       select(all_of(col), range_size) %>%
#       drop_na()
#     colnames(temp)[1] <-  sub("log\\(", "", colnames(temp)[1]) # rename columns to remove any transformation notation (e.g. log())
#     colnames(temp)[1] <-  sub(" \\+ 1\\)", "", colnames(temp)[1])
#     colnames(temp)[1] <-  sub("\\)", "", colnames(temp)[1])
#     
#     
#     mod <- cor.test(temp[,1], log(temp[,2])) # run correlation test
#     
#     out[col-1,1:2] <- mod$conf.int[1:2] # confidence interval
#     out[col-1,3] <- mod$estimate # correlation estimate
#     out[col-1,4] <- mod$statistic # t-stat
#     out[col-1,5] <- names(temp)[1] # driver
#     out[col-1,6] <- ifelse(mod$conf.int[1] > 0 & mod$conf.int[2] > 0 | mod$conf.int[1] < 0 & mod$conf.int[2] < 0, "*", NA)
#   }
#   return(out)
# }
# # 
# # cor_test(col_dat)
# # cor_test(cc_dat)
# # cor_test(ext_dat)
# # 
## # Density structured models -----------------------------------------------
# 
# 
## ## Bison -------------------------------------------------------------------
# density_out <- read_csv("./output/4_density-model-output-bison.csv") %>% 
#   left_join(ranges) %>% 
#   mutate(effect = abs(effect)) %>% 
#   filter(effect != Inf) # we have a few species that have denominators of zero (dont occure in ungrazed), leading to Inf values
# 
# density_out %>% 
#   ggplot(aes(x = log(range_size), y = effect)) +
#   geom_hline(yintercept = 0, color = "gray80") +
#   geom_point(color = "darkgray") +
#   # geom_smooth(method = "lm", se = F) +
#   # scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
#   labs(x = "Range size (log)", y = "Effect size (log)")
# ggsave("./figs/5_mod-results-SUP-density-mods-bison.png", height = 2, width = 2)
# 
# mod <- lm(effect ~ log(range_size), data = density_out, na.action = "na.fail") # make model
# dredge_out <- MuMIn::dredge(mod) # test with AICc
# MuMIn::get.models(dredge_out, subset = 1)[[1]]$coefficients # get best model
# 
# #
## ### multiply coeff values by sd ---------------------------------------------
# # dvr_sd <- read_csv("./output/4_driver-sd.csv") %>% 
# #   mutate(driver_sd = case_when(driver == "n_grasshop" ~ log(driver_sd), # format to match models
# #                                driver == "n_grasshop_previous" ~ log(driver_sd),
# #                                driver == "years_since_last_burn" ~ log(driver_sd+1),
# #                                .default = driver_sd),
# #          driver = case_when(driver == "n_grasshop" ~ 'log(n_grasshop)',
# #                             driver == 'n_grasshop_previous' ~ 'log(n_grasshop_previous)',
# #                             driver == 'years_since_last_burn' ~ 'log(years_since_last_burn + 1)',
# #                             .default = driver
# #                             ))
# # 
# # 
# # col_mod_test_abs <- col_dat %>% 
# #   mutate(`log(n_grasshop)` = `log(n_grasshop)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop)',2]),
# #          `log(n_grasshop_previous)` = `log(n_grasshop_previous)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop_previous)',2]),
# #          `log(years_since_last_burn + 1)` = `log(years_since_last_burn + 1)`*pull(dvr_sd[dvr_sd[,1] == 'log(years_since_last_burn + 1)',2]),
# #          spei = spei*pull(dvr_sd[dvr_sd[,1] == 'spei',2])) %>% 
# #   mod_test_abs()
# # cc_mod_test_abs <- cc_dat %>% 
# #   mutate(`log(n_grasshop)` = `log(n_grasshop)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop)',2]),
# #          `log(n_grasshop_previous)` = `log(n_grasshop_previous)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop_previous)',2]),
# #          `log(years_since_last_burn + 1)` = `log(years_since_last_burn + 1)`*pull(dvr_sd[dvr_sd[,1] == 'log(years_since_last_burn + 1)',2]),
# #          spei = spei*pull(dvr_sd[dvr_sd[,1] == 'spei',2])) %>% 
# #   mod_test_abs()
# # ext_mod_test_abs <- ext_dat %>% 
# #   mutate(`log(n_grasshop)` = `log(n_grasshop)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop)',2]),
# #          `log(n_grasshop_previous)` = `log(n_grasshop_previous)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop_previous)',2]),
# #          `log(years_since_last_burn + 1)` = `log(years_since_last_burn + 1)`*pull(dvr_sd[dvr_sd[,1] == 'log(years_since_last_burn + 1)',2]),
# #          spei = spei*pull(dvr_sd[dvr_sd[,1] == 'spei',2])) %>% 
# #   mod_test_abs()
# # # which models have range size as a predictor?
# # # colonization
# # temp <- unlist(col_mod_test_abs['coef_dredge'])
# # names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# # temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# # names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# # col_range_coeff_abs <- temp
# # # cc change
# # temp <- unlist(cc_mod_test_abs['coef_dredge'])
# # names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# # temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# # names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# # cc_range_coeff_abs <- temp
# # # extinction
# # temp <- unlist(ext_mod_test_abs['coef_dredge'])
# # names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# # temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# # names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# # ext_range_coeff_abs <- temp
# # # range_size coefficient values for the given driver models
# # col_range_coeff_abs
# # cc_range_coeff_abs
# # ext_range_coeff_abs
# # 
# # # New facet label names for supp variable
# # names(col_range_coeff_abs) <- case_when(names(col_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
# #                                         names(col_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
# #                                         names(col_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
# #                                         .default = names(col_range_coeff_abs))
# # names(cc_range_coeff_abs) <- case_when(names(cc_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
# #                                        names(cc_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
# #                                        names(cc_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
# #                                        .default = names(cc_range_coeff_abs))
# # names(ext_range_coeff_abs) <- case_when(names(ext_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
# #                                         names(ext_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
# #                                         names(ext_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
# #                                         .default = names(ext_range_coeff_abs))
# # f_labels <- bind_rows(mutate(as_tibble(col_range_coeff_abs, rownames = "var"), # add range coefficients
# #                              mod = "colonization"),
# #                       mutate(as_tibble(cc_range_coeff_abs, rownames = "var"),
# #                              mod = "midpoint_change"),
# #                       mutate(as_tibble(ext_range_coeff_abs, rownames = "var"),
# #                              mod = "extinction")) %>%   
# #   rename(range_coeff = value) %>% 
# #   mutate(
# #     range_coeff = round(range_coeff, 2),
# #     label = paste0(range_coeff) #"Range coef.:", 
# #   ) %>% 
# #   filter(var %in% c("log(n_grasshop)", "log(n_grasshop_previous)", "log(years_since_last_burn + 1)", "spei")) # remove unmodified drivers
# # 
# # dat %>% 
# #   left_join(rename(dvr_sd, var = driver)) %>% 
# #   filter(var %in% dvr_sd$driver) %>% 
# #   mutate(sig = case_when(mod == 'colonization' & var %in% names(col_range_coeff_abs) ~ T, # add column to indicate if the driver is significant
# #                          mod == 'midpoint_change' & var %in% names(cc_range_coeff_abs) ~ T,
# #                          mod == 'extinction' & var %in% names(ext_range_coeff_abs) ~ T),
# #          modified_val = val*driver_sd
# #   ) %>% 
# # # plot
# #   ggplot(aes(x = log(range_size), y = modified_val)) +
# #   geom_hline(yintercept = 0, color = "gray80") +
# #   geom_point(color = "darkgray") +
# #   geom_smooth(method = "lm", aes(color = sig), se = F) +
# #   scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
# #   facet_grid(var~factor(mod, levels = c("colonization", "midpoint_change", "extinction"), labels = c("Colonization", "Change in cover", "Extinction")), # for whatever reason when I try to order 'mod' this is the only setup where i can get them in order AND with the pretty names
# #              scales = "free_y",
# #              labeller = labeller(mod = mod.labs, var = var.labs)
# #   ) +
# #   geom_text(
# #     data    = f_labels,
# #     mapping = aes(x = -Inf, y = Inf, label = label),
# #     hjust   = -0.25,
# #     vjust   = 2
# #   ) +
# #   labs(x = "Range size (log)", y = "Effect size")
# # ggsave("./figs/5_mod-results-SUP-sd.png", height = 8, width = 6)
# 
# 
# 
## ## GH -------------------------------------------------------------------
# density_out <- read_csv("./output/4_density-model-output-gh.csv") %>% 
#   left_join(ranges) %>% 
#   mutate(effect = abs(effect)) %>% 
#   filter(effect != Inf) # we have a few species that have denominators of zero (dont occure in ungrazed), leading to Inf values
# 
# density_out %>% 
#   ggplot(aes(x = log(range_size), y = effect)) +
#   geom_hline(yintercept = 0, color = "gray80") +
#   geom_point(color = "darkgray") +
#   # geom_smooth(method = "lm", se = F) +
#   # scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
#   labs(x = "Range size (log)", y = "Effect size (log)")
# ggsave("./figs/5_mod-results-SUP-density-mods-gh.png", height = 2, width = 2)
# 
# mod <- lm(effect ~ log(range_size), data = density_out, na.action = "na.fail") # make model
# dredge_out <- MuMIn::dredge(mod) # test with AICc
# MuMIn::get.models(dredge_out, subset = 1)[[1]]$coefficients # get best model
# 
## # multiply coeff values by sd AND loging---------------------------------------------
# dvr_sd <- read_csv("./output/4_driver-sd.csv") %>% 
#   mutate(driver_sd = case_when(driver == "n_grasshop" ~ log(driver_sd), # format to match models
#                                driver == "n_grasshop_previous" ~ log(driver_sd),
#                                driver == "years_since_last_burn" ~ log(driver_sd+1),
#                                .default = driver_sd),
#          driver = case_when(driver == "n_grasshop" ~ 'log(n_grasshop)',
#                             driver == 'n_grasshop_previous' ~ 'log(n_grasshop_previous)',
#                             driver == 'years_since_last_burn' ~ 'log(years_since_last_burn + 1)',
#                             .default = driver
#          ))
# 
# 
# col_mod_test_abs <- col_dat %>% 
#   mutate(`log(n_grasshop)` = `log(n_grasshop)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop)',2]),
#          `log(n_grasshop_previous)` = `log(n_grasshop_previous)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop_previous)',2]),
#          `log(years_since_last_burn + 1)` = `log(years_since_last_burn + 1)`*pull(dvr_sd[dvr_sd[,1] == 'log(years_since_last_burn + 1)',2]),
#          spei = spei*pull(dvr_sd[dvr_sd[,1] == 'spei',2])) %>% 
#   mutate(across(1:7, abs),
#          across(1:7, log)) %>% 
#   mod_test()
# cc_mod_test_abs <- cc_dat %>% 
#   mutate(`log(n_grasshop)` = `log(n_grasshop)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop)',2]),
#          `log(n_grasshop_previous)` = `log(n_grasshop_previous)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop_previous)',2]),
#          `log(years_since_last_burn + 1)` = `log(years_since_last_burn + 1)`*pull(dvr_sd[dvr_sd[,1] == 'log(years_since_last_burn + 1)',2]),
#          spei = spei*pull(dvr_sd[dvr_sd[,1] == 'spei',2])) %>% 
#   mutate(across(1:7, abs),
#          across(1:7, log)) %>% 
#   mod_test()
# ext_mod_test_abs <- ext_dat %>% 
#   mutate(`log(n_grasshop)` = `log(n_grasshop)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop)',2]),
#          `log(n_grasshop_previous)` = `log(n_grasshop_previous)`*pull(dvr_sd[dvr_sd[,1] == 'log(n_grasshop_previous)',2]),
#          `log(years_since_last_burn + 1)` = `log(years_since_last_burn + 1)`*pull(dvr_sd[dvr_sd[,1] == 'log(years_since_last_burn + 1)',2]),
#          spei = spei*pull(dvr_sd[dvr_sd[,1] == 'spei',2])) %>% 
#   mutate(across(1:7, abs),
#          across(1:7, log)) %>% 
#   mod_test()
# # which models have range size as a predictor?
# # colonization
# temp <- unlist(col_mod_test_abs['coef_dredge'])
# names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# col_range_coeff_abs <- temp
# # cc change
# temp <- unlist(cc_mod_test_abs['coef_dredge'])
# names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# cc_range_coeff_abs <- temp
# # extinction
# temp <- unlist(ext_mod_test_abs['coef_dredge'])
# names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# ext_range_coeff_abs <- temp
# # range_size coefficient values for the given driver models
# col_range_coeff_abs
# cc_range_coeff_abs
# ext_range_coeff_abs
# 
# # New facet label names for supp variable
# names(col_range_coeff_abs) <- case_when(names(col_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
#                                         names(col_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
#                                         names(col_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
#                                         .default = names(col_range_coeff_abs))
# names(cc_range_coeff_abs) <- case_when(names(cc_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
#                                        names(cc_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
#                                        names(cc_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
#                                        .default = names(cc_range_coeff_abs))
# names(ext_range_coeff_abs) <- case_when(names(ext_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
#                                         names(ext_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
#                                         names(ext_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
#                                         .default = names(ext_range_coeff_abs))
# f_labels <- bind_rows(mutate(as_tibble(col_range_coeff_abs, rownames = "var"), # add range coefficients
#                              mod = "colonization"),
#                       mutate(as_tibble(cc_range_coeff_abs, rownames = "var"),
#                              mod = "midpoint_change"),
#                       mutate(as_tibble(ext_range_coeff_abs, rownames = "var"),
#                              mod = "extinction")) %>%   
#   rename(range_coeff = value) %>% 
#   mutate(
#     range_coeff = round(range_coeff, 2),
#     label = paste0(range_coeff) #"Range coef.:", 
#   ) %>% 
#   filter(var %in% c("log(n_grasshop)", "log(n_grasshop_previous)", "log(years_since_last_burn + 1)", "spei")) # remove unmodified drivers
# 
# dat %>% 
#   left_join(rename(dvr_sd, var = driver)) %>% 
#   filter(var %in% dvr_sd$driver) %>% 
#   mutate(sig = case_when(mod == 'colonization' & var %in% names(col_range_coeff_abs) ~ T, # add column to indicate if the driver is significant
#                          mod == 'midpoint_change' & var %in% names(cc_range_coeff_abs) ~ T,
#                          mod == 'extinction' & var %in% names(ext_range_coeff_abs) ~ T),
#          modified_val = val*driver_sd
#   ) %>% 
#   # plot
#   ggplot(aes(x = log(range_size), y = log(modified_val))) +
#   geom_hline(yintercept = 0, color = "gray80") +
#   geom_point(color = "darkgray") +
#   geom_smooth(method = "lm", aes(color = sig), se = F) +
#   scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
#   facet_grid(var~factor(mod, levels = c("colonization", "midpoint_change", "extinction"), labels = c("Colonization", "Change in cover", "Extinction")), # for whatever reason when I try to order 'mod' this is the only setup where i can get them in order AND with the pretty names
#              scales = "free_y",
#              labeller = labeller(mod = mod.labs, var = var.labs)
#   ) +
#   geom_text(
#     data    = f_labels,
#     mapping = aes(x = -Inf, y = Inf, label = label),
#     hjust   = -0.25,
#     vjust   = 2
#   ) +
#   labs(x = "Range size (log)", y = "Effect size (log)")
# ggsave("./figs/5_mod-results-SUP-sd-log.png", height = 8, width = 6)
# 
## # # correlation between drivers ---------------------------------------------
# # panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
# #   usr <- par("usr")
# #   on.exit(par(usr))
# #   par(usr = c(0, 1, 0, 1))
# #   Cor <- abs(cor(x, y)) # Remove abs function if desired
# #   txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
# #   if(missing(cex.cor)) {
# #     cex.cor <- 0.4 / strwidth(txt)
# #   }
# #   text(0.5, 0.5, txt,
# #        cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
# # }
# # 
# # 
# # temp <- col_dat %>% 
# #   mutate(mod = "colonization") %>% 
# #   full_join(cc_dat %>% mutate(mod = "midpoint_change"))%>% 
# #   full_join(ext_dat %>% mutate(mod = "extinction")) %>% 
# #   select(mod, 2:7) 
# # # col
# # temp %>% 
# #   filter(mod == 'colonization') %>% 
# #   select(!mod) %>% 
# #   mutate_all(abs) %>% # make all values positive
# #   pairs(.,                     # Data frame of variables
# #       labels = colnames(.),  # Variable names
# #       upper.panel = panel.cor,    # Correlation panel
# #       main = "Colonization"
# #       # lower.panel = panel.smooth) # Smoothed regression lines                 
# #   )
# # # cc
# # temp %>% 
# #   filter(mod == 'midpoint_change') %>% 
# #   select(!mod) %>% 
# #   mutate_all(abs) %>% # make all values positive
# #   pairs(.,                     # Data frame of variables
# #         labels = colnames(.),  # Variable names
# #         upper.panel = panel.cor,    # Correlation panel
# #         main = "Change in cover"
# #         # lower.panel = panel.smooth) # Smoothed regression lines     
# #   )
# # # ext
# # temp %>% 
# #   filter(mod == 'extinction') %>% 
# #   select(!mod) %>% 
# #   mutate_all(abs) %>% # make all values positive
# #   pairs(.,                     # Data frame of variables
# #         labels = colnames(.),  # Variable names
# #         upper.panel = panel.cor,    # Correlation panel
# #         main = "Extinction"
# #         # lower.panel = panel.smooth) # Smoothed regression lines   
# #   )
## # only doing the same species across model types --------------------------
# sp <- Reduce(intersect, list(col_dat$species,cc_dat$species,ext_dat$species)) #overwrite sp
# 
# col_mod_test_abs <- col_dat %>% 
#   filter(species %in% sp) %>% 
#   mod_test_abs()
# cc_mod_test_abs <- cc_dat %>% 
#   filter(species %in% sp) %>% 
#   mod_test_abs()
# ext_mod_test_abs <- ext_dat %>% 
#   filter(species %in% sp) %>% 
#   mod_test_abs()
# # which models have range size as a predictor?
# # colonization
# temp <- unlist(col_mod_test_abs['coef_dredge'])
# names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# col_range_coeff_abs <- temp
# # cc change
# temp <- unlist(cc_mod_test_abs['coef_dredge'])
# names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# cc_range_coeff_abs <- temp
# # extinction
# temp <- unlist(ext_mod_test_abs['coef_dredge'])
# names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
# temp <- temp[str_detect(names(temp), "range")] # keep only those with range size coefficient
# names(temp) <- gsub("\\..*", "", names(temp)) # clean names to only include driver
# ext_range_coeff_abs <- temp
# # range_size coefficient values for the given driver models
# col_range_coeff_abs
# cc_range_coeff_abs
# ext_range_coeff_abs
# 
# # New facet label names for supp variable
# names(col_range_coeff_abs) <- case_when(names(col_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
#                                         names(col_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
#                                         names(col_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
#                                         .default = names(col_range_coeff_abs))
# names(cc_range_coeff_abs) <- case_when(names(cc_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
#                                        names(cc_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
#                                        names(cc_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
#                                        .default = names(cc_range_coeff_abs))
# names(ext_range_coeff_abs) <- case_when(names(ext_range_coeff_abs) == "n_grasshop" ~ "log(n_grasshop)",
#                                         names(ext_range_coeff_abs) == "n_grasshop_previous" ~ "log(n_grasshop_previous)",
#                                         names(ext_range_coeff_abs) == "years_since_last_burn" ~ "log(years_since_last_burn + 1)",
#                                         .default = names(ext_range_coeff_abs))
# 
# f_labels <- bind_rows(mutate(as_tibble(col_range_coeff_abs, rownames = "var"), # add range coefficients
#                              mod = "colonization"),
#                       mutate(as_tibble(cc_range_coeff_abs, rownames = "var"),
#                              mod = "midpoint_change"),
#                       mutate(as_tibble(ext_range_coeff_abs, rownames = "var"),
#                              mod = "extinction")) %>%   
#   rename(range_coeff = value) %>% 
#   mutate(
#     range_coeff = round(range_coeff, 2),
#     label = paste0(range_coeff) #"Range coef.:", 
#   )
# 
# 
# dat <- dat %>% 
#   mutate(sig = case_when(mod == 'colonization' & var %in% names(col_range_coeff_abs) ~ T, # add column to indicate if the driver is significant
#                          mod == 'midpoint_change' & var %in% names(cc_range_coeff_abs) ~ T,
#                          mod == 'extinction' & var %in% names(ext_range_coeff_abs) ~ T),
#   )
# dat %>% 
#   ggplot(aes(x = log(range_size), y = val)) +
#   geom_hline(yintercept = 0, color = "gray80") +
#   geom_point(color = "darkgray") +
#   geom_smooth(method = "lm", aes(color = sig), se = F) +
#   scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
#   facet_grid(var~factor(mod, levels = c("colonization", "midpoint_change", "extinction"), labels = c("Colonization", "Change in cover", "Extinction")), # for whatever reason when I try to order 'mod' this is the only setup where i can get them in order AND with the pretty names
#              scales = "free_y",
#              labeller = labeller(mod = mod.labs, var = var.labs)
#   ) +
#   labs(x = "Range size (log)", y = "Effect size") +
#   geom_text(
#     data    = f_labels,
#     mapping = aes(x = -Inf, y = Inf, label = label),
#     hjust   = -0.25,
#     vjust   = 2
#   )
# ggsave("./figs/5_mod-results-SUP-sameSP.png", height = 11, width = 6)

## # drivers by abiotic or biotic --------------------------------------------
# # abiotic
# dat_for_mod <- drop_na(filter(dat, type == "abiotic"), val) %>% mutate(species = as.factor(species))
# mod <- lmer(log(val) ~ log(range_size) + (1|species), data = dat_for_mod, na.action = "na.fail") # make model. random effect barely changes model
# dredge_out <- MuMIn::dredge(mod) # test with AICc
# dredge_out
# # biotic
# dat_for_mod <- drop_na(filter(dat, type == "biotic"), val) %>% mutate(species = as.factor(species))
# mod <- lmer(log(val) ~ log(range_size) + (1|species), data = dat_for_mod, na.action = "na.fail") # make model. Including random effect does not impact model
# dredge_out <- MuMIn::dredge(mod) # test with AICc
# dredge_out
# # plot
# dat %>% 
#   mutate(type = case_when(type == "abiotic" ~ "Abiotic",
#                           type == "biotic" ~ "Biotic")) %>% 
#   ggplot(aes(x = log(range_size), y = val)) +
#   geom_point() +
#   facet_wrap(~type) +
#   geom_hline(yintercept = 0, color = "gray80") +
#   geom_point(color = "darkgray") +
#   labs(x = "Range size (log)", y = "Effect size")
# ggsave("./figs/5_mod-results-combined-by-type.png", height = 3, width = 5)
# 
