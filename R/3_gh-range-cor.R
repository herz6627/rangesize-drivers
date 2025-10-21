# Herzog, Fall 2024

# Load library ------------------------------------------------------------
library(tidyverse)

# load files --------------------------------------------------------------
range_sizes <- read_csv("./output/2_range-size-df.csv") %>% # range sizes
  rename(species = sp,
         range_est = area) %>% 
  select(!range_size_class) # extra info

# grasshopper network ---
# welti et al. 2018 data
path <- "./data/doi_10_5061_dryad_s6j1822__v20190107/Welti et al. 2018_Online data_Plant-grasshopper networks.xlsx"
sheetnames <- readxl::excel_sheets(path)
sheetnames <- sheetnames[-1] # remove info sheet
mylist <- lapply(sheetnames, readxl::read_excel, path = path) # read in each sheet
names(mylist) <- sheetnames # name the dataframes

# combine matrices into one large matrix with all grasshopper sp and plant sp
sum <- matrix(0, 0, 0)
for (i in seq_along(mylist)) {
  Mat1 <- sum
  Mat2 <- mylist[[i]] %>% 
    column_to_rownames("...1") %>% 
    as.matrix()
  # Get the row and column names 
  rn1 <- rownames(Mat1)
  rn2 <- rownames(Mat2) # since its a tibble gotta extract column
  cn1 <- colnames(Mat1)
  cn2 <- colnames(Mat2)
  # Construct row and column names for the sum matrix
  rnsum <- unique(c(rn1, rn2))
  cnsum <- unique(c(cn1, cn2))
  # Make the matrix of zeros
  sum <- matrix(0, length(rnsum), length(cnsum),
                dimnames = list(rnsum, cnsum))
  # Put all indices of each matrix into a matrix
  # with column 1 being the row name, column 2 being the 
  # column name, and add the results into the sum
  ind <- cbind(rn1[row(Mat1)], cn1[col(Mat1)])
  sum[ind] <- sum[ind] + Mat1[ind]
  ind <- cbind(rn2[row(Mat2)], cn2[col(Mat2)])
  sum[ind] <- sum[ind] + Mat2[ind]
} 

# replace values with 1/0 as values in matrix are now meaningless
sum[sum > 0] <- 1 

# filter to grasshopper species with at least 10 grasshopper samples
# based off Table S1 in Supp. material of Welti 2018
sum <- sum[,!colnames(sum) %in% c("A conspersa", "B auriventris", "B magna", "M differentialis", "M sanguinipus", "P brunneri")]

# species found to be eaten by grasshoppers:
pl_gh <-  str_to_sentence(rownames(sum))
# how many grasshopper species eat a plant species
pl_gh_links <- rowSums(sum)
names(pl_gh_links) <- str_to_sentence(names(pl_gh_links)) # format species names to match



# Begin -------------------------------------------------------------------

## correlation between links and range size --------------------------------
pl_gh_links <- pl_gh_links %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(species = 1, gh_links = 2) %>% 
  as_tibble() %>% 
  separate(species, c("genus", "epithet"), remove = F) %>% 
  mutate(id_level = case_when(
    epithet == "spp" ~ "genus",
    epithet != "spp" ~ "species"
  )) %>% 
  left_join(range_sizes) %>% 
  filter(gh_links != 0) # dont care about the sp that gh dont eat
temp1 <- range_sizes %>% #filter data to just species and range size
  separate(species, c("genus", "epithet"), remove = F) %>% 
  distinct() %>% 
  filter(genus %in% unique(pl_gh_links$genus)) %>%  # filter genera to those found in grasshoppers
  left_join(select(pl_gh_links, species, gh_links, id_level))
temp2 <- temp1 %>% # add link data for genus level observations
  filter(is.na(gh_links)) %>% # filter to species missing link data
  select(-c(gh_links, id_level)) %>% 
  left_join(select(
    filter(pl_gh_links, epithet == "spp") # filter to only those observations at the genus level
    , genus, gh_links, id_level), by = "genus") %>%  # NOTE there may be discrepancies due to multiple genus matches
  filter(!is.na(gh_links)) # filter to species *not* missing link data
range_gh_dat <- temp1 %>% # combine data
  filter(!is.na(gh_links)) %>% # filter to species *not* missing link data
  bind_rows(temp2) %>%  # add to first table
  filter(epithet != "x") # remove hybrids


# model selection
mod <- lm(log(gh_links+1) ~ log(range_est), data = range_gh_dat, na.action = "na.fail") # make model. random effect barely changes model
dredge_out <- MuMIn::dredge(mod) # test with AICc
dredge_out
DHARMa::simulateResiduals(fittedModel = mod, plot = T) # check model fit
summary(MuMIn::get.models(dredge_out, 1)[[1]]) # best model

# with summary stats
MuMIn::dredge(mod, m.lim = c(NA, 1), extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
    
  })
)


# graph -------------------------------
theme_pretty <- theme_set(
  theme_light()
)
theme_pretty <- theme_update(
  strip.background = element_rect(
    color="black", fill="gray80", size=1.5, linetype="blank"),
  strip.text = element_text(colour = 'black')
)


range_gh_dat %>% 
  ggplot(aes(x = range_est/(10^7), y = log(gh_links))) + # transform range sizes?
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  theme_pretty +
  # labs(x = "Plant range size, log", y = "Number consumer grasshopper species, log")
  labs(x = expression(paste("Plant range size * ", 10^7, " ", km^2)), y = "Number consumer grasshopper species, log")
ggsave("./figs/3_gh-range-correlation.png", units = "in", width = 4, height = 4) # save


# save data
write_csv(range_gh_dat, "./output/3_range-gh-dat.csv")
