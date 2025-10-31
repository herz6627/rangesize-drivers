# extract coefficient names and their estimates for each species
mod_list_to_df <- function(mod){ # helper function
  # extract coefficient names and their estimates for each species
  out <- data.frame()
  for (i in 1:length(mod)) {
    if (class(mod[[i]])[1] %in% c("lmerMod", "glmerMod")) {
      temp <- fixef(mod[[i]])
    } else {
      if (class(mod[[i]])[1] == "glm") {
        temp <- coefficients(mod[[i]]) # get coefficients
      } else {
        if(class(mod[[i]])[1] == "numeric") { # when model coeffs have already been extracted into list
          temp <- mod[[i]]
        } else {stop("unknown type of model")}
      }
    }
    # print(i)
    out <- bind_rows(out,temp)
  }
  # out <- as.data.frame(out[-1,]) # drop first row of missing values
  out[,'species'] <- names(mod) # add column with species names
  return(out)
}

mod_test <- function(dat, include_interaction = F, include_lifespan = T, plot.title = NULL, include_weights = F, transform_coeffs){ # helper function when including functional groups
  out1 <- list() # list to fill with full dredge output
  out2 <- list() # list to fill with just coefficients from dredge output
  out3 <- list() # list to fill with top model
  
  params <- (ncol(dat)-3) # number of parameters # ! this value (3) includes a column of n (+species+range_size)
  message(paste("detected ", params, " parameters"))
  if (transform_coeffs == T) {require(glmmTMB)}
  
  for (col in 2:params) { # iterate over columns (drivers)
    temp <- dat %>% 
      left_join(sp_info) %>% 
      mutate(across(lifespan:photo, as.factor)) %>% 
      select(all_of(col), range_size, lifespan, growthform, n) %>% 
      drop_na()
    
    if (nrow(temp) < 5) {print (paste('less than 5 coefficients found, skipping',  colnames(temp)[1])) ; next}
    
    message(paste0("\n",colnames(temp)[1], " has n = ", nrow(temp) ," observations for the model")) # helpful info
    colnames(temp)[1] <-  sub("log\\(", "", colnames(temp)[1]) # rename columns to remove any transformation notation (e.g. log())
    colnames(temp)[1] <-  sub(" \\+ 1\\)", "", colnames(temp)[1]) 
    colnames(temp)[1] <-  sub("\\)", "", colnames(temp)[1]) 
    colnames(temp)[1] <-  sub("\\:", ".", colnames(temp)[1]) # replace ':' for interaction effects
    
    # make model
    if(include_interaction == F & include_lifespan == T) form <- paste0(names(temp)[1], "~range_size + lifespan") #mod <- lm(as.formula(paste0(names(temp)[1], "~range_size + lifespan")), data = temp, na.action = "na.fail") # + growthform
    if(include_interaction == T & include_lifespan == T) form <- paste0(names(temp)[1], "~range_size*lifespan") # mod <- lm(as.formula(paste0(names(temp)[1], "~range_size*lifespan")), data = temp, na.action = "na.fail") #*growthform
    if(include_lifespan == F) {
       form <- as.formula(paste0(names(temp)[1], "~range_size")) # mod <- lm(as.formula(paste0(names(temp)[1], "~range_size")), data = temp, na.action = "na.fail")
      print("ignoring 'include_interaction' parameter")
    }
    
    if (length(unique(temp$lifespan)) < 2 & include_lifespan == T) {
      form <- as.formula(paste0(names(temp)[1], "~range_size"))
        print("not enough lifespan levels, dropping parameter from model")
        
    }
    
    
    if(include_weights == T) mod <- lm(as.formula(form), data = temp, na.action = "na.fail", weights = n)
    if(include_weights == F) mod <- lm(as.formula(form), data = temp, na.action = "na.fail")

    # if(include_weights == T & transform_coeffs == F) mod <- lm(as.formula(form), data = temp, na.action = "na.fail", weights = n)
    # if(include_weights == F & transform_coeffs == F) mod <- lm(as.formula(form), data = temp, na.action = "na.fail")
    # 
    # if(include_weights == T & transform_coeffs == T) mod <- glmmTMB::glmmTMB(as.formula(form), family = ziGamma, data = temp, na.action = "na.fail", weights = n)
    # if(include_weights == F & transform_coeffs == T) mod <- glmmTMB::glmmTMB(as.formula(form), family = ziGamma, data = temp, na.action = "na.fail")
    #
    # residuals
    # simulationOutput <- DHARMa::simulateResiduals(fittedModel = mod, plot = F)
    # plot(simulationOutput, title = paste0(plot.title, " ", colnames(temp)[1]))
    # hist(simulationOutput)
    
    # plot(mod, which = 1)
    # title(main = paste0(plot.title, " ", colnames(temp)[1]))
    
    # histogram of residuals
    print(
      ggplot(data = temp, aes(x = mod$residuals)) +
      geom_histogram(fill = 'steelblue', color = 'black') +
      labs(title = paste('Histogram of Residuals\n', plot.title, " ", colnames(temp)[1]), x = 'Residuals', y = 'Frequency')
    )
    
    
    # save output
    temp_name = names(temp)[1] # name object in list as driver
    
    out1[[col]] <- MuMIn::dredge(mod) # test with AICc
    names(out1)[col] <- temp_name
    
    top_mod = MuMIn::get.models(out1[[col]], subset = 1)[[1]] # top model
    
    out2[[col]] <- top_mod$coefficients # top model coefficients
    names(out2)[col] <- temp_name
    
    out3[[col]] <- temp %>%
      rename(val = 1) %>% # rename first column to show it is just the driver response
      mutate(pred = predict(top_mod),
             var = temp_name) # driver name column
    names(out3)[col] <- temp_name
    
    

  }
  out <- list(out1, out2, out3
              )
  names(out) <- c("full_dredge_out", "coef_dredge", "predicted")
  return(out)
}

# extract which models have range as a significant predictor for drivers
get_sig_coeffs <- function(col_test_results, cc_test_results, ext_test_results){
  # which models have range size as a predictor?
  out <- list()
  # colonization
  temp <- unlist(col_test_results['coef_dredge'])
  names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
  temp <- temp[str_detect(names(temp), "(?<!:)range_size(?!:)")] # keep only those with range size coefficient # when log transformed: "(?<!:)log\\(range_size\\)(?!:)"
  names(temp) <- gsub("\\.range_size", "", names(temp)) # clean names to only include driver
  col_range_coeff_abs <- temp
  # cc change
  temp <- unlist(cc_test_results['coef_dredge'])
  names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
  temp <- temp[str_detect(names(temp), "(?<!:)range_size(?!:)")] # keep only those with range size coefficient
  names(temp) <- gsub("\\.range_size", "", names(temp)) # clean names to only include driver
  cc_range_coeff_abs <- temp
  # extinction
  temp <- unlist(ext_test_results['coef_dredge'])
  names(temp) <- str_remove(names(temp), "coef_dredge.") # remove extra strings
  temp <- temp[str_detect(names(temp), "(?<!:)range_size(?!:)")] # keep only those with range size coefficient, no interaction (searches for ':' and drops log(range_size) if present)
  names(temp) <- gsub("\\.range_size", "", names(temp)) # clean names to only include driver
  ext_range_coeff_abs <- temp
  # range_size coefficient values for the given driver models
  out[["col"]] <- col_range_coeff_abs
  out[["cc"]] <- cc_range_coeff_abs
  out[["ext"]] <- ext_range_coeff_abs
  
  return(out)
}



# make consistent figures ---------------

# coefficient value figures ----------------
plot_coefficients <- function(
    dat, 
    filterXL = filterXL, # whether to remove largest range sizes (>2)
    coeff_type = c("all", "primary", "secondary") # plot all or a subset of coefficients
) {
  coeff_type <- match.arg(coeff_type)
  
  # Map variable names
  temp <- dat %>%
    mutate(var = case_when(
      var == 'bison1' ~ "Bison",
      var == 'fri' ~ "Burn~interval",
      var == 'n_grasshop' ~ "Grasshopper",
      var == 'n_grasshop_previous' ~ "Grasshopper[t-1]",
      var == 'years_since_last_burn' ~ "Time~since~fire",
      var == 'spei' ~ "Climate",
      var == "bison1:n_grasshop" ~ "Bison:GH",
      var == 'bison1_mod' ~ "Bison_mod",
      var == 'n_grasshop_mod' ~ "GH_mod"
    ),
    var = fct_relevel(var,c("Climate", "Time~since~fire", "Burn~interval", "Bison", "Grasshopper", "Grasshopper[t-1]", "Bison:GH", "Bison_mod", "GH_mod"))) %>% 
    {if (filterXL) filter(., range_size < 2) else .}
  # 
  # dat_text <- data.frame(
  #   label = c("A", "B", "C"),
  #   mod   = c("colonization", "midpoint_change", "extinction")
  # )
  
  # Select variables based on coeff_type
  vars_to_plot <- switch(
    coeff_type,
    all = unique(temp$var),
    primary = c("Burn~interval","Grasshopper[t-1]","Climate","Time~since~fire","Bison","Grasshopper"),
    secondary = c(
      
      "Bison_mod","GH_mod","Bison:GH"
      )
  )
  
  temp_plot <- temp %>% filter(var %in% vars_to_plot)
  
  p <- ggplot(temp_plot, aes(x = var, y = val)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_jitter(color = "gray", size = 0.4, alpha = 0.9) +
    geom_violin(alpha = 0.6, fill = "gray") +
    labs(x = "Driver", y = "Effect size") +
    scale_x_discrete("Driver", labels = parse(text = as.character(sort(unique(temp_plot$var))))) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, lineheight = 0.5)
    ) +
    facet_wrap(~factor(mod, levels = c("colonization", "midpoint_change", "extinction"),
                       labels = c("Colonization", "Change in cover", "Extirpation"))) 
    # label on right
    # geom_text(
    #   data = dat_text,
    #   mapping = aes(x = Inf, y = Inf, label = label),
    #   hjust = 2, vjust = 2, label.size = 0.5, fontface = "bold"
    # )
    # label on left
    # geom_text(
    #   data = dat_text,
    #   mapping = aes(x = -Inf, y = Inf, label = label),
    #   hjust = -0.5, vjust = 2, label.size = 0.5, fontface = "bold"
    # )
  
  return(p)
}

# linear trend figures ------------
make_fig <- function(dat_df = dat, pred_df = pred_dat, sig_coeffs_list,
                     height = 7, width = 4, fld = "./figs/5_mod-results.png",
                     filter.sp = FALSE, filter.sp.vec = NULL, filter.var.vec = NULL,
                     val.transform = TRUE) {
  
  d <- sig_coeffs_list
  
  # Lookup tables for nicer labels
  var_lookup <- c(
    bison1 = "Bison",
    fri = "Burn~interval",
    spei = "Climate",
    n_grasshop = "Grasshopper",
    n_grasshop_previous = "Grasshopper[t-1]",
    years_since_last_burn = "Time~since~fire",
    `bison1.n_grasshop` = "Bison:GH",
    `bison1:n_grasshop` = "Bison:GH",
    bison1_mod = "Bison_mod",
    n_grasshop_mod = "GH_mod"
  )
  
  mod_lookup <- c(
    colonization = "Colonization",
    midpoint_change = "Change~'in'~cover",
    extinction = "Extirpation"
  )
  
  # Combine coefficients for labeling
  f_labels <- bind_rows(
    as_tibble(d[[1]], rownames = "var") %>% mutate(mod = "Colonization"),
    as_tibble(d[[2]], rownames = "var") %>% mutate(mod = "Change~'in'~cover"),
    as_tibble(d[[3]], rownames = "var") %>% mutate(mod = "Extirpation")
  ) %>%
    filter(!var %in% filter.var.vec) %>%
    rename(range_coeff = value) %>%
    mutate(
      range_coeff = round(range_coeff, 2),
      label = as.character(range_coeff),
      var = recode(var, !!!var_lookup),
      var = fct_relevel(var,c("Climate", "Time~since~fire", "Burn~interval", "Bison", "Grasshopper", "Grasshopper[t-1]", "Bison:GH", "Bison_mod", "GH_mod"))
      
    )
  
  # Helper to apply filters and recode
  prepare_df <- function(df, filter_species = filter.sp) {
    # Only filter if the species column exists and filtering is requested
    if(filter_species && "species" %in% names(df)) {
      df <- df %>% filter(species %in% filter.sp.vec)
    }
    
    df %>%
      filter(!var %in% filter.var.vec) %>%
      mutate(
        mod = recode(mod, !!!mod_lookup), # !!! unpacks the named vector so each named element becomes its own argument for recode()
        var = recode(var, !!!var_lookup),
        var = fct_relevel(var,c("Climate", "Time~since~fire", "Burn~interval", "Bison", "Grasshopper", "Grasshopper[t-1]", "Bison:GH", "Bison_mod", "GH_mod")),
        sig = mapply(
          function(m, v) v %in% f_labels$var[f_labels$mod == m],
          mod, var
        )
      )
  }
  
  temp_dat <- prepare_df(dat_df) %>% 
    mutate(var = fct_drop(var)) # drop any extra factor levels (i.e. (Intercept))
  temp_pred_df <- prepare_df(pred_df)
  
  # Panel annotation setup
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dat_text <- data.frame(
    label = LETTERS702[1:(3*length(unique(temp_dat$var)))],
    mod   = rep(mod_lookup, times = length(unique(temp_dat$var))),
    var   = rep(levels(temp_dat$var), each = 3)
  )
  
  # Plot
  p <- ggplot(temp_dat, aes(x = range_size, y = val)) +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_point(color = "darkgray", alpha = 0.75) +
    geom_line(data = filter(temp_pred_df, lifespan == "p"), # only plotting line for perennial species. This factor has been acounted for in our models
              aes(y = pred, color = sig), linewidth = 1) +
    scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
    facet_grid(factor(var, levels = levels(temp_dat$var)) ~ factor(mod, levels = mod_lookup), # have to force the facet order
               labeller = label_parsed) +
    geom_text(data = f_labels, aes(x = Inf, y = -Inf, label = label),
              hjust = 1.25, vjust = -0.5) +
    labs(x = expression(paste("Range size * ", 10^7, " ", km^2)),
         y = "Effect size") +
    geom_text(data = dat_text, aes(x = Inf, y = Inf, label = label),
              hjust = 1.5, vjust = 1.5, label.size = 0.5, fontface = "bold")
  
  print(p)
  ggsave(fld, height = height, width = width)
  message(paste0("image saved to ", fld))
}

# OLD
# make_fig <- function(dat_df = dat, pred_df = pred_dat, sig_coeffs_list, height = 7, width = 4, fld = "./figs/5_mod-results.png", filter.sp = F, filter.sp.vec = NULL, filter.var.vec = NULL, val.transform = T){
#   
#   d <- sig_coeffs_list
#   
#   # New facet label names for supp variable
#   f_labels <- bind_rows(mutate(as_tibble(d[[1]], rownames = "var"), # add range coefficients
#                                mod = "Colonization"),
#                         mutate(as_tibble(d[[2]], rownames = "var"),
#                                mod = "Change~'in'~cover"),
#                         mutate(as_tibble(d[[3]], rownames = "var"),
#                                mod = "Extirpation")) %>%
#     filter(!var %in% filter.var.vec) %>%  # remove any unwanted drivers
#     rename(range_coeff = value) %>%
#     mutate(
#       range_coeff = round(range_coeff, 2),
#       label = paste0(range_coeff),  #"Range coef.:",
#       var = case_when( # make pretty labels
#         var == "bison1" ~ "Bison",
#         var == "fri" ~ "FRI",
#         var == "n_grasshop" ~ "GH",
#         var == "n_grasshop_previous" ~ "GH[t-1]",
#         var == "years_since_last_burn" ~ "TSF",
#         var == "spei" ~ "SPEI",
#         var == "bison1.n_grasshop" ~ "Bison:GH",
#         
#         var == 'bison1_mod' ~ "Bison_mod",
#         var == 'n_grasshop_mod' ~ "GH_mod"
#       ))
#   
#   temp_dat <- dat_df %>% 
#     {if(filter.sp == T) {dplyr::filter(., species %in% filter.sp.vec)} else .} %>% # remove any species if specified
#     dplyr::filter(!var %in% filter.var.vec) %>% # remove any drivers if specified
#     mutate(
#       mod = case_when( # make pretty labels
#         mod == "colonization" ~ "Colonization",
#         mod == "midpoint_change" ~ "Change~'in'~cover",
#         mod == "extinction" ~ "Extirpation"
#       ),
#       var = case_when( # make pretty labels
#         var == "bison1" ~ "Bison",
#         var == "fri" ~ "FRI",
#         var == "n_grasshop" ~ "GH",
#         var == "n_grasshop_previous" ~ "GH[t-1]",
#         var == "years_since_last_burn" ~ "TSF",
#         var == "spei" ~ "SPEI",
#         var == "bison1:n_grasshop" ~ "Bison:GH",
#         
#         var == 'bison1_mod' ~ "Bison_mod",
#         var == 'n_grasshop_mod' ~ "GH_mod"
#       ),
#       sig = case_when(mod == 'Colonization' & var %in% filter(f_labels, mod == "Colonization")$var ~ T, # add column to indicate if the driver is significant
#                       mod == "Change~'in'~cover" & var %in% filter(f_labels, mod == "Change~'in'~cover")$var ~ T,
#                       mod == 'Extirpation' & var %in% filter(f_labels, mod == "Extirpation")$var ~ T)
#         
#     ) 
#   
#   temp_pred_df = pred_df %>% 
#     # {if(filter.sp == T) {dplyr::filter(., species %in% filter.sp.vec)} else .} %>% # remove any species if specified
#     filter(!var %in% filter.var.vec) %>% # remove any drivers if specified
#     mutate(
#       mod = case_when( # make pretty labels
#         mod == "colonization" ~ "Colonization",
#         mod == "midpoint_change" ~ "Change~'in'~cover",
#         mod == "extinction" ~ "Extirpation"
#       ),
#       var = case_when( # make pretty labels
#         var == "bison1" ~ "Bison",
#         var == "fri" ~ "FRI",
#         var == "n_grasshop" ~ "GH",
#         var == "n_grasshop_previous" ~ "GH[t-1]",
#         var == "years_since_last_burn" ~ "TSF",
#         var == "spei" ~ "SPEI",
#         var == "bison1.n_grasshop" ~ "Bison:GH",
#         
#         var == 'bison1_mod' ~ "Bison_mod",
#         var == 'n_grasshop_mod' ~ "GH_mod"
#       ),
#       sig = case_when(mod == 'Colonization' & var %in% filter(f_labels, mod == "Colonization")$var ~ T, # add column to indicate if the driver is significant
#                       mod == "Change~'in'~cover" & var %in% filter(f_labels, mod == "Change~'in'~cover")$var ~ T,
#                       mod == 'Extirpation' & var %in% filter(f_labels, mod == "Extirpation")$var ~ T)
#       
#     ) 
#     
#   # add panel annotations
#   LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS))) # need more than 26 letters sometimes
#   dat_text <- data.frame( # add panel labels
#     label = LETTERS702[1:(3*length(unique(temp_dat$var)))],
#     mod   = c("Colonization", "Change~'in'~cover", "Extirpation"),
#     var = rep(sort(unique(temp_dat$var)), each=3)
#   )
#   
#   # plot
#   p <- 
#     temp_dat %>% 
#     ggplot(aes(x = range_size, y = val)) +
#     geom_hline(yintercept = 0, color = "gray80") +
#     geom_point(color = "darkgray", alpha = 0.75) +
#     # geom_smooth(method = "lm", aes(color = sig), se = F) +
#     
#     geom_line(data = filter(temp_pred_df, lifespan == "p"), aes(y = pred, color = sig), linewidth = 1) + # NOTE: we are only including perennial values for the figures
#     
#     
#     scale_color_manual(values = c("TRUE" = "darkslateblue"), na.value = NA, guide = "none") +
#     facet_grid(var~factor(mod, levels = c("Colonization", "Change~'in'~cover", "Extirpation")), # for whatever reason when I try to order 'mod' this is the only setup where i can get them in order AND with the pretty names
#                # scales = "free_y",
#                labeller = label_parsed) +
#     geom_text(
#       data    = f_labels,
#       # # top left
#       # mapping = aes(x = -Inf, y = Inf, label = label),
#       # hjust   = -0.25,
#       # vjust   = 2
#       # bottom right
#       mapping = aes(x = Inf, y = -Inf, label = label),
#       hjust   = 1.25,
#       vjust   = -0.5
#     ) +
#     labs(x =expression(paste("Range size * ", 10^7, " ", km^2)), y = "Effect size") +
#     geom_text( # add panel annotations
#       data    = dat_text,
#       mapping = aes(x = Inf, y = Inf, label = label),
#       hjust   = 1.5,
#       vjust   = 1.5,
#       label.size = 0.5,
#       fontface="bold"
#     )
#   # save
#   print(p)
#   ggsave(fld, height = height, width = width)
#   # finish
#   message(paste0("image saved to", fld))
# }

