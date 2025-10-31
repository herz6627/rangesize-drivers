
### Herzog 2024

big_func <- function(dat, 
                     param = c("bison + spei + n_grasshop + n_grasshop_previous + years_since_last_burn + fri"), 
                     file_name = NULL, 
                     dir_out = "./output/4_sp-mod-out/", 
                     model_selection = F) {
  
  # make directory 
  dir_out <- paste0(dir_out, file_name, "/")
  
  # create directories if not already present
  dir.create(dir_out)
  dir.create(file.path(dir_out, "sp-colonization"))
  dir.create(file.path(dir_out, "sp-cc-change"))
  dir.create(file.path(dir_out, "sp-extinction"))
  # empty folders if data is already in there
  unlink(file.path(dir_out, "sp-colonization/*"))
  unlink(file.path(dir_out, "sp-cc-change/*"))
  unlink(file.path(dir_out, "sp-extinction/*"))
  
  for (s in sp) { # break data down by species 
    start_time <- Sys.time() # document time  
    sp_dat <- dat %>% # filter to species of interest
      filter(species == s)
    parameters <- str_replace_all(param, "[\r\n]" , "") # model parameters. remove extra lines if needed
    parameters <- gsub(" ", "", parameters, fixed = TRUE) # remove whitespace
    
    ## test if data is missing any factors. If so, drop it from the global model
    if(length(unique(sp_dat$bison)) < 2){ # if lacking multiple levels of bison, remove parameter
      parameters <- str_remove(parameters, "bison\\+")
      parameters <- str_remove(parameters, "\\+bison:n_grasshop")
    }
    if(length(unique(sp_dat$fri)) < 2){ # if lacking multiple levels of fri, remove parameter
      parameters <- str_remove(parameters, "\\+fri")
    }

    
    ##### make model for colonization -----------------------
    col_mod_dat <-
      sp_dat %>% 
      filter(midpoint == 0 & midpoint_next == 0 | midpoint == 0 & midpoint_next != 0)
    parameters_col <- parameters
    if(length(unique(col_mod_dat$bison)) < 2){ # if lacking multiple levels of bison, remove parameter
      parameters_col <- str_remove(parameters_col, "bison\\+")
      parameters_col <- str_remove(parameters_col, "\\+bison:n_grasshop")
    }
    if(length(unique(col_mod_dat$fri)) < 2){ # if lacking multiple levels of fri, remove parameter
      parameters_col <- str_remove(parameters_col, "\\+fri")
    }
    if(length(unique(col_mod_dat$watershed)) < 2){ # if lacking multiple levels of watersheds, remove parameter
      parameters_col <- str_remove(parameters_col, "\\+watershed")
    }
    
    
    if(sum(col_mod_dat$colonization_event) > min_dat){ # arbitrary cutoff of number colonization events to run model
      col_mod <- tryCatch(
        glm(data = col_mod_dat,
            as.formula(paste0("colonization_event ~ ", parameters_col)), 
            na.action = "na.fail", 
            family = "binomial"),
        error = function(e)  {
          message(paste0("error with colonization model for species ",s,". Passing it."))
          NA # if model fails just produce NA
        },
        warning = function(w) {
          message(paste0("warning with colonization model for species ",s,". Passing it."))
          NA # if model fails just produce NA
        }
      )
      # save result
      if (sum(is.na(col_mod)) == 0) { # dont bother saving if it is NA
        coef_out <- col_mod$coefficients
        if (model_selection == T) {
          col_mod <- summary(get.models(MuMIn::dredge(col_mod), 1)[[1]])
          coef_out <- col_mod$coefficients[,1]
          if(length(coef_out)==1) names(coef_out) = "(Intercept)" # if it's an intercept only model, give it a name
        }
        
        n <- nrow(col_mod_dat) #get number of observations used in this model
        names(n) <- "n" # give it a name
        out <- append(coef_out, n) # add on to our model output
        
        saveRDS(out, paste0(dir_out ,"sp-colonization/4_", s, "-colonization", ".RDS"))
      }
    }
    print("done with colonization")
    
    ##### make model for cover changes ----------------------
    cc_mod_dat <-
      sp_dat %>% 
      filter(midpoint != 0 & midpoint_next != 0)
    parameters_cc <- parameters
    if(length(unique(cc_mod_dat$bison)) < 2){ # if lacking multiple levels of bison, remove parameter
      parameters_cc <- str_remove(parameters_cc, "bison\\+")
      parameters_cc <- str_remove(parameters_cc, "\\+bison:n_grasshop")
    }
    if(length(unique(cc_mod_dat$fri)) < 2){ # if lacking multiple levels of fri, remove parameter
      parameters_cc <- str_remove(parameters_cc, "\\+fri")
    }
    if(length(unique(cc_mod_dat$watershed)) < 2){ # if lacking multiple levels of fri, remove parameter
      parameters_cc <- str_remove(parameters_cc, "\\+watershed") # if present as fixed
      parameters_cc <- str_remove(parameters_cc, "\\+\\(1\\|watershed\\)") # if present as random
    }
    
    if(nrow(cc_mod_dat) > min_dat & length(unique(cc_mod_dat$midpoint)) > 1){ # need enough data and more than 1 coverclass
      if (str_detect(parameters_cc, "\\(1\\|")) { # if random effects are present
        cc_mod <- glmer(data = cc_mod_dat,
                        as.formula(paste0("(midpoint_next/midpoint) ~ ", parameters_cc)), 
                        , na.action = "na.fail")
      } else {
        if(!(str_detect(parameters_cc, "\\(1\\|"))) { # if random effects are not present
          cc_mod <- glm(data = cc_mod_dat,
                        as.formula(paste0("(midpoint_next/midpoint) ~ ", parameters_cc)), 
                        , na.action = "na.fail")
        }
      }
      
      if (sum(is.na(cc_mod)) == 0) { # dont bother saving if it is NA. have to use sum for when the model objects have length > 1
        coef_out <- cc_mod$coefficients
        if (model_selection == T) {
          cc_mod <- summary(get.models(MuMIn::dredge(cc_mod), 1)[[1]])
          coef_out <- cc_mod$coefficients[,1]
          if(length(coef_out)==1) names(coef_out) = "(Intercept)" # if it's an intercept only model, give it a name
        }
        
        n <- nrow(cc_mod_dat) #get number of observations used in this model
        names(n) <- "n" # give it a name
        out <- append(coef_out, n) # add on to our model output
        
        saveRDS(out, paste0(dir_out ,"sp-cc-change/4_", s, "-cc-change", ".RDS"))
      }
    }
    print("done with cc changes")

    
    ##### make model for extinction -------------------------
    # filter data
    ext_mod_dat <-
      sp_dat %>% 
      filter(midpoint != 0 & midpoint_next != 0 | midpoint != 0 & midpoint_next == 0)
    # clean up parameters
    parameters_ext <- parameters # modify parameters if needed for this model
    if(length(unique(ext_mod_dat$bison)) < 2){ # if lacking multiple levels of bison, remove parameter
      parameters_ext <- str_remove(parameters_ext, "bison\\+")
      parameters_ext <- str_remove(parameters_ext, "\\+bison:n_grasshop")
    }
    if(length(unique(ext_mod_dat$fri)) < 2){ # if lacking multiple levels of fri, remove parameter
      parameters_ext <- str_remove(parameters_ext, "fri\\+")
    }
    if(length(unique(ext_mod_dat$watershed)) < 2){ # if lacking multiple levels of fri, remove parameter
      parameters_ext <- str_remove(parameters_ext, "\\+watershed")
    }
    
    # run model
    if(sum(ext_mod_dat$extinction_event) > min_dat){ # only run if we have enough extinction events
      ext_mod <- tryCatch(
        glm(data = ext_mod_dat,
            as.formula(paste0("extinction_event ~ ", parameters_ext))
            , na.action = "na.fail", 
            family = "binomial"
        ),
        error = function(e)  {
          message(paste0("error with extinction model for species ",s,". Passing it."))
          NA # if model fails just produce NA
        },
        warning = function(w) {
          message(paste0("warning with extinction model for species ",s,". Passing it."))
          NA # if model fails just produce NA
        }
      )

      # save result
      if (sum(is.na(ext_mod)) == 0) { # dont bother saving if it is NA
        coef_out <- ext_mod$coefficients
        if (model_selection == T) {
          ext_mod <- summary(get.models(MuMIn::dredge(ext_mod), 1)[[1]])
          coef_out <- ext_mod$coefficients[,1]
          if(length(coef_out)==1) names(coef_out) = "(Intercept)" # if it's an intercept only model, give it a name
        }
        
        n <- nrow(ext_mod_dat) #get number of observations used in this model
        names(n) <- "n" # give it a name
        out <- append(coef_out, n) # add on to our model output
        
        saveRDS(out, paste0(dir_out ,"sp-extinction/4_", s, "-extinction", ".RDS"))  
      }
      
    }
    print("done with extinction")
    print(paste0("done with ", s))
  }
}
