#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Systematic review on herbivory diversity in Arctic tundra
#                      Meta-regression models
#     Jonas Trepel (jonas.trepel@gmail.com / jonas.trepel@bio.au.dk)
#     Erick J Lundgren (erick.lundgren@gmail.com / ejlundgren@bio.au.dk)
#                        28/29-April-2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# what happens in this script: 
# 0. get an overview of the responses and data distribution 
# 1. intercept only models 
# 1.1 create model guide
# 1.2. loop models
# 2. plot

# clean the environment first
rm(list=ls())

# libraries --------------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)
library(harrypotter) # colour palettes for consistency with systematic map
library(ggnewscale)  # to add several scales to ggplot
library(gridExtra)
library(scales)
library(ggpubr)
library(data.table)
library(tictoc) 
library(metafor)
library(stringr)
library(clubSandwich)  # needed to do the sensitivity analysis
library(flextable)    # to make beautiful tables :)
library(orchaRd)  # meta-analysis visualization package
                  # needed to calculate measures of heterogeneity for multilevel meta-analysis, marginal R2 etc
  # devtools::install_github("daniel1noble/orchaRd", force = TRUE)


# load functions ---------------------------------------------------------------
# this function (written by Erick) allows to make predictions with the models 
# https://stackoverflow.com/questions/63554740/predict-rma-onto-complex-new-data-in-metafor-polynomials-and-factor-levels

rma_predictions <- function(m, newgrid){
  if(!is.data.frame(newgrid)){errorCondition("ERROR newgrid must be a data frame")}
  
  # create the new model matrix 
  if(!all(unlist(lapply(names(newgrid), # lapply through names of newgrid to check that they're in the model formula
                        function(x) grepl(pattern=x, x = as.character(m$formula.mods)[-1]))))){
    errorCondition("ERROR: variables in newgrid are not in model formula")
  }
  
  predgrid <- (model.matrix(m$formula.mods, data = newgrid))
    predgrid
  
  if(any(grepl("intercept", colnames(predgrid), 
               ignore.case = TRUE))){
    #if intercept is present, remove it?
    predgrid <- predgrid[, -1]
  }
  
  # predict onto the new model matrix
  pred.out <- as.data.frame(predict(m, newmods=predgrid))
  
  # attach predictions to variables for plotting
  final.pred <- cbind(newgrid, pred.out)
    final.pred
}

# defining the parameters for our graphs
theme_systrev <- function(){ #create a new theme function for the style of graphs
  theme_bw() +                  #use a predefined theme as a base
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10, color="black"), 
          axis.title = element_text(size = 14, color="black"),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.title = element_text(size = 18, vjust = 1, hjust = 0, color="black"),
          legend.text = element_text(size = 12, color="black"),          
          legend.title = element_text(size = 10, color="black"),
          legend.position = "bottom",
          legend.key.width = unit(2.5, "cm"),
          legend.key.height = unit(0.8, "lines"),
          legend.key.size = unit(0.8, "lines"),
          legend.box = "vertical", 
          legend.margin = margin())
}


# load data --------------------------------------------------------------------
# the data loaded here was built with script 4_calculate_effect_sizes
# and includes the coding database with the corresponding effect sizes
# file: data_with_effect_sizes.csv

# set the working directory to the folders where the data files are stored
setwd("./data")

# load the dataset
dt <- fread("effect_sizes/data_with_effect_sizes.csv")


# clean up a bit (note that dataset dt_f includes more variables than in script 5)
dt_f <- dt %>% filter(!is.na(yi_smd)) %>%   # remove missing values (287 studies) 
            # remove also studies for which there is no contrast (133 studies)
            filter(!change.long_f == "no_contrast") %>%  
            # remove extreme values (36 studies)
            filter(yi_smd > quantile(dt$yi_smd, .005, na.rm = TRUE) & 
                   yi_smd < quantile(dt$yi_smd, .995, na.rm = TRUE)) %>% 
            mutate(article_ID = as.character(article_ID), # make sure article_ID is character
                   # define study length
                   study_length = suppressWarnings(as.numeric(year_end))- 
                                   suppressWarnings(as.numeric(year_start))) %>% 
            # select relevant columns
            select("article_ID", "study_ID", 
                   "herbivore_fgr_ID_higher_diversity", "herbivore_fgr_ID_lower_diversity",
                   "herbivore_fgr_nr_high", "herbivore_fgr_nr_low", "change_f", "change.long_f",
                   "var_class", "MA.value", "outcome_type", 
                   # include potential moderators
                   "overall_criterion", "error_type", "spatial_resolution", "year",
                   "sample_size_higher_diversity", "sample_size_lower_diversity", "study_length",
                   # include environmental context
                   "permafrost_D", "habitat_type_D", "extent_of_recent_change", "recent_greening", 
                   "recent_warming", "productivity", "soil_type_D", "bioclimatic_zone", 
                   "distance_from_coast", "growing_season", "precipitation", "temperature",
                   "distance_to_treeline", "elevation_DEM", "soil_chemistry", "soil_texture",
                   "soil_moisture", `soil_type `, "permafrost", "habitat_type",
                   "yi_smd", "vi_smd") %>% 
            rename(soil_type = `soil_type `) %>% 
            # create a variable for "exclusion studies" (lower diversity = zero)
            mutate(exclusion = ifelse(herbivore_fgr_ID_lower_diversity == "zero", "complete", "partial")) %>% 
            # calculate numerical variable for the change in herbivore group richness
            mutate(herb_fgr_change = herbivore_fgr_nr_high - herbivore_fgr_nr_low)

# count the number of articles/studies for each response variable
resp_var_f <- dt_f %>% group_by(var_class, MA.value, article_ID) %>% 
                summarize(ns = n()) %>%
                  group_by(var_class, MA.value) %>% 
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
              arrange(desc(nr_articles))



# 1. build the models including one moderator at a time -----------------------------------------------
# (these are 475 models! can skip this step and go directly to step 2)

# here we will use studies reporting variables measured by at least 10 articles
length(subset(resp_var_f, nr_articles >= 10)$MA.value) # 25 variables have 10 articles or more

# filter the data to include only responses with more than 10 articles  
dt_10_f <- dt_f %>% filter(MA.value %in% subset(resp_var_f, nr_articles >= 10)$MA.value) # 2231 studies

# check correlations between environmental variables that 
# we want to consider as moderators
moderators <- dt_10_f %>% 
                # we make our bioclimatic_zone variable numeric
                mutate(bioclim_z = as.numeric(unclass(as.factor(bioclimatic_zone)))) %>% 
                select("elevation_DEM", "distance_to_treeline", "distance_from_coast", 
                       "bioclim_z", "temperature", "precipitation", "growing_season", 
                       "productivity", "recent_warming", "recent_greening",  
                       "extent_of_recent_change","permafrost_D")
summary(moderators)

ggcorrplot(cor(moderators, use = "complete.obs"), 
           hc.order = TRUE, outline.col = "white", type = "upper", lab = TRUE)
# recent_greening and extent of recent change are strongly correlated (0.93) --> 
# we keep recent_greening; productivity, temperature, growing season and bioclimatic 
# zone are also correlated --> we keep temperature (less NAs)

# categorical moderators
dt_10_f %>%  group_by(soil_chemistry) %>% summarize(n = n()) # there are a lot of missing values
dt_10_f %>%  group_by(soil_texture) %>% summarize(n = n()) # there are a lot of missing values
dt_10_f %>%  group_by(soil_moisture) %>% summarize(n = n()) # there are a lot of missing values
dt_10_f %>%  group_by(soil_type) %>% summarize(n = n()) %>% print(n = Inf) # there are a lot of missing values
dt_10_f %>%  group_by(soil_type_D) %>% summarize(n = n()) # we can include this :)
  ggplot(dt_10_f) + geom_boxplot(aes(soil_type_D, permafrost_D), na.rm = T)
dt_10_f %>%  group_by(permafrost) %>% summarize(n = n()) # there are a lot of missing values
dt_10_f %>%  group_by(habitat_type) %>% summarize(n = n()) # we can include this :)
dt_10_f %>%  group_by(habitat_type_D) %>% summarize(n = n()) # there are a lot of missing values



## 1.1 create model guide ------------------------------------------------------
model_guide <- CJ(eco_response = unique(dt_10_f$MA.value),
                  # define the list of variables that we want to include as moderators
                  # for categorical variables -1 removes intercept from model, which makes interpretation easier
                  variable = c("change_f -1", "herb_fgr_change",
                               "ess.se", "year.c", "study_length", "overall_criterion -1", 
                               "error_type -1", "spatial_resolution -1", "exclusion -1",
                               "permafrost_D", "recent_greening", 
                               "recent_warming", "soil_type_D -1", "distance_from_coast", 
                               "precipitation", "temperature", "distance_to_treeline", 
                               "elevation_DEM", "habitat_type -1")) %>% 
                # create 'selection' formula (tells the loop which MA.value to take)
                # and create formula for model
                mutate(select = paste0("MA.value == ", "'" , eco_response, "'"),
                       formula = paste("yi_smd ~ ", variable),
                       # add a unique ID for each model (remove -1 from name)
                       model_id = gsub(" -1", "", paste(eco_response, variable, sep = "_")))
set.seed(161)


## 1.2. model for loop ---------------------------------------------------------
# first create empty data tables to store the outputs
tmp <- data.table()
model.results <- data.table(eco_response = NA, variable = NA,  model_id = NA, model_call = NA,
                          LRT = NA, pval = NA, BIC = NA, AIC = NA, AICc = NA,
                          estimate = NA, ci.lb = NA, ci.ub = NA, estimate_pval = NA,
                          r2_marginal = NA)

# set the working directory to data (where the folder "builds" is)
# setwd("./data")

# we have to include different variables as moderators: 
# - publication bias: small study effect and decline effect
# - critical appraisal (overall_criterion)
# - scale dependence (spatial_resolution)
# - ecological modifiers

# to assess small study effect we need a helper function 
# to calculate adapted sampling variance based on effective sampling size (tilde n)
ess.var_cal <- function(dat){1/dat$sample_size_higher_diversity + 1/dat$sample_size_lower_diversity} 

tic()
for(i in 1:nrow(model_guide)){
  result <- tryCatch({ 
    # take subset of data for each response variable 
    data.sub = dt_10_f[eval(parse(text = model_guide[i, ]$select)), ]
    
    # calculate additional variables
    data.sub$ess.se <- sqrt(ess.var_cal(data.sub)) # to assess small study effect 
    data.sub$year.c <- data.sub$year - mean(data.sub$year) # center year to estimate decline effect
     
    m0 <- rma.mv(yi_smd ~ 1, # intercept only model
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data.sub, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID),
                 data = data.sub, 
                 method = "ML",  test = "t", dfs = "contain") 
    
    m1 <- rma.mv(as.formula(model_guide[i, ]$formula), # model with the variables for which we want to test the effect
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data.sub, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID), 
                 data = data.sub, 
                 method = "ML",  test = "t", dfs = "contain") 
    
    an <- anova(m0, m1) # this tells you if adding the variable improves the model (if p < 0.05, it does)
    
    # store the outputs of each model (m1) in a temporary object
        # info about the model
        tmp$eco_response <- model_guide[i, ]$eco_response
        tmp$variable <- model_guide[i, ]$variable
        tmp$model_id <- model_guide[i, ]$model_id
        tmp$model_call <- paste0("builds/models/variable_models/", model_guide[i, ]$model_id,".rds")

        # comparisons of the models to the null model
        tmp$LRT <- an$LRT
        tmp$pval <- an$pval
        tmp$BIC <- an$fit.stats.f[4]
        tmp$AIC <- an$fit.stats.f[3]
        tmp$AICc <- an$fit.stats.f[5]
        
        # model results
        tmp$estimate <- m1$b[2] # estimate, ci.lb, ci.ub and estimate_pval are more or less useless for categorical variables but nice for continuous ones 
        tmp$ci.lb <- m1$ci.lb[2]
        tmp$ci.ub <- m1$ci.ub[2] 
        tmp$estimate_pval <- m1$pval[2] 
        tmp$r2_marginal <- r2_ml(m1)[1] # goodness-of-fit (R^2): marginal R2
        
    # combine all model results into a single object
    model.results <- rbind(model.results, tmp)
    
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(result)) {
    next
  }
  write_rds(m1, paste0("builds/models/variable_models/", model_guide[i, ]$model_id, ".rds"))
  cat(i,"/", nrow(model_guide), "\r")
}
toc()

model.results <- model.results[!is.na(model.results$eco_response)]
model.results$row_number <- c(1:nrow(model.results)) 

sig <- model.results[pval < 0.05,] # note that this is the anova/likelihood ratio test
sig # these are the models where adding the variable improved the model (28)

fwrite(model.results, "builds/model_results/variable_model_results.csv")



# 2. multi-moderator models --------------------------------------------------

# (if skipping step 1, then need to run this first:
# dt_10_f <- dt_f %>% filter(MA.value %in% subset(resp_var_f, nr_articles >= 10)$MA.value) # 2232 studies
# set the working directory to data (where the folder "builds" is) and
# load model results from step 1:
# model.results <- fread("builds/model_results/variable_model_results.csv")
# sig <- model.results[pval < 0.05,] # check which variables are significant (28)

# for those moderators that resulted significant in the univariate models
# we will build multi-moderator models. Same as above, we start by defining
# a model guide:
model_guide_multi <- CJ(eco_response = unique(sig$eco_response)) %>% 
                # create 'selection' formula (tells the loop which moderators to take)
                # and create formula for model (we do it manually...)
                mutate(select = paste0("MA.value == ", "'" , eco_response, "'"),
                       formula = c("yi_smd ~ change_f + exclusion + spatial_resolution -1",  # plant CN ratio
                                   "yi_smd ~ overall_criterion -1",                          # plant C content
                                   "yi_smd ~ exclusion -1",                                  # bryophyte abundance
                                   "yi_smd ~ herb_fgr_change + exclusion + 
                                      spatial_resolution -1",                                # dwarf shrub abundance
                                   "yi_smd ~ change_f -1",                                   # graminoid abundance
                                   "yi_smd ~ change_f -1",                                   # lichen abundance
                                   "yi_smd ~ error_type -1",                                 # tall shrub abundance
                                   "yi_smd ~ distance_to_treeline + habitat_type -1",        # plant diversity
                                   "yi_smd ~ error_type -1",                                 # plant fitness
                                   "yi_smd ~ change_f + error_type + exclusion + 
                                      habitat_type + overall_criterion + study_length -1",   # plant height
                                   "yi_smd ~ permafrost_D + recent_warming + temperature",   # plant productivity
                                   "yi_smd ~ permafrost_D",                                  # plant species richness
                                   "yi_smd ~ year.c",                                        # plant structure
                                   "yi_smd ~ error_type -1",                                 # soil C labile
                                   "yi_smd ~ error_type + year.c -1"))                       # soil moisture
set.seed(161)

# create empty data tables to store outputs
tmp <- data.table()
model.results.M <- data.table(eco_response = NA, formula = NA, model_call = NA,
                          LRT = NA, pval = NA, BIC = NA, AIC = NA, AICc = NA, deltaAICc = NA)
# and create an empty list to store models
model_list <- list()

# we run a loop to check if the multi-moderator models improve the models relative 
# to the intercept-only models

tic()
for(i in 1:nrow(model_guide_multi)){
  result <- tryCatch({ 
    # take subset of data for each response variable 
    data.sub = dt_10_f[eval(parse(text = model_guide_multi[i, ]$select)), ]
    
    # calculate additional variables
    data.sub$year.c <- data.sub$year - mean(data.sub$year) # center year to estimate decline effect
    
    m0 <- rma.mv(yi_smd ~ 1, # intercept only model
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data.sub, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID),
                 data = data.sub, 
                 method = "ML",  test = "t", dfs = "contain") 
    
    m1 <- rma.mv(as.formula(model_guide_multi[i, ]$formula), # model with the variables for which we want to test the effect
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data.sub, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID), 
                 data = data.sub, 
                 method = "ML",  test = "t", dfs = "contain") 
    
    an <- anova(m0, m1) # this tells you if adding the variables improves the model (if p < 0.05, it does)
    
    # store the outputs of each model (m1) in a temporary object
    # info about the model
        tmp$eco_response <- model_guide_multi[i, ]$eco_response
        tmp$formula <- model_guide_multi[i, ]$formula
        tmp$model_call <- paste0("builds/models/variable_models_M/", model_guide_multi[i, ]$eco_response, ".rds")
    # comparisons of the models to the null model
        tmp$LRT <- an$LRT
        tmp$pval <- an$pval
        tmp$BIC <- an$fit.stats.f[4]
        tmp$AIC <- an$fit.stats.f[3]
        tmp$AICc <- an$fit.stats.f[5]
        tmp$deltaAICc <- an$fit.stats.f[5] - an$fit.stats.r[5] # calculate change in AICc
       
    # combine all model results into a single object
    model.results.M <- rbind(model.results.M, tmp)
    
    model_list[[i]] <- m1 # store each multi-moderator model in our list 
    names(model_list)[i] <- model_guide_multi[i, ]$eco_response
    
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(result)) {
    next
  }
  write_rds(m1, paste0("builds/models/variable_models_M/", model_guide_multi[i, ]$eco_response, ".rds"))
  cat(i,"/", nrow(model_guide_multi), "\r")
}
toc()

model.results.M <- model.results.M[!is.na(model.results.M$eco_response)]
model.results.M$row_number <- c(1:nrow(model.results.M)) 

sig_multi <- model.results.M[pval < 0.05,] # note that this is the anova/likelihood ratio test
sig_multi # these are the models where the multi-moderator models were significantly better
          # than the intercept-only models -- all of them :)

fwrite(model.results.M, "builds/model_results/variable_model_M_results.csv")

model.table <- model.results.M %>% 
                mutate(eco_response = gsub("_", " ", eco_response),
                       pval = ifelse(pval < 0.001, "<0.001", round(pval, 4)),
                       LRT = round(LRT, 2), BIC = round(BIC, 2), 
                       AIC = round(AIC, 2), AICc = round(AICc, 2), 
                       deltaAICc = round(deltaAICc, 2)) %>% 
                select(eco_response, formula, LRT, pval,
                       BIC, AIC, AICc, deltaAICc) %>% 
                # rename columns
                rename("outcome variable" = "eco_response")
               
model.t <- flextable(model.table) %>% 
            autofit() 
model.t

save_as_docx(model.t, path = "./tables/model_table.docx")


# 3. prediction / estimate plots -----------------------------------------------
# for each model above we will plot the effects of the moderators that significantly
# affect the effect of herbivore diversity

## 3.1 plant CN ratio ----------------------------------------------------------
m <- model_list[["plant_CN_ratio"]]
  summary(m) # exclusion and change_f have a significant effect
data <- dt_10_f[MA.value == "plant_CN_ratio", ]

# spatial resolution was dropped from the models because it was redundant,
# so we need to re-run the model without it
m <- rma.mv(yi_smd ~ change_f + exclusion -1, # model with the variables for which we want to test the effect
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID), 
                 data = data, 
                 method = "ML",  test = "t", dfs = "contain") 

# plot for each moderator:
# when we have two moderators, to make predictions for one moderator 
# we have to keep the values of the other moderator constant. In the case
# of categorical moderators, to calculate predictions across levels 
# we need to take into account the proportion of observations in each level
data %>% group_by(change_f) %>% summarize(n = n()/length(data$habitat_type))
data %>% group_by(exclusion) %>% summarize(n = n()/length(data$habitat_type))

newgrid <- CJ(change_f = unique(data$change_f), 
              exclusion = unique(factor(data$exclusion))) 
  pred1 <- rma_predictions(m, newgrid) %>% 
            # create columns for weights
            mutate(weights.exclusion = case_when(exclusion == "partial" ~ 0.0571,
                                                 exclusion == "complete" ~ 0.943), 
                   weights.change_f = case_when(change_f == "F1" ~ 0.486,
                                                change_f == "F2_F3" ~ 0.171,
                                                change_f == "F3" ~ 0.229,
                                                change_f == "inv" ~ 0.114))

#### change_f ------------------------------------------------------------------
# calculate predictions at average levels of the categorical moderators
pred_change_f <- pred1 %>% group_by(change_f) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.exclusion),
                              ci.lb = weighted.mean(ci.lb, w = weights.exclusion),
                              ci.ub = weighted.mean(ci.ub, w = weights.exclusion),
                              pi.lb = weighted.mean(pi.lb, w = weights.exclusion),
                              pi.ub = weighted.mean(pi.ub, w = weights.exclusion)) %>% 
                    # add significance direction 
                    mutate(sig = case_when(ci.lb > 0 & pred > 0 ~ "Positive effect",
                                           ci.ub < 0 & pred < 0 ~ "Negative effect",
                                           ci.ub > 0.05 & ci.lb < 0 ~ "Non-significant"))

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(change_f, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(change_f) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
change_plantCN <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(change_f, 
                                         levels = c("inv", "F1", "F2", "F3", 
                                                    "F1_F2", "F1_F3", "F2_F3", 
                                                    "F1_F2_F3")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              # relabel categories on x-axis
              scale_x_discrete(labels = c("inv" = "Inv", "F1" = "F1", "F2" = "F2", 
                                          "F3" = "F3", "F1+F2", "F1_F3" = "F1+F3", 
                                          "F2_F3" = "F2+F3", "F1_F2_F3" = "F1+F2+F3")) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred_change_f, 
                             aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = change_f)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data =  pred_change_f, 
                              aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = change_f, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Plant CN ratio", x = "Herbivore functional groups", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = change_f, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(change_plantCN) # print plot to screen


#### exclusion -----------------------------------------------------------------
# calculate predictions at average levels of the categorical moderators
pred_exclusion <- pred1 %>% group_by(exclusion) %>% 
                    # calculate weighted values
                    summarise(pred = weighted.mean(pred, w = weights.change_f),
                              ci.lb = weighted.mean(ci.lb, w = weights.change_f),
                              ci.ub = weighted.mean(ci.ub, w = weights.change_f),
                              pi.lb = weighted.mean(pi.lb, w = weights.change_f),
                              pi.ub = weighted.mean(pi.ub, w = weights.change_f)) %>% 
                    # add significance direction 
                    mutate(sig = case_when(ci.lb > 0 & pred > 0 ~ "Positive effect",
                                           ci.ub < 0 & pred < 0 ~ "Negative effect",
                                           ci.ub > 0.05 & ci.lb < 0 ~ "Non-significant"))

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(exclusion, article_ID) %>% 
                        summarize(ns = n()) %>% group_by(exclusion) %>% 
                          summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))

exclusion_plantCN <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = exclusion,  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred_exclusion, 
                             aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                 x = exclusion)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred_exclusion, 
                              aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                  x = exclusion, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, values = c("Negative effect" = "#006699FF", 
                                                 "Non-significant" = "#B3B8B3FF",
                                                 "Positive effect" = "#B35900FF")) +
              labs(title = "Plant CN ratio", x = "Exclusion", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = exclusion, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(exclusion_plantCN) # print plot to screen



## 3.2 plant C content ---------------------------------------------------------
m <- model_list[["plant_C_content"]]
  summary(m) # overall_criterion has a significant effect
data <- dt_10_f[MA.value == "plant_C_content", ]

#### bias ----------------------------------------------------------------------
newgrid <- CJ(overall_criterion = unique(data$overall_criterion)) 
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0.05 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(overall_criterion, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(overall_criterion) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
bias_plantC <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(overall_criterion, 
                                         levels = c("Overall low risk of bias", 
                                                    "Overall medium risk of bias", 
                                                    "Overall high risk of bias")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              # relabel categories on x-axis
              scale_x_discrete(labels = c("Overall low risk of bias" = "low", 
                                          "Overall medium risk of bias" = "medium",
                                          "Overall high risk of bias" = "high")) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = overall_criterion)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = overall_criterion, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Plant C content", x = "Overall risk of bias", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = overall_criterion, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(bias_plantC) # print plot to screen

## model without influential study
data2 <- dt_10_f[MA.value == "plant_C_content", ] %>% filter(!study_ID == "200_c")
m2 <- rma.mv(yi_smd ~ overall_criterion -1, # model with the variables for which we want to test the effect
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data2, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID), 
                 data = data2, 
                 method = "ML",  test = "t", dfs = "contain") 
  summary(m2) # overall_criterion no longer has a significant effect

# compare to intercept-only model to assess significance of overall_criterion  
m0b <- rma.mv(yi_smd ~ 1, # intercept only model
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data2, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID),
                 data = data2, 
                 method = "ML",  test = "t", dfs = "contain") 
anova(m0b, m2) 
  
  


## 3.3 bryophyte abundance -----------------------------------------------------
m <- model_list[["plant_abundance_bryophytes"]]
  summary(m) # exclusion did NOT have a significant effect
data <- dt_10_f[MA.value == "plant_abundance_bryophytes", ]



## 3.4 dwarf shrub abundance ---------------------------------------------------
m <- model_list[["plant_abundance_dwarf_shrubs"]]
  summary(m) # none of the moderators had a significant effect
data <- dt_10_f[MA.value == "plant_abundance_dwarf_shrubs", ]



## 3.5 graminoid abundance -----------------------------------------------------
m <- model_list[["plant_abundance_graminoids"]]
  summary(m) # change_f has a significant effect
data <- dt_10_f[MA.value == "plant_abundance_graminoids", ]

# plot for each moderator:
#### change_f ------------------------------------------------------------------
newgrid <- CJ(change_f = unique(data$change_f))
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0.05 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(change_f, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(change_f) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
change_graminoids <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(change_f, 
                                         levels = c("inv", "F1", "F2", "F3", 
                                                    "F1_F2", "F1_F3", "F2_F3", 
                                                    "F1_F2_F3")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              # relabel categories on x-axis
              scale_x_discrete(labels = c("inv" = "Inv", "F1" = "F1", "F2" = "F2", 
                                          "F3" = "F3", "F1+F2", "F1_F3" = "F1+F3", 
                                          "F2_F3" = "F2+F3", "F1_F2_F3" = "F1+F2+F3")) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = change_f)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = change_f, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Graminoid abundance", x = "Herbivore functional groups", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = change_f, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(change_graminoids) # print plot to screen



## 3.6 lichen abundance --------------------------------------------------------
m <- model_list[["plant_abundance_lichens"]]
  summary(m) # change_f has a significant effect
data <- dt_10_f[MA.value == "plant_abundance_lichens", ]

# plot for each moderator:
#### change_f ------------------------------------------------------------------
newgrid <- CJ(change_f = unique(data$change_f))
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0.05 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(change_f, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(change_f) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
change_lichens <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(change_f, 
                                         levels = c("inv", "F1", "F2", "F3", 
                                                    "F1_F2", "F1_F3", "F2_F3", 
                                                    "F1_F2_F3")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              # relabel categories on x-axis
              scale_x_discrete(labels = c("inv" = "Inv", "F1" = "F1", "F2" = "F2", 
                                          "F3" = "F3", "F1+F2", "F1_F3" = "F1+F3", 
                                          "F2_F3" = "F2+F3", "F1_F2_F3" = "F1+F2+F3")) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = change_f)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = change_f, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Graminoid abundance", x = "Herbivore functional groups", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = change_f, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(change_lichens) # print plot to screen


## 3.7 tall shrub abundance ----------------------------------------------------
m <- model_list[["plant_abundance_tall_shrubs"]]
  summary(m) # error_type has a significant effect
data <- dt_10_f[MA.value == "plant_abundance_tall_shrubs", ]

#### error type ----------------------------------------------------------------
newgrid <- CJ(error_type = unique(data$error_type)) 
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(error_type, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(error_type) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
error_tallshrubs <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(error_type, levels = c("reported", "CI", 
                                                                "IQR", "estimated")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = error_type)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = error_type, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Dwarf shrub abundance", x = "Error type", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = error_type, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(error_tallshrubs) # print plot to screen



## 3.8 plant diversity ---------------------------------------------------------
m <- model_list[["plant_diversity"]]
  summary(m) # habitat type has a significant effect
data <- dt_10_f[MA.value == "plant_diversity", ]

#### habitat type --------------------------------------------------------------
newgrid <- CJ(habitat_type = unique(data$habitat_type), 
              distance_to_treeline = mean(data$distance_to_treeline)) # keep constant
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0.05 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(habitat_type, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(habitat_type) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
hab_diversity <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(habitat_type, levels = c("barren", "graminoid tundra", 
                                                                  "prostrate-shrub tundra (shrubs <40 cm)",
                                                                  "erect-shrub tundra (shrubs >40 cm)",
                                                                  "wetland", "other", "not reported")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              # relabel categories on x-axis
              scale_x_discrete(labels = c("barren" = "barren", "graminoid tundra" = "graminoid", 
                                          "prostrate-shrub tundra (shrubs <40 cm)" = "dwarf shrub",
                                          "erect-shrub tundra (shrubs >40 cm)" = "tall shrub",
                                          "wetland" = "wetland", "other" = "other", 
                                          "not reported" = "not reported")) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, aes(y = pred, ymin = pi.lb, ymax = pi.ub, x = habitat_type)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, aes(y = pred, ymin = ci.lb, ymax = ci.ub, x = habitat_type, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Plant diversity", x = "Habitat type", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = habitat_type, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(hab_diversity) # print plot to screen


## 3.9 plant fitness -----------------------------------------------------------
m <- model_list[["plant_fitness"]]
  summary(m) # error type has a significant effect
data <- dt_10_f[MA.value == "plant_fitness", ]

#### error type ----------------------------------------------------------------
newgrid <- CJ(error_type = unique(data$error_type)) 
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(error_type, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(error_type) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
error_fitness <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(error_type, levels = c("reported", "CI", 
                                                                "IQR", "estimated")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = error_type)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = error_type, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Plant fitness", x = "Error type", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = error_type, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(error_fitness) # print plot to screen



## 3.10 plant height ----------------------------------------------------------
m <- model_list[["plant_height"]]
  summary(m) # study length and error type have a significant effect
data <- dt_10_f[MA.value == "plant_height", ]

# need to run the model again (drop overall criterion and exclusion)
m <- rma.mv(yi_smd ~ change_f + error_type + 
                     habitat_type + study_length -1, # model with the variables for which we want to test the effect
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID), 
                 data = data, 
                 method = "ML",  test = "t", dfs = "contain") 
  summary(m) # change_f, study length and error type have a significant effect

#### study length --------------------------------------------------------------
# proportion of observations in each level of the other categorical moderators
data %>% group_by(change_f) %>% summarize(n = n()/length(data$habitat_type))
data %>% group_by(error_type) %>% summarize(n = n()/length(data$habitat_type))
data %>% group_by(habitat_type) %>% summarize(n = n()/length(data$habitat_type))

newgrid <- CJ(change_f = unique(data$change_f),
              error_type = unique(data$error_type),
              habitat_type = unique(data$habitat_type),
              study_length = seq(min(data$study_length, na.rm = TRUE), 
                                 max(data$study_length, na.rm = TRUE), by = .1)) 
pred2 <- rma_predictions(m, newgrid) %>% 
            # create columns for weights; because we have three other moderators
            # we have to do it in 3 steps; we start with change_f
            mutate(weights.change = case_when(change_f == "F1" ~ 0.09,
                                             change_f == "F1_F3" ~ 0.11, 
                                             change_f == "F2" ~ 0.08,
                                             change_f == "F2_F3" ~ 0.13,
                                             change_f == "F3" ~ 0.59)) %>% 
            # now make our predictions for each value of study length,
            # but keeping the levels of the other moderators separate
            group_by(study_length, error_type, habitat_type) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.change),
                              ci.lb = weighted.mean(ci.lb, w = weights.change),
                              ci.ub = weighted.mean(ci.ub, w = weights.change),
                              pi.lb = weighted.mean(pi.lb, w = weights.change),
                              pi.ub = weighted.mean(pi.ub, w = weights.change)) %>% 
            # add weights for second moderator
            mutate(weights.error = case_when(error_type == "estimated" ~ 0.06,
                                             error_type == "IQR" ~ 0.02,
                                             error_type == "reported" ~ 0.92)) %>% 
            group_by(study_length, habitat_type) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.error),
                              ci.lb = weighted.mean(ci.lb, w = weights.error),
                              ci.ub = weighted.mean(ci.ub, w = weights.error),
                              pi.lb = weighted.mean(pi.lb, w = weights.error),
                              pi.ub = weighted.mean(pi.ub, w = weights.error)) %>% 
            # add weights for third moderator
            mutate(weights.habitat = case_when(habitat_type == "erect-shrub tundra (shrubs >40 cm)" ~ 0.01,
                                             habitat_type == "graminoid tundra" ~ 0.31,
                                             habitat_type == "prostrate-shrub tundra (shrubs <40 cm)" ~ 0.56,
                                             habitat_type == "wetland" ~ 0.12)) %>% 
            group_by(study_length) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.habitat),
                              ci.lb = weighted.mean(ci.lb, w = weights.habitat),
                              ci.ub = weighted.mean(ci.ub, w = weights.habitat),
                              pi.lb = weighted.mean(pi.lb, w = weights.habitat),
                              pi.ub = weighted.mean(pi.ub, w = weights.habitat))

# plot the figure
length_height <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 6),
                          aes(x = study_length, y = yi_smd, size = 1/vi_smd, color = yi_smd), 
                          width = 0.1, alpha = 0.5) +
              # remove decimal places on x axis
              scale_x_continuous(labels = number_format(accuracy = 1)) +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # plot predicted values for dwarf shrub tundra
              geom_ribbon(data = pred2, 
                          aes(x = study_length, ymin = ci.lb, ymax = ci.ub), 
                          fill = "grey90", alpha = 0.7) +
              geom_line(data = pred2, 
                        aes(x = study_length, y = pred), color = "black") +
              labs(title = "Plant height", x = "Study length", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              scale_colour_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              theme_systrev() + 
                theme(legend.position = "none") 
  print(length_height)  # print plot to screen


#### change_f ------------------------------------------------------------------
# proportion of observations in each level of the other categorical moderators
data %>% group_by(error_type) %>% summarize(n = n()/length(data$habitat_type))
data %>% group_by(habitat_type) %>% summarize(n = n()/length(data$habitat_type))

newgrid <- CJ(change_f = unique(data$change_f),
              error_type = unique(data$error_type),
              habitat_type = unique(data$habitat_type),
              study_length = mean(data$study_length, na.rm = TRUE)) # keep constant
              
pred2 <- rma_predictions(m, newgrid) %>% 
            # create columns for weights; because we have two other moderators
            # we have to do it in 2 steps; we start with error_type
            mutate(weights.error = case_when(error_type == "estimated" ~ 0.06,
                                             error_type == "IQR" ~ 0.02,
                                             error_type == "reported" ~ 0.92)) %>% 
            group_by(change_f, habitat_type) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.error),
                              ci.lb = weighted.mean(ci.lb, w = weights.error),
                              ci.ub = weighted.mean(ci.ub, w = weights.error),
                              pi.lb = weighted.mean(pi.lb, w = weights.error),
                              pi.ub = weighted.mean(pi.ub, w = weights.error)) %>% 
            # add weights for third moderator
            mutate(weights.habitat = case_when(habitat_type == "erect-shrub tundra (shrubs >40 cm)" ~ 0.01,
                                             habitat_type == "graminoid tundra" ~ 0.31,
                                             habitat_type == "prostrate-shrub tundra (shrubs <40 cm)" ~ 0.56,
                                             habitat_type == "wetland" ~ 0.12)) %>% 
            group_by(change_f) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.habitat),
                              ci.lb = weighted.mean(ci.lb, w = weights.habitat),
                              ci.ub = weighted.mean(ci.ub, w = weights.habitat),
                              pi.lb = weighted.mean(pi.lb, w = weights.habitat),
                              pi.ub = weighted.mean(pi.ub, w = weights.habitat)) %>% 
            # add significance direction 
            mutate(sig = case_when(ci.lb > 0 & pred > 0 ~ "Positive effect",
                                   ci.ub < 0 & pred < 0 ~ "Negative effect",
                                   ci.ub > 0.05 & ci.lb < 0 ~ "Non-significant"))

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(change_f, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(change_f) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))

change_plantheight <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(change_f, 
                                         levels = c("inv", "F1", "F2", "F3", 
                                                    "F1_F2", "F1_F3", "F2_F3", 
                                                    "F1_F2_F3")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              # relabel categories on x-axis
              scale_x_discrete(labels = c("inv" = "Inv", "F1" = "F1", "F2" = "F2", 
                                          "F3" = "F3", "F1+F2", "F1_F3" = "F1+F3", 
                                          "F2_F3" = "F2+F3", "F1_F2_F3" = "F1+F2+F3")) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred2, 
                             aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = change_f)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data =  pred2, 
                              aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = change_f, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Plant height", x = "Herbivore functional groups", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = change_f, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(change_plantheight)


#### error type ----------------------------------------------------------------
# proportion of observations in each level of the other categorical moderators
data %>% group_by(change_f) %>% summarize(n = n()/length(data$habitat_type))
data %>% group_by(habitat_type) %>% summarize(n = n()/length(data$habitat_type))

newgrid <- CJ(change_f = unique(data$change_f),
              error_type = unique(data$error_type),
              habitat_type = unique(data$habitat_type),
              study_length = mean(data$study_length, na.rm = TRUE)) # keep constant
              
pred2 <- rma_predictions(m, newgrid) %>% 
            # create columns for weights; because we have two other moderators
            # we have to do it in 2 steps; we start with change_f
            mutate(weights.change = case_when(change_f == "F1" ~ 0.09,
                                             change_f == "F1_F3" ~ 0.11, 
                                             change_f == "F2" ~ 0.08,
                                             change_f == "F2_F3" ~ 0.13,
                                             change_f == "F3" ~ 0.59)) %>% 
            # now make our predictions for each value of error_type
            group_by(error_type, habitat_type) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.change),
                              ci.lb = weighted.mean(ci.lb, w = weights.change),
                              ci.ub = weighted.mean(ci.ub, w = weights.change),
                              pi.lb = weighted.mean(pi.lb, w = weights.change),
                              pi.ub = weighted.mean(pi.ub, w = weights.change)) %>% 
            # add weights for second moderator
            mutate(weights.habitat = case_when(habitat_type == "erect-shrub tundra (shrubs >40 cm)" ~ 0.01,
                                             habitat_type == "graminoid tundra" ~ 0.31,
                                             habitat_type == "prostrate-shrub tundra (shrubs <40 cm)" ~ 0.56,
                                             habitat_type == "wetland" ~ 0.12)) %>% 
            group_by(error_type) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.habitat),
                              ci.lb = weighted.mean(ci.lb, w = weights.habitat),
                              ci.ub = weighted.mean(ci.ub, w = weights.habitat),
                              pi.lb = weighted.mean(pi.lb, w = weights.habitat),
                              pi.ub = weighted.mean(pi.ub, w = weights.habitat)) %>% 
                    # add significance direction 
                    mutate(sig = case_when(ci.lb > 0 & pred > 0 ~ "Positive effect",
                                           ci.ub < 0 & pred < 0 ~ "Negative effect",
                                           ci.ub > 0.05 & ci.lb < 0 ~ "Non-significant"))

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(error_type, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(error_type) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
error_plantheight <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(error_type, levels = c("reported", "CI", 
                                                                "IQR", "estimated")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred2, aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = error_type)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred2, aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = error_type, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Plant height", x = "Error type", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = error_type, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(error_plantheight) # print plot to screen


## 3.11 plant productivity -----------------------------------------------------
m <- model_list[["plant_productivity"]]
  summary(m) # permafrost has a significant effect
data <- dt_10_f[MA.value == "plant_productivity", ]

#### permafrost ----------------------------------------------------------------
newgrid <- CJ(permafrost_D = seq(min(data$permafrost_D, na.rm = TRUE), 
                                 max(data$permafrost_D, na.rm = TRUE), by = .1),
              temperature = mean(data$temperature),       # keep temperature constant
              recent_warming = mean(data$recent_warming)) # keep recent warming constant
pred2 <- rma_predictions(m, newgrid)

permafrost_plantprod <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 6),
                          aes(x = permafrost_D, y = yi_smd, size = 1/vi_smd, 
                              color = yi_smd), width = 0.1, alpha = 0.5) +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # plot predicted values for dwarf shrub tundra
              geom_ribbon(data = pred2, 
                          aes(x = permafrost_D, ymin = ci.lb, ymax = ci.ub), 
                          fill = "grey90", alpha = 0.7) +
              geom_line(data = pred2, 
                        aes(x = permafrost_D, y = pred), color = "black") +
              labs(title = "Plant productivity", x = "Permafrost (%)", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              scale_colour_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              theme_systrev() + 
                theme(legend.position = "none") 
  print(permafrost_plantprod)  # print plot to screen

  

## 3.12 plant species richness -------------------------------------------------
m <- model_list[["plant_species_richness"]]
  summary(m) # permafrost does NOT have a significant effect
data <- dt_10_f[MA.value == "plant_species_richness", ]



## 3.13 plant structure ----------------------------------------------------------
m <- model_list[["plant_structure"]]
  summary(m) # year.c has a significant effect
data <- dt_10_f[MA.value == "plant_structure", ]

# calculate additional variables
data$year.c <- data$year - mean(data$year) # center year to estimate decline effect

#### year ----------------------------------------------------------------------
newgrid <- CJ(year.c = seq(min(data$year.c, na.rm = TRUE), 
                           max(data$year.c, na.rm = TRUE), by = .1)) 
pred2 <- rma_predictions(m, newgrid)
    
year_plantstruct <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 6),
                          aes(x = year.c, y = yi_smd, size = 1/vi_smd, color = yi_smd), 
                          width = 0.1, alpha = 0.5) +
              # remove decimal places on x axis
              scale_x_continuous(labels = number_format(accuracy = 1)) +
              geom_hline(yintercept = 0, linetype = "dashed") +
              geom_ribbon(data = pred2, 
                          aes(x = year.c, ymin = ci.lb, ymax = ci.ub), 
                          fill = "grey90", alpha = 0.7) +
              geom_line(data = pred2, 
                        aes(x = year.c, y = pred), color = "black") +
              labs(title = "Plant structure", x = "Publication year (centered)", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              scale_colour_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              theme_systrev() + 
                theme(legend.position = "none") 
  print(year_plantstruct)  # print plot to screen
  
  

## 3.14 soil C labile ---------------------------------------------------------
m <- model_list[["soil_C_labile"]]
  summary(m) # error type does NOT have a significant effect
data <- dt_10_f[MA.value == "soil_C_labile", ]



## 3.15 soil moisture -----------------------------------------------------
m <- model_list[["soil_moisture"]]
  summary(m) # year and error type have a significant effect
data <- dt_10_f[MA.value == "soil_moisture", ]

# calculate additional variables
    data$year.c <- data$year - mean(data$year) # center year to estimate decline effect

#### year ----------------------------------------------------------------------
newgrid <- CJ(year.c = seq(min(data$year.c, na.rm = TRUE), 
                           max(data$year.c, na.rm = TRUE), by = .1), 
              error_type = unique(data$error_type)) 
    
# proportion of observations in each level of the categorical moderator
data %>% group_by(error_type) %>% summarize(n = n()/length(data$habitat_type))

pred2 <- rma_predictions(m, newgrid) %>% 
            # create columns for weights
            mutate(weights.error = case_when(error_type == "CI" ~ 0.0312,
                                             error_type == "estimated" ~ 0.0312, 
                                             error_type == "reported" ~ 0.938)) %>% 
            # now make our predictions for each value of publication year
            group_by(year.c) %>% 
                    summarise(pred = weighted.mean(pred, w = weights.error),
                              ci.lb = weighted.mean(ci.lb, w = weights.error),
                              ci.ub = weighted.mean(ci.ub, w = weights.error),
                              pi.lb = weighted.mean(pi.lb, w = weights.error),
                              pi.ub = weighted.mean(pi.ub, w = weights.error))

# plot the figure
year_soilmoist <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 6),
                          aes(x = year.c, y = yi_smd, size = 1/vi_smd, color = yi_smd), 
                          width = 0.1, alpha = 0.5) +
              # remove decimal places on x axis
              scale_x_continuous(labels = number_format(accuracy = 1)) +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # plot predicted values for dwarf shrub tundra
              geom_ribbon(data = pred2, 
                          aes(x = year.c, ymin = ci.lb, ymax = ci.ub), 
                          fill = "grey90", alpha = 0.7) +
              geom_line(data = pred2, 
                        aes(x = year.c, y = pred), color = "black") +
              labs(title = "Soil moisture", x = "Publication year (centered)", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              scale_colour_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              theme_systrev() + 
                theme(legend.position = "none") 
  print(year_soilmoist)  # print plot to screen

#### error type ----------------------------------------------------------------
newgrid <- CJ(year.c = mean(data$year.c), # keep year constant
              error_type = unique(data$error_type)) 
  pred1 <- rma_predictions(m, newgrid)
    # add significance direction 
    pred1[ci.lb > 0 & pred > 0, sig := "Positive effect"]
    pred1[ci.ub < 0 & pred < 0, sig := "Negative effect"]
    pred1[ci.ub > 0 & ci.lb < 0, sig := "Non-significant"]

# create data table for sample sizes (for categorical moderators only)
sample_sizes <- data %>% group_by(error_type, article_ID) %>% 
                              summarize(ns = n()) %>% group_by(error_type) %>% 
                                summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
# plot the figure
error_soilmoist <- ggplot() +
              geom_jitter(data = subset(data, abs(yi_smd) < 20), # this is all data
                          aes(x = factor(error_type, levels = c("reported", "CI", 
                                                                "IQR", "estimated")),  
                              y = yi_smd, size = 1/vi_smd, fill = yi_smd), 
                          shape = 21, colour= "lightgrey", stroke = 0, width = 0.1, alpha = 0.9) +
              scale_fill_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              # add the model results
              new_scale_fill() +
              # add prediction intervals
              geom_linerange(data = pred1, 
                             aes(y = pred, ymin = pi.lb, ymax = pi.ub, 
                                               x = error_type)) +
              # add predicted effect sizes and confidence intervals
              geom_pointrange(data = pred1, 
                              aes(y = pred, ymin = ci.lb, ymax = ci.ub, 
                                                x = error_type, fill = sig),
                              shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
              scale_fill_manual(NULL, #name = "Effect size",
                                values = c("Negative effect" = "#006699FF", 
                                           "Non-significant" = "#B3B8B3FF",
                                           "Positive effect" = "#B35900FF")) +
              labs(title = "Soil moisture", x = "Error type", 
                   y = "Effect size (Hedges g)", fill = "Effect size") +
              ylim(-6, 6) +
              guides(size = FALSE) +
              geom_text(data = sample_sizes, # indicate number of articles(studies)
                        aes(x = error_type, y = -6, label = sample_size),
                        size = 3, color = "grey10") +
              theme_systrev() + 
                theme(legend.position = "none") 
print(error_soilmoist) # print plot to screen





## 4. figures for manuscript ---------------------------------------------------
#### Figure 6. herbivore diversity ---------------------------------------------
# subset of outcome variables for which change_f and exclusion 
# were significant in multi-moderator models
fig6 <- ggarrange(change_plantCN + ggtitle("a. Plant C:N"), 
          change_graminoids + ggtitle("b. Graminoid abundance"), 
          change_lichens + ggtitle("c. Lichen abundance"),
          change_plantheight + ggtitle("d. Plant height"), 
          exclusion_plantCN + ggtitle("e. Plant C:N"), 
          ncol = 3, nrow = 2, 
          legend = "bottom", common.legend = T)
fig6

# save the figure to a separate file
ggsave(fig6, file = "figures/Fig6.png", dpi = 600)


#### Figure 7. methodological moderators ---------------------------------------
# subset of outcome variables for which risk of bias, publication year,
# spatial resolution, study length and error type were significant in multi-moderator models
fig7 <- ggarrange(bias_plantC + ggtitle("a. Plant C content"), 
          year_plantstruct + ggtitle("b. Plant structure"), 
          year_soilmoist + ggtitle("c. Soil moisture"), 
          length_height + ggtitle("d. Plant height"), 
          error_tallshrubs + ggtitle("e. Dwarf shrub abundance"),
          error_fitness + ggtitle("f. Plant fitness"), 
          error_plantheight + ggtitle("g. Soil labile C"), 
          error_soilmoist + ggtitle("h. Soil moisture"), 
          ncol = 4, nrow = 2, 
          legend = "bottom", common.legend = T)
fig7

# save the figure to a separate file
ggsave(fig7, file = "figures/Fig7.png", dpi = 600)


#### Figure 8. ecological moderators ---------------------------------------
# subset of outcome variables for which habitat_type and permafrost were significant in multi-moderator models
fig8 <- ggarrange(hab_diversity + ggtitle("a. Plant diversity"), 
          permafrost_plantprod + ggtitle("b. Plant productivity"), 
          ncol = 2, nrow = 1, 
          legend = "bottom", common.legend = T)
fig8

# save the figure to a separate file
ggsave(fig8, file = "figures/Fig8.png", dpi = 600)


## Figures for presentations
fig6b <- ggarrange(change_graminoids + ggtitle("a. Graminoid abundance"), 
                   change_lichens + ggtitle("b. Lichen abundance"),
                   change_plantheight + ggtitle("c. Plant height"), 
                   ncol = 3, nrow = 1, 
                   legend = "bottom", common.legend = T)
fig6b


## 5. sensitivity analysis (leave-one-out) -------------------------------------
cooks.distance_results <- cooks.distance(m0, progbar=TRUE,  reestimate=TRUE, 
                                         parallel="snow", ncpus=4, cluster = study_ID)

# add these values into a dataframe
cooks_distance_df <- as.data.frame(cooks.distance_results)
cooks_distance_df$study_ID <- rownames(cooks_distance_df)
cooks_distance_df <- cooks_distance_df[, c("study_ID", "cooks.distance_results")]
colnames(cooks_distance_df)[2] <- "cooks_distance"

cooks_distance_df %>% filter(cooks_distance > 1) #study 395_e

plot(cooks_plant_C_content, type = "o", pch = 19, xlab = "article", ylab = "Cook's Distance", xaxt = "n")
axis(side = 1, at = seq_along(cooks.distance_results), labels = as.numeric(names(cooks.distance_results)))


# assuming MA.value is a column in your dataframe dt_10_f
unique_ma_values <- unique(dt_10_f$MA.value)

cooks_list <- list()  # create an empty list to store results

for (i in seq_along(unique_ma_values)) {
  
  ma_value <- unique_ma_values[i]
  
  # filter the dataframe for the current MA.value
  data.sub <- dt_10_f[dt_10_f$MA.value == ma_value, ]
  
  m0 <- rma.mv(
    yi_smd ~ 1,
    V = vcalc(vi_smd,
              cluster = article_ID,
              obs = study_ID,
              data = data.sub,
              rho = 0.5),
        random = list(~ 1 | article_ID / study_ID),
        data = data.sub,
        method = "ML",
        test = "t",
        dfs = "contain")
  
  # calculate Cook's distance
  cooks_distance_results <- cooks.distance(
    m0,
    progbar = TRUE,
    reestimate = TRUE,
    parallel = "snow",
    ncpus = 4,
    cluster = study_ID)
  
  # Store the results in the list
  cooks_list[[i]] <- list(
    model_id = ma_value,
    cooks_distance_results = cooks_distance_results)
}

# combine results into a dataframe
cooks_results_df <- do.call(rbind, lapply(cooks_list, data.frame))
cooks_results_df$study_ID <- rownames(cooks_results_df)

# check values higher than 1 (influential studies)
cooks_results_df %>% filter(cooks_distance_results>1) #plant C content study 200_c
cooks_plant_C_content <- cooks_results_df %>% filter(model_id == "plant_C_content" )

plot(x = rep(seq_along(cooks_plant_C_content$study_ID), each = ncol(cooks_plant_C_content$cooks_distance_results)),
     y = as.vector(cooks_plant_C_content$cooks_distance_results),
     type = "o", pch = 19,
     xlab = "article", ylab = "Cook's Distance", xaxt = "n")

# add x-axis labels
axis(side = 1, at = seq_along(cooks_plant_C_content$study_ID), labels = cooks_plant_C_content$study_ID)


# plant C content plot
plot(y = rep(seq_along(cooks_plant_C_content$study_ID), each = ncol(cooks_plant_C_content$cooks_distance_results)),
    x = as.vector(cooks_plant_C_content$cooks_distance_results),
    type = "o", pch = 19,
    col = ifelse(cooks_plant_C_content$cooks_distance_results > 1, "#B35900FF", "#B3B8B3FF"),  
    xlab = "Cook's Distance", ylab = "study", yaxt = "n")

# labels
axis(side = 2, at = seq_along(cooks_plant_C_content$study_ID), labels = cooks_plant_C_content$study_ID, las = 1, cex.axis = 0.7)  

influential <- dt %>% filter(study_ID=="200_c")

names(coded_data)

