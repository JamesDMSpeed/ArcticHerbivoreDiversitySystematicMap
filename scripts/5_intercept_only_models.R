#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Systematic review on herbivory diversity in Arctic tundra
#                       Intercept only models
#     Jonas Trepel (jonas.trepel@gmail.com / jonas.trepel@bio.au.dk)
#                        19/20-April-2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# what happens in this script: 
# 1. summary of studies included in quantitative synthesis
# 2. intercept only models 
# 2.1 create model guide
# 2.2. loop models
# 3. plot
# 4. herbivore diversity or herbivore exclusion?

# in this script we take into account the grouping of herbivores
# into FUNCTIONAL GROUPS (sensu Speed et al 2019)

# libraries----
library(data.table)
library(tidyverse)
library(ggplot2)
library(harrypotter) # colour palettes for consistency with systematic map
library(ggnewscale)  # to be able to add several scales to ggplot
library(gridExtra)
library(metafor)
library(tictoc)
library(orchaRd)  # meta-analysis visualization package
                  # needed to calculate measures of heterogeneity for multilevel meta-analysis, marginal R2 etc
  # devtools::install_github("daniel1noble/orchaRd", force = TRUE)


# customised functions ---------------------------------------------------------
theme_systrev <- function(){ #create a new theme function for the style of graphs
  theme_bw() +                  #use a predefined theme as a base
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 10, color = "black"), 
          axis.title = element_text(size = 14, color = "black"),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.title = element_text(size = 18, vjust = 1, hjust = 0, color="black"),
          legend.text = element_text(size = 12, color = "black"),          
          legend.title = element_text(size = 10, color = "black"),
          legend.position = "bottom",
          legend.key.width = unit(2.5, "cm"),
          legend.key.height = unit(0.8, "lines"),
          legend.key.size = unit(0.8, "lines"),
          legend.box = "vertical", 
          legend.margin = margin())
}

# colour palette (we will use these colours later)
hp(n = 15, house = "Ravenclaw")[1] # colour for negative effects
hp(n = 15, house = "Ravenclaw")[9] # #B3B8B3FF colour for non-significant
hp(n = 15, house = "Ravenclaw")[15] # colour for positive effects

 
# load data --------------------------------------------------------------------
# clean the environment first
rm(list=ls())

# set the working directory to the folders where the data files are stored
setwd("./data/effect_sizes")

# the data loaded here was built with script 4_calculate_effect_sizes
# and includes the coding database with the corresponding effect sizes
# file: data_with_effect_sizes.csv

# load the dataset
dt <- fread("data_with_effect_sizes.csv")


# 1. summary of studies included in quantitative synthesis --------------------
# add information for studies included in quantitative synthesis (lower part of ROSES): 
# - studies for which we could not calculate effect size (var=0, n=1, missing SD or non-convertible effect size)
# - excluded because no contrast according to our herbivore grouping
# - excluded because outliers generated when imputing missing values
# - included because outcome variable was reported by at least 5 or 10 articles

# nr of articles with no effect size
dt %>% filter(is.na(yi_smd)) %>% summarize(n = n()) # 287 studies

# nr of articles with no herbivore diversity contrast
dt %>% filter(change_f == "zero") %>% summarize(n = n()) # 133 studies

# clean up a bit
dt_f <- dt %>% filter(!is.na(yi_smd)) %>%   # remove missing values (286 studies) 
            # remove also studies for which there is no contrast (133 studies)
            filter(!change.long_f == "no_contrast") %>%
            # remove extreme values (36 studies)
            filter(yi_smd > quantile(dt$yi_smd, .005, na.rm = TRUE) & 
                   yi_smd < quantile(dt$yi_smd, .995, na.rm = TRUE)) %>% 
            # make sure article_ID is character
            mutate(article_ID = as.character(article_ID)) # dt_f has 3263 observations

length(dt$study_ID) - length(dt_f$study_ID) # excluding 450 studies so far

# count the number of articles/studies for each response variable
resp_var_f <- dt_f %>% group_by(var_class, MA.value, article_ID) %>% 
                summarize(ns = n()) %>%
                  group_by(var_class, MA.value) %>% 
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
              arrange(desc(nr_articles))

# list of studies measuring responses with more than 5 or 10 articles  
var_5_f <- dt_f %>% filter(MA.value %in% subset(resp_var_f, nr_articles >= 5)$MA.value) %>% 
  select(study_ID) %>% mutate(MA.5 = "included (5 articles)") # 2846 studies
var_10_f <- dt_f %>% filter(MA.value %in% subset(resp_var_f, nr_articles >= 10)$MA.value) %>% 
  select(study_ID) %>% mutate(MA.10 = "included (10 articles)") # 2232 studies
length(var_5_f$study_ID) - length(var_10_f$study_ID) 
    # 614 excluded studies from meta-regressions because the outcome was reported
    # by less than 10 articles

# nr of excluded studies because they were reported by less than five articles
dt_f %>% filter(MA.value %in% subset(resp_var_f, nr_articles < 5)$MA.value) %>% 
  summarize(n = n()) # 417 studies

# nr of studies excluded from quantitative synthesis
length(dt$study_ID) - length(var_5_f$study_ID)


## DO THIS FOR THE FINAL ANALYSES to add to Additional file 3
# (add row to the coded_data file indicating whether a study was included or not in quantitative synthesis and why)
dt_f_all <- dt %>% left_join(var_5_f, by = "study_ID") %>% left_join(var_10_f, by = "study_ID") %>% 
             mutate(MA = ifelse(is.na(MA.10), MA.5, MA.10),
                    quantitative_synthesis = case_when(is.na(yi_smd) ~ "not possible to estimate effect size",
                                                       change.long_f == "no_contrast" ~ "no diversity contrast",
                                                       !is.na(MA) ~ MA,
                                                       TRUE ~ "excluded (<5 articles)"))
dt_f_all %>% distinct(quantitative_synthesis)

quant_synthesis <- dt_f_all %>% select(study_ID, quantitative_synthesis)
fwrite(quant_synthesis, "quant_synthesis.csv")


# overview of outcome variables
outcome.vars <- dt_f %>% 
                  select("article_ID", "study_ID", 
                         "var_class", "MA.value", "measured_response_variable_new") %>% 
                  group_by(var_class, MA.value, article_ID) %>% 
                    summarize(ns = n()) %>%
                  group_by(var_class, MA.value) %>% 
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
              arrange(desc(nr_articles)) %>% 
              mutate(MA.value = case_when(MA.value == "plant_productivity" ~ "primary_productivity",
                                          MA.value == "gross_primary_productivity" ~ "gross_primary_production",
                                          MA.value == "plant_CN_ratio" ~ "plant_C:N",
                                          MA.value == "soil_CN_ratio" ~ "soil_C:N",
                                          TRUE ~ MA.value)) %>% 
              # remove underscores for plotting
              mutate(MA.value = str_replace_all(MA.value, "_", " ")) 
length(subset(outcome.vars, nr_articles >= 5)$MA.value) # 47 variables have 5 articles or more
length(subset(outcome.vars, nr_articles >=10)$MA.value) # 25 variables have 10 articles or more

length(subset(outcome.vars, nr_articles == 1)$MA.value) # 23 variables were reported in only one article


### Figure 4 -------------------------------------------------------------------
# plot the number of articles for each outcome variable reported by > 5 articles
# define our colour palette 
pal <- hp(n = 8, house = "Ravenclaw")
# link factor names to the colours 
names(pal) <- c("ecosystem", "disease", "fungal", "invertebrate", "microbial", "soil", "herbivory", "plant")

fig4 <- ggplot(outcome.vars[outcome.vars$nr_articles > 4,]) +
          geom_bar(aes(x = reorder(MA.value, nr_articles), y = nr_articles, 
                       fill = var_class), stat = "identity") +
          labs(x = " ", y = "Number of articles") +
          geom_hline(yintercept = 10, linetype = "dashed", size = 0.6) +
          annotate(geom = "text", y = 10.5, x = 5,
                     label = "10 articles", angle = 90, vjust = 1) +
          geom_text(aes(x = reorder(MA.value, nr_articles), y = nr_articles + 0.6, 
                        label = nr_studies), size = 2.5) +
          scale_fill_manual(values = pal,  # using our palette
                            name = "Outcome variable") +
          theme_systrev() +
            theme(axis.text = element_text(size = 8)) +
           coord_flip() 

# save the figure to a separate file
# set the working directory to the main data folder
setwd("./data")
ggsave(fig4, file = "figures/Fig4.png", dpi = 600)


# pie chart for presentations:
# we need to have a summary of the data, build a stacked barchart with only one bar
# and then make it circular
outcome.class.vars <- dt_f %>% 
                  select("article_ID", "study_ID", 
                         "var_class") %>% 
                  group_by(var_class, article_ID) %>% 
                    summarize(ns = n()) %>%
                  group_by(var_class) %>% 
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
              arrange(desc(nr_articles))

# define our colour palette 
pal <- hp(n = 8, house = "Ravenclaw")
# link factor names to the colours 
names(pal) <- c("ecosystem", "disease", "fungal", "invertebrate", "microbial", "soil", "herbivory", "plant")

ggplot(outcome.class.vars, aes(x="", y = nr_studies, fill = var_class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = pal,  # using our palette
                    name = "Outcome variable") +
  theme_void()

# percentages of each var_class
data.frame(var_class = outcome.class.vars$var_class,
           percentage = round(outcome.class.vars$nr_studies/sum(outcome.class.vars$nr_studies)*100,2))



# 2. intercept only models -----------------------------------------------------
# here we will use studies reporting variables measured by at least 5 articles
length(subset(resp_var_f, nr_articles >= 5)$MA.value) # 47 variables have 5 articles or more

# filter the data to include only responses with more than 5 articles
dt_5_f <- dt_f %>% filter(MA.value %in% subset(resp_var_f, nr_articles >= 5)$MA.value) %>% 
       # create a grouping variable for plotting results of meta-analysis
       mutate(var_class = case_when(MA.value %in% c("plant_abundance_bryophytes", "plant_abundance_graminoids",
                                                    "plant_abundance_total", "plant_abundance_dwarf_shrubs",
                                                    "plant_abundance_forbs", "plant_abundance_lichens", "plant_abundance_woody_species",
                                                    "plant_abundance_tall_shrubs", "plant_abundance_belowground") ~ "plant abundance",
                                    MA.value %in% c("plant_quality", "plant_N_content", "plant_C_content", 
                                                    "plant_CN_ratio", "plant_defense", "plant_P_content") ~ "plant chemistry",
                                    MA.value %in% c("moss_depth", "plant_structure_change", 
                                                    "plant_population_dynamics", "plant_fitness", 
                                                    "plant_structure", "plant_leaf_size", "vegetation_structure",
                                                    "plant_height", 
                                                    "plant_species_richness", "litter_abundance", 
                                                    "plant_diversity") ~ "plant community",
                                    MA.value %in% c("microorganism_abundance", "fungal_abundance", "microbial_C_cycle",
                                                    "microbial_abundance", "microbial_N_content", "microbial_C_content") ~ "soil microorganisms",
                                    MA.value %in% c("soil_temperature", "soil_N_labile", "soil_moisture", "redox_conditions",
                                                    "soil_compaction", "soil_CN_ratio", "soil_C_total", "soil_C_labile",
                                                    "soil_N_total") ~ "soil properties",
                                    MA.value %in% c("gross_primary_productivity", "ecosystem_respiration", "plant_productivity", "C_turnover",
                                                    "herbivory_marks", "net_ecosystem_exchange") ~ "ecosystem processes", 
                                    TRUE ~ "missing!")) 
# just to check if we are missing any variable from our grouping
# (there should be 6 var_class, and no "missing!")
dt_5_f %>% distinct(var_class) %>% print(n = Inf)

# number of articles(studies) using size selective exclosures included in quantitative synthesis
dt_5_f %>% filter(size_selective_exclosures == "yes") %>% 
              group_by(article_ID) %>% 
                    summarize(ns = n()) %>%
              ungroup() %>% summarize(nr_articles = n(), nr_studies = sum(ns)) 
# 13 articles, 496 studies


## 2.1 create model guide ------------------------------------------------------
model_guide <- dt_5_f %>% distinct(MA.value, var_class) %>% 
                rename(eco_response = MA.value) %>% 
                # create 'selection' formula (tells the loop which MA.value to take)
                mutate(select = paste0("MA.value == ", "'" , eco_response, "'"),
                       # add a unique ID for each model
                       model_id = paste(eco_response, "intercept_only",sep = "_"))
set.seed(161)


## 2.2. model for loop ---------------------------------------------------------
# first create empty data tables to store model outputs
intercept.tmp <- data.table(var_class = NA, eco_response = NA, 
                            model_id = NA, model_call = NA,
                            estimate = NA, ci.lb = NA, ci.ub = NA, 
                            pi.lb = NA, pi.ub = NA, pval = NA)

eco.intercepts <- data.table(var_class = NA, eco_response = NA, 
                             model_id = NA, model_call = NA,
                             estimate = NA, ci.lb = NA, ci.ub = NA, 
                             pi.lb = NA, pi.ub = NA, pval = NA)

# and create an empty list to store orchard plots
orchard.plot_list <- list()

# make sure that there is a folder called "builds" in the working directory
# with the folders "models" and "null_models" in it
setwd("./data")

# run the models
tic()
for(i in 1:nrow(model_guide)){
  result <- tryCatch({ 
    # select subset of data for each response variable
    data.sub = dt_5_f[eval(parse(text = model_guide[i, ]$select)),]
    # build multi-level meta-analysis (only intercept)
    m0 <- rma.mv(yi_smd ~ 1, # intercept only model
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data.sub, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID),
                 data = data.sub, 
                 method = "REML",  test = "t", dfs = "contain") # method has to be changed to ML when comparing models with likelihood tests 
    # extract model results and add to intercept.tmp
        intercept.tmp$var_class <- model_guide[i, ]$var_class
        intercept.tmp$eco_response <- model_guide[i, ]$eco_response
        intercept.tmp$model_id <- model_guide[i, ]$model_id
        intercept.tmp$estimate <- m0$b[1]
        intercept.tmp$ci.lb <- m0$ci.lb[1]
        intercept.tmp$ci.ub <- m0$ci.ub[1] 
        intercept.tmp$pi.lb <- predict(m0)$pi.lb[1]
        intercept.tmp$pi.ub <- predict(m0)$pi.ub[1]
        intercept.tmp$pval <- m0$pval[1] 
    intercept.tmp$model_call <- paste0("builds/models/null_models/", model_guide[i, ]$model_id,".rds")
    eco.intercepts <- rbind(eco.intercepts, intercept.tmp)
    
    # orchard plot
    orchard.plot <- orchard_plot(m0, mod = "1", 
             xlab = "Effect size (Hedges g)", legend.pos = "none",
             group = "article_ID", k = TRUE, g = TRUE) + 
             scale_x_discrete(labels = c("Overall effect")) +
             labs(title = gsub("_", " ", model_guide[i, ]$eco_response))
    orchard.plot_list[[i]] <- orchard.plot # store each plot in our list 
    names(orchard.plot_list)[i] <- model_guide[i, ]$eco_response
    
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(result)) {
    next
  }
  write_rds(m0, paste0("builds/models/null_models/", model_guide[i, ]$model_id, ".rds"))
  cat(i,"/", nrow(model_guide), "\r")
}
toc()

model_master <- eco.intercepts[!is.na(eco.intercepts$eco_response)] # remove first row (all NAs)

sig <- model_master[pval < 0.05,]
sig # 9 intercept-only models are significant

fwrite(model_master, "builds/model_results/null_model_results.csv")

# check the orchard plots
do.call(grid.arrange, c(orchard.plot_list[1:12], ncol=4, nrow = 3))
do.call(grid.arrange, c(orchard.plot_list[13:24], ncol=4, nrow = 3))
do.call(grid.arrange, c(orchard.plot_list[25:36], ncol=4, nrow = 3))
do.call(grid.arrange, c(orchard.plot_list[37:47], ncol=4, nrow = 3))


# 3. plotting -------------------------
model_master <- model_master %>% 
                  ## add significance stars and direction
                  mutate(stars = case_when(pval < 0.05 & pval > 0.01 ~ "*",
                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                          pval <= 0.001 ~ "***",
                                          TRUE ~ ""),
                         sig = case_when(pval < 0.05 & estimate < 0 ~ "Negative effect",
                                        pval < 0.05 & estimate > 0 ~ "Positive effect",
                                        pval > 0.05 ~ "Non-significant")) 

# rename 
dt_5_fb <- dt_5_f %>% rename(eco_response = MA.value) %>% 
            # manually change names of eco_response for plotting
            mutate(eco_response = case_when(eco_response == "plant_abundance_belowground" ~ "belowground",
                                            eco_response == "plant_abundance_bryophytes" ~ "bryophytes",
                                            eco_response == "plant_abundance_dwarf_shrubs" ~ "dwarf_shrubs",
                                            eco_response == "plant_abundance_forbs" ~ "forbs",
                                            eco_response == "plant_abundance_graminoids" ~ "graminoids",
                                            eco_response == "plant_abundance_lichens" ~ "lichens",
                                            eco_response == "plant_abundance_tall_shrubs" ~ "tall_shrubs",
                                            eco_response == "plant_abundance_total" ~ "total",
                                            eco_response == "plant_abundance_woody_species" ~ "woody_species",
                                            eco_response == "plant_CN_ratio" ~ "plant_C:N",
                                            eco_response == "soil_CN_ratio" ~ "soil_C:N",
                                            TRUE ~ eco_response),
                   eco_response_ord = factor(gsub("_", " ", eco_response), # gsub messes up the order, so need to do it here already
                                             # reorder manually so that rows are in a logical way (from bottom to top)
                                             levels = c(
                                               # ecosystem processes
                                               "herbivory marks", "C turnover", "ecosystem respiration", "net ecosystem exchange", 
                                               "gross primary productivity", "plant productivity",
                                               # plant abundance
                                               "belowground", "bryophytes", "lichens", "woody species", "tall shrubs", "dwarf shrubs",  
                                               "forbs", "graminoids", "total",  
                                               # plant chemistry
                                               "plant quality", "plant defense", "plant C content", "plant C:N", 
                                               "plant N content", "plant P content", 
                                               # plant community
                                               "litter abundance", "moss depth", "plant species richness", "plant diversity", 
                                               "plant population dynamics", "plant fitness","plant height", "plant leaf size", 
                                               "plant structure change", "plant structure", "vegetation structure",
                                               #soil microorganisms
                                               "microbial N content", "microbial C content", "microbial C cycle", 
                                               "microbial abundance", "fungal abundance", "microorganism abundance", 
                                               # soil properties
                                               "redox conditions", "soil C:N", "soil N labile", "soil N total", 
                                               "soil C labile", "soil C total", 
                                               "soil compaction", "soil moisture", "soil temperature")))     
                   

model_masterb <- model_master %>% 
            # manually change names of model_master to match
            mutate(eco_response = case_when(eco_response == "plant_abundance_belowground" ~ "belowground",
                                            eco_response == "plant_abundance_bryophytes" ~ "bryophytes",
                                            eco_response == "plant_abundance_dwarf_shrubs" ~ "dwarf_shrubs",
                                            eco_response == "plant_abundance_forbs" ~ "forbs",
                                            eco_response == "plant_abundance_graminoids" ~ "graminoids",
                                            eco_response == "plant_abundance_lichens" ~ "lichens",
                                            eco_response == "plant_abundance_tall_shrubs" ~ "tall_shrubs",
                                            eco_response == "plant_abundance_total" ~ "total",
                                            eco_response == "plant_abundance_woody_species" ~ "woody_species",
                                            eco_response == "plant_CN_ratio" ~ "plant_C:N",
                                            eco_response == "soil_CN_ratio" ~ "soil_C:N",
                                            TRUE ~ eco_response),
                   eco_response_ord = factor(gsub("_", " ", eco_response), # gsub messes up the order, so need to do it here already
                                             # reorder manually so that rows are in a logical way (from bottom to top)
                                             levels = c(
                                               # ecosystem processes
                                               "herbivory marks", "C turnover", "ecosystem respiration", "net ecosystem exchange", 
                                               "gross primary productivity", "plant productivity",
                                               # plant abundance
                                               "belowground", "bryophytes", "lichens", "woody species", "tall shrubs", "dwarf shrubs",  
                                               "forbs", "graminoids", "total",  
                                               # plant chemistry
                                               "plant quality", "plant defense", "plant C content", "plant C:N", 
                                               "plant N content", "plant P content", 
                                               # plant community
                                               "litter abundance", "moss depth", "plant species richness", "plant diversity", 
                                               "plant population dynamics", "plant fitness","plant height", "plant leaf size", 
                                               "plant structure change", "plant structure", "vegetation structure",
                                               #soil microorganisms
                                               "microbial N content", "microbial C content", "microbial C cycle", 
                                               "microbial abundance", "fungal abundance", "microorganism abundance", 
                                               # soil properties
                                               "redox conditions", "soil C:N", "soil N labile", "soil N total", 
                                               "soil C labile", "soil C total", 
                                               "soil compaction", "soil moisture", "soil temperature")))
levels(as.factor(dt_5_fb$eco_response_ord))


# create data table for sample sizes 
sample_sizes <- dt_5_fb %>% group_by(var_class, eco_response_ord, article_ID) %>% 
                           summarize(ns = n()) %>% group_by(var_class, eco_response_ord) %>% 
                             summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
  
## Figure 5 ----
fig5 <- ggplot()+
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed")+
  geom_jitter(data = dt_5_fb, # could subset data to dt[abs(yi_smd) < 6,]
              aes(x = yi_smd, y = eco_response_ord, 
                  size = 1/vi_smd, color = yi_smd), #size = 1/vi_smd indicates model weight (the bigger the more influential)
              height = .1, inherit.aes = F, alpha = 0.5) +
  scale_colour_gradient2(name = "effect size", 
                         low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
  # add the model predictions
  new_scale_fill() +
  # add prediction intervals
  geom_linerange(data = model_masterb[], 
                 aes(x = estimate, xmin = pi.lb, xmax = pi.ub, 
                     y = eco_response_ord)) +
  # add predicted effect sizes and confidence intervals
  geom_pointrange(data = model_masterb[], 
                  aes(x = estimate, xmin = ci.lb, xmax = ci.ub,
                      fill = sig, y = eco_response_ord),
                  shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
  scale_fill_manual(NULL, values = c("Negative effect" = "#006699FF", 
                                     "Non-significant" = "#B3B8B3FF",
                                     "Positive effect" = "#B35900FF")) +
  geom_text(data = sample_sizes, # indicate number of articles and studies
            aes(x = -18, y = eco_response_ord, label = sample_size),
            size = 2.5, color = "grey10") +
  geom_text(data = model_masterb[],
            aes(y = eco_response_ord, x= 6, label = stars), 
            color = "grey10", size = 5) +
  facet_wrap(vars(var_class), scales = "free_y", ncol = 2) +
  xlab("Effect size") + ylab(NULL) +
  xlim(-20, 15) + 
  guides(size = FALSE) +
  theme_systrev() +
  theme(legend.position = "bottom",
    legend.key.width = unit(2.5, "cm"),
    legend.box = "vertical", legend.margin = margin(),
    legend.key.height = unit(0.8, "lines"),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 9))
fig5 

# save the figure to a separate file
ggsave(fig5, file = "figures/Fig5.png", dpi = 600)


# 4. herbivore diversity or herbivore exclusion? ------------------------------
# this is a nice plot for the intercept-only models
# however, since in many studies the contrast in herbivore diversity is given
# by exclusion of herbivores, these effects may simply reflect effect of herbivores
dt_5_f %>% filter(herbivore_fgr_ID_lower_diversity == "zero") %>% summarize(n = n())
(2337/2846)*100 # 82.1% of studies have zero as lower diversity!

# we can run the analyses again, but for studies not comparing to zero herbivores
dt_f0 <- dt_f %>% # remove studies where low diversity is zero
            filter(!herbivore_fgr_ID_lower_diversity == "zero") 

# count the number of articles/studies for each response variable
resp_var_f0 <- dt_f0 %>% group_by(var_class, MA.value, article_ID) %>% 
                summarize(ns = n()) %>%
                  group_by(var_class, MA.value) %>% 
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
              arrange(desc(nr_articles))

# we will use studies reporting variables measured by at least 5 articles
length(subset(resp_var_f0, nr_articles >= 5)$MA.value) # 5 variables have 5 articles or more

# filter the data to include only responses with more than 5 articles
dt_5_f0 <- dt_f0 %>% filter(MA.value %in% subset(resp_var_f0, nr_articles >= 5)$MA.value) 
       # here we do not really need to define grouping variables for plotting results of meta-analysis
       
model_guide <- dt_5_f0 %>% distinct(MA.value, var_class) %>% 
                rename(eco_response = MA.value) %>% 
                # create 'selection' formula (tells the loop which MA.value to take)
                mutate(select = paste0("MA.value == ", "'" , eco_response, "'"),
                       # add a unique ID for each model
                       model_id = paste(eco_response, "intercept_only",sep = "_"))

set.seed(161)

# loop models
# first create empty data tables to store model outputs
intercept.tmp <- data.table(var_class = NA, eco_response = NA, 
                            model_id = NA, model_call = NA,
                            estimate = NA, ci.lb = NA, ci.ub = NA, 
                            pi.lb = NA, pi.ub = NA, pval = NA)

eco.intercepts <- data.table(var_class = NA, eco_response = NA, 
                             model_id = NA, model_call = NA,
                             estimate = NA, ci.lb = NA, ci.ub = NA, 
                             pi.lb = NA, pi.ub = NA, pval = NA)

# and create an empty list to store forest and orchard plots
forest.plot_list <- list()
orchard.plot_list <- list()

# set the working directory to data (where the folder "extra_analyses">"builds_zero" is)
# make sure the folder "builds" contains the folders "models" and "null_models"
setwd("./data")

# run the models
tic()
for(i in 1:nrow(model_guide)){
  result <- tryCatch({ 
    # select subset of data for each response variable
    data.sub = dt_5_f0[eval(parse(text = model_guide[i, ]$select)),]
    # build multi-level meta-analysis (only intercept)
    m0 <- rma.mv(yi_smd ~ 1, # intercept only model
                 # add variance covariance matrix to account for statistical non-independence
                 V = vcalc(vi = vi_smd, # sampling variances that are correlated with each within the same study;
                           cluster = article_ID,  # study identity - clustering variable;
                           obs = study_ID, # different effect sizes corresponding to the same response/dependent variable
                           data = data.sub, 
                           rho = 0.5), 
                 random = list(~ 1 | article_ID / study_ID),
                 data = data.sub, 
                 method = "REML",  test = "t", dfs = "contain") # method has to be changed to ML when comparing models with likelihood tests 
    # extract model results and add to intercept.tmp
        intercept.tmp$var_class <- model_guide[i, ]$var_class
        intercept.tmp$eco_response <- model_guide[i, ]$eco_response
        intercept.tmp$model_id <- model_guide[i, ]$model_id
        intercept.tmp$estimate <- m0$b[1]
        intercept.tmp$ci.lb <- m0$ci.lb[1]
        intercept.tmp$ci.ub <- m0$ci.ub[1] 
        intercept.tmp$pi.lb <- predict(m0)$pi.lb[1]
        intercept.tmp$pi.ub <- predict(m0)$pi.ub[1]
        intercept.tmp$pval <- m0$pval[1] 
    intercept.tmp$model_call <- paste0("extra_analyses/builds_zero/models/null_models/", model_guide[i, ]$model_id,".rds")
    eco.intercepts <- rbind(eco.intercepts, intercept.tmp)
    
    # forest plot
    forest.plot <- forest(m0, xlab = "Hedges g", addpred = T, slab = study_ID)
                     title(model_guide[i, ]$eco_response)
    forest.plot_list[[i]] <- forest.plot # store each plot in our list 
    names(forest.plot_list)[i] <- model_guide[i, ]$eco_response
    
    # orchard plot
    orchard.plot <- orchard_plot(m0, mod = "1", 
             xlab = "Effect size (Hedges g)", legend.pos = "none",
             group = "article_ID", k = TRUE, g = TRUE) + 
             scale_x_discrete(labels = c("Overall effect")) +
             labs(title = gsub("_", " ", model_guide[i, ]$eco_response))
    orchard.plot_list[[i]] <- orchard.plot # store each plot in our list 
    names(orchard.plot_list)[i] <- model_guide[i, ]$eco_response
    
  }, error = function(e) {
    return(NULL)
  })
  if (is.null(result)) {
    next
  }
  write_rds(m0, paste0("extra_analyses/builds_zero/models/null_models/", model_guide[i, ]$model_id, ".rds"))
  cat(i,"/", nrow(model_guide), "\r")
}
toc()

model_master <- eco.intercepts[!is.na(eco.intercepts$eco_response)] # remove first row (all NAs)

sig <- model_master[pval < 0.05,]
sig # empty table -- no models are significant

fwrite(model_master, "extra_analyses/builds_zero/model_results/null_model_results_zero.csv")

# check the orchard plots
do.call(grid.arrange, c(orchard.plot_list[1:5], ncol=3, nrow = 2))

model_master <- model_master %>% 
                  ## add significance stars and direction
                  mutate(stars = case_when(pval < 0.05 & pval > 0.01 ~ "*",
                                          pval <= 0.01 & pval > 0.001 ~ "**",
                                          pval <= 0.001 ~ "***",
                                          TRUE ~ ""),
                         sig = case_when(pval < 0.05 & estimate < 0 ~ "Negative effect",
                                        pval < 0.05 & estimate > 0 ~ "Positive effect",
                                        pval > 0.05 ~ "Non-significant")) 

#rename 
dt_5_f0 <- dt_5_f0 %>% rename(eco_response = MA.value)

# create data table for sample sizes 
sample_sizes <- dt_5_f0 %>% group_by(var_class, eco_response, article_ID) %>% 
                           summarize(ns = n()) %>% group_by(var_class, eco_response) %>% 
                             summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                  # create a variable with nr of articles and studies for plotting
                  mutate(sample_size = paste(nr_articles, "(", nr_studies, ")", sep=""))
  
## plot
p <- ggplot()+
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed")+
  geom_jitter(data = dt_5_f0, # could subset data to dt[abs(yi_smd) < 6,]
              aes(x = yi_smd, y = gsub("_", " ", eco_response), 
                  size = 1/vi_smd, color = yi_smd), #size = 1/vi_smd indicates model weight (the bigger the more influential)
              height = .1, inherit.aes = F, alpha = 0.5) +
  scale_colour_gradient2(low = "#006699FF", mid = "#B3B8B3FF", high = "#B35900FF",
                  midpoint = 0, space = "Lab") +
  # add the model predictions
  new_scale_fill() +
  # add prediction intervals
  geom_linerange(data = model_master[], 
                 aes(x = estimate, xmin = pi.lb, xmax = pi.ub, 
                     y = gsub("_", " ", eco_response))) +
  # add predicted effect sizes and confidence intervals
  geom_pointrange(data = model_master[], 
                  aes(x = estimate, xmin = ci.lb, xmax = ci.ub,
                      fill = sig, y = gsub("_", " ", eco_response)),
                  shape = 21, size = 1, linewidth = 1.2, stroke = 0.5) +
  scale_fill_manual(NULL, values = c("Negative effect" = "#006699FF", 
                                     "Non-significant" = "#B3B8B3FF",
                                     "Positive effect" = "#B35900FF")) +
  # add asterisks beside axis names to indicate the ones that changed
  scale_y_discrete(labels = c("plant height" = "plant height*", 
                              "plant abundance total" = "plant abundance total*",
                              "plant abundance lichens" = "plant abundance lichens",
                              "plant abundance graminoids" = "plant abundance graminoids*",
                              "plant abundance dwarf shrubs" = "plant abundance dwarf shrubs")) +
  geom_text(data = sample_sizes, # indicate number of articles and studies
            aes(x = -18, y = gsub("_", " ", eco_response), label = sample_size),
            size = 3, color = "grey10") +
  geom_text(data = model_master[],
            aes(y = gsub("_", " ", eco_response), x= 6, label = stars), 
            color = "grey10", size = 5) +
  xlab("Effect size") + ylab(NULL) +
  xlim(-20, 15) + 
  guides(size = FALSE) +
  theme_systrev() +
  theme(legend.position = "bottom",
    legend.key.width = unit(2.5, "cm"),
    legend.box = "vertical", legend.margin = margin(),
    legend.key.height = unit(0.8, "lines"),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 9))
p 

ggsave(p, file= "builds_zero/plots/estimate_plots/intercept_plot.png", dpi = 600)
