#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Systematic review on herbivory diversity in Arctic tundra
#                 Isabel C Barrio (isabel@lbhi.is)
#                              1-July-2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# script to run some descriptive statistics for the full set of eligible 
# studies included in the systematic review

# libraries --------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyverse)
library(stringr) 
library(flextable)    # to make beautiful tables :)
library(harrypotter)  # colour palette


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



# 1. load datasets -------------------------------------------------------------
# clean the environment first
rm(list=ls())

# set the working directory to the folders where the data files are stored
setwd("./data")

# in this script we will use the file coded_data.csv (from script 1), which
# contains coded raw data extracted from studies, as well as the lists of
# herbivores, herbivore groups and plants (from script 1)

coded_data <- read_csv("coded_data.csv") 
herbivore_groups <- read_csv("herbivore_groups.csv") 
new_values_herbivores <- read_csv("new_values_herbivores.csv") 
new_values_plants <- read_csv("new_values_plants.csv") 


# 2. descriptive statistics ----------------------------------------------------
# for the first part of this section presented in the manuscript (i.e. ROSES 
# diagram and repeatability of screening) see script 2_screening_process

# range of publication years for the articles included?
coded_data %>% summarize(min = min(year), max = max(year)) 

# language of the articles included?
coded_data %>% group_by(language, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(language) %>% 
                        summarize(nr_articles = n(), nr_studies = sum(ns)) 

# geographical location
coded_data %>% group_by(country, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(country) %>% 
                        summarize(nr_articles = n(), nr_studies = sum(ns)) %>% 
                    arrange(desc(nr_articles))
# the number of articles is actually highest for Canada, but in that part
# of the text we are reporting number of studies


# management and conservation focus of the article?
coded_data %>% group_by(management_focus) %>%  
                  summarize(n = n(),
                            percentage = (n/length(coded_data$management_focus))*100)
coded_data %>% group_by(management_focus, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(management_focus) %>% 
                        summarize(nr_articles = n(), nr_studies = sum(ns)) 
70/201 # 34.8% of articles had a management focus

coded_data %>% group_by(conservation_focus) %>%  
                  summarize(n = n(),
                            percentage = (n/length(coded_data$conservation_focus))*100)
coded_data %>% group_by(conservation_focus, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(conservation_focus) %>% 
                        summarize(nr_articles = n(), nr_studies = sum(ns)) 
13/201 # 6.5% of articles had a conservation focus

# study design and diversity contrast?
coded_data %>% group_by(study_design) %>%  
                  # report by study here
                  summarize(n = n(),
                            percentage = (n/length(coded_data$study_design))*100)

# what kind of experimental design?
coded_data %>% filter(study_design == "experimental") %>%  
                  group_by(diversity_contrast, article_ID) %>%  
                  # report by study here
                  summarize(ns = n()) %>%
                      group_by(diversity_contrast) %>% 
                        summarize(nr_articles = n(), nr_studies = sum(ns))
coded_data %>% filter(study_design == "experimental" & diversity_contrast == "other") %>% 
                select(article_ID)


coded_data %>% group_by(study_design, diversity_contrast) %>%  
                  summarize(n = n())

# how many articles(studies) used size selective exclosures?
coded_data %>% filter(size_selective_exclosures == "yes") %>% 
                  group_by(size_selective_exclosures, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(size_selective_exclosures) %>% 
                        summarize(nr_articles = n(), nr_studies = sum(ns))

# we will make a table of these, to include in Additional file 4
sizeselective.table <- coded_data %>% 
                    filter(size_selective_exclosures == "yes") %>% 
                    select(article_ID, study_ID, author_list, year, country, locality,
                           herbivore_ID_higher_diversity, herbivore_ID_lower_diversity,
                           herbivore_fgr_ID_higher_diversity, herbivore_fgr_ID_lower_diversity,
                           contrast_f, change_f, change.long_f,
                           MA.value, measured_response_variable_new, quantitative_synthesis) %>% 
                    group_by(article_ID, author_list, year, locality, quantitative_synthesis) %>% summarize(n = n()) 

sizeselective.t <- flextable(sizeselective.table) %>% 
            merge_v(j = ~ `article_ID` + `author_list` + `year`) %>% 
            autofit() 
sizeselective.t

save_as_docx(sizeselective.t, path = "./tables/sizeselective_table.docx")


# how many articles(studies) studied invertebrate herbivory?
coded_data %>%
  filter(herbivore_fgr_ID_higher_diversity %in% "inv" | herbivore_fgr_ID_lower_diversity %in% "inv") %>%
  group_by(article_ID) %>% summarize(ns = n()) %>%
    ungroup() %>%  summarize(nr_articles = n(), nr_studies = sum(ns))

(18/201)*100 # 8.9% of articles considered invertebrates


# herbivore data?
coded_data %>% group_by(herbivore_data, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(herbivore_data) %>% 
                        summarize(nr_articles = n(), 
                                  percent = (nr_articles/length(unique(coded_data$article_ID)))*100, 
                                  nr_studies = sum(ns)) 


# ecological context reporting?
coded_data %>% group_by(soil_chemistry, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(soil_chemistry) %>% 
                        summarize(nr_articles = n(), 
                                  percent = (nr_articles/length(unique(coded_data$article_ID)))*100, 
                                  nr_studies = sum(ns)) 
((201-172)/201)*100 # reported for 14.4% of included articles

coded_data %>% group_by(soil_texture, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(soil_texture) %>% 
                        summarize(nr_articles = n(), 
                                  percent = (nr_articles/length(unique(coded_data$article_ID)))*100, 
                                  nr_studies = sum(ns)) 
((201-184)/201)*100 # reported for 8.5% of included articles

coded_data %>% group_by(soil_moisture, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(soil_moisture) %>% 
                        summarize(nr_articles = n(), 
                                  percent = (nr_articles/length(unique(coded_data$article_ID)))*100, 
                                  nr_studies = sum(ns)) 
((201-128)/201)*100 # reported for 36.3% of included articles

coded_data %>% group_by(permafrost, article_ID) %>%  
                    summarize(ns = n()) %>%
                      group_by(permafrost) %>% 
                        summarize(nr_articles = n(), 
                                  percent = (nr_articles/length(unique(coded_data$article_ID)))*100, 
                                  nr_studies = sum(ns)) 
((201-155)/201)*100 # reported for 22.9% of included articles




# 3. summary tables ------------------------------------------------------------
# define general formatting properties for tables using the package flextable
set_flextable_defaults(font.size = 10, theme_fun = theme_vanilla,
                       padding = 4, background.color = "#EFEFEF")
# init_flextable_defaults() # to restore defaults

## 3.1 outcome variable data ---------------------------------------------------
outcome.table <- coded_data %>% 
                    select(article_ID, study_ID, 
                           var_class, MA.value, measured_response_variable_new) %>% 
                    group_by(var_class, MA.value) %>% 
                      distinct(measured_response_variable_new) %>% # remove duplicate values
                      arrange(var_class, MA.value, measured_response_variable_new) %>% 
                     # remove underscores before printing table
                     mutate(MA.value = str_replace_all(MA.value, "_", " "),
                            measured_response_variable_new = str_replace_all(measured_response_variable_new, "_", " ")) %>% 
                     # rename columns
                     rename("outcome variable" = MA.value, 
                            "group" = var_class,
                            "reported outcome variable" = measured_response_variable_new) 

outcome.table %>% ungroup() %>% summarize(n= length(`reported outcome variable`)) 
  # 319 unique values for outcome variables
outcome.table %>% ungroup() %>% distinct(`outcome variable`) %>% 
  summarize(n= length(`outcome variable`)) 
  # 100 unique groups (MA.value)

outcome.t <- flextable(outcome.table) %>% 
            merge_v(j = ~ `group` + `outcome variable`) %>% 
            autofit() 
outcome.t

save_as_docx(outcome.t, path = "./tables/outcome_table.docx")


# a plot of outcome variables is provided in script 5


## 3.2 herbivore data ----------------------------------------------------------
# make a table of herbivore groups and which values of herbivore ID each one includes 

### 3.2.1 body size groups -----------------------------------------------------
herbivore.table <- new_values_herbivores %>% group_by(herbivore_group, type) %>% 
                     distinct(new.value) %>% # remove duplicate values
                     select(herbivore_group, type, new.value) %>% 
                     arrange(herbivore_group, desc(type), new.value) %>% 
                     # remove underscores before printing table
                     mutate(herbivore_group = str_replace_all(herbivore_group, "_", " "),
                            new.value = str_replace_all(new.value, "_", " ")) %>% 
                     filter(!herbivore_group == "zero") %>% # remove row for zero
                     # rename columns
                     rename("body size group" = "herbivore_group", 
                            "reported herbivore ID" = new.value)

herb.t <- flextable(herbivore.table) %>% 
            merge_v(j = ~ `body size group` + type) %>% 
            italic(i = ~ type == "species", j = "reported herbivore ID") %>% 
            autofit() 
herb.t

save_as_docx(herb.t, path = "./tables/herbivore_table.docx")

# herbivore diversity contrasts (body size)
herbivore.gr <- herbivore_groups %>% group_by(change.num, change.long, change, contrast, article_ID) %>% 
                  summarize(ns = n()) %>% group_by(change.num, change.long, change, contrast) %>% 
                    # get nr of articles and studies for each contrast
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% ungroup() %>% 
                  mutate(change.z = paste0(str_replace_all(change.long, "_", " "), " (",
                                           change, ")" ),
                         contrast = str_replace_all(contrast, "_", " "),
                         records = paste0(nr_articles, "(", nr_studies, ")")) %>% 
                  select(change.num, change.z, contrast, records) %>% 
                  # make prettier labels for table header
                  rename("numerical change" = change.num,
                         "identity of change" = change.z,
                         "herbivore diversity contrast" = contrast,
                         "nr of records" = records)

herb.contrasts.t <- flextable(herbivore.gr) %>% 
            merge_v(j = ~ `numerical change` + `identity of change`) %>% 
            align(j = c("numerical change", "nr of records"), align = "center") %>% 
            autofit() 
herb.contrasts.t 

save_as_docx(herb.contrasts.t, path = "./tables/herbivore_contrasts.docx")


### 3.2.2 functional groups ----------------------------------------------------
herbivore.f.table <- new_values_herbivores %>% group_by(herbivore_fct_group, type) %>% 
                     distinct(new.value) %>% # remove duplicate values
                     select(herbivore_fct_group, type, new.value) %>% 
                     arrange(herbivore_fct_group, desc(type), new.value) %>% 
                     # remove underscores before printing table
                     mutate(herbivore_fct_group = str_replace_all(herbivore_fct_group, "_", " "),
                            new.value = str_replace_all(new.value, "_", " ")) %>% 
                     filter(!herbivore_fct_group == "zero") %>% # remove row for zero
                     # rename columns
                     rename("functional group" = herbivore_fct_group, 
                            "reported herbivore ID" = new.value) 
                     

herb.f.t <- flextable(herbivore.f.table) %>% 
            merge_v(j = ~ `functional group` + type) %>% 
            italic(i = ~ type == "species", j = "reported herbivore ID") %>% 
            autofit() 
herb.f.t

save_as_docx(herb.f.t, path = "./tables/herbivore_f_table.docx")

# herbivore diversity contrasts
herbivore.f.gr <- herbivore_groups %>% group_by(change.num.f, change.long_f, change_f, contrast_f, article_ID) %>% 
                  summarize(ns = n()) %>% group_by(change.num.f, change.long_f, change_f, contrast_f) %>% 
                    # get nr of articles and studies for each contrast
                    summarize(nr_articles = n(), nr_studies = sum(ns)) %>% ungroup() %>% 
                  mutate(change.z_f = str_replace_all(change.long_f, "_", " "),
                         contrast = str_replace_all(contrast_f, "_", " "),
                         records = paste0(nr_articles, "(", nr_studies, ")")) %>% 
                  select(change.num.f, change.z_f, contrast_f, records) %>% 
                  # make prettier labels for table header
                  rename("numerical change" = change.num.f,
                         "identity of change" = change.z_f,
                         "herbivore diversity contrast" = contrast_f,
                         "nr of records" = records)

# which article has the contrast F2 | F3?
herbivore_groups %>% filter(contrast_f == "F2 | F3") %>% distinct(article_ID)

herb.contrasts.f.t <- flextable(herbivore.f.gr) %>% 
            merge_v(j = ~ `numerical change` + `identity of change`) %>% 
            align(j = c("numerical change", "nr of records"), align = "center") %>% 
            autofit() 
herb.contrasts.f.t 

save_as_docx(herb.contrasts.f.t, path = "./tables/herbivore_contrasts_f.docx")





## 3.3 plant data --------------------------------------------------------------
plant.table <- new_values_plants %>% group_by(plant_group, type) %>% 
                     distinct(new.value) %>% # remove duplicate values
                     select(plant_group, type, new.value) %>% 
                     arrange(plant_group, desc(type), new.value) %>% 
                     # remove underscores before printing table
                     mutate(plant_group = str_replace_all(plant_group, "_", " "),
                            new.value = str_replace_all(new.value, "_", " ")) %>% 
                     # rename columns
                     rename("plant group" = plant_group, 
                            "reported plant ID" = new.value) 

plant.t <- flextable(plant.table) %>% 
            merge_v(j = ~ `plant group` + type) %>% 
            italic(i = ~ type == "species", j = "reported plant ID") %>% 
            autofit() 
plant.t

save_as_docx(plant.t, path = "./tables/plant_table.docx")



  
