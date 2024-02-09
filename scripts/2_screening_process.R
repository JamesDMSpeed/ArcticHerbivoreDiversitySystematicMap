#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Systematic review on herbivory diversity in Arctic tundra
#                        Screening process
#                 Isabel C Barrio (isabel@lbhi.is)
#                         10-March-2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# script to provide an overview of the screening process (ROSES diagram)
# and calculate agreement between observers at different stages 
# of article screening, including assessment of risk of bias

# libraries --------------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(irr)         # to calculate Cohen's kappa


# load data --------------------------------------------------------------------
# clean the environment first
rm(list=ls())

# set the working directory to the folders where the data files are stored
setwd("./data")

# in this script we will use two files saved as csv after running script 1: 
# 1) overview.csv: contains bibliographic information about the articles and 
#    the screening process
# 2) coded_data.csv: contains coded raw data extracted from studies

overview <- read_csv("overview.csv")
coded_data <- read_csv("coded_data.csv") 

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


# 1. ROSES diagram -------------------------------------------------------------
# to build the upper part of the ROSES diagram we need information on number of 
# articles found from different sources and reasons for exclusion at each stage. 
# Here we get the nr of documents that we will then use to populate the upper
# part of the ROSES diagram; the lower part we will get in script 5.

# total number of documents
names(overview)
  length(levels(as.factor(overview$ID))) # 5944 articles

# total number of documents by source
overview %>%
  mutate(search_source = case_when(search_source %in% c("WoS", "Scopus") ~ "WoS_Scopus",
                                   search_source == "Google Scholar" ~ "Google Scholar",
                                   search_source == "snowballing" ~ "snowballing",
                                   TRUE ~ "Local search")) %>% 
  group_by(search_source) %>% 
  count()  # Google Scholar = 280, Local search = 413, snowballing = 70, WoS_Scopus = 5181

# how many articles retrieved through snowballing were included in the end?
overview %>% filter(search_source == "snowballing") %>% 
              filter(full_text_score_1 == "include") %>% count() # 14 articles

# number of articles included or excluded after title screening
score_title <- overview %>%
            group_by(title_score) %>% count()

# number of duplicate documents            
score_title[score_title$title_score == "duplicate", ]$n  #1997
  # number of documents after removing duplicates
  length(levels(as.factor(overview$ID))) - score_title[score_title$title_score == "duplicate", ]$n  #3947

# number of documents excluded at the title stage            
score_title[score_title$title_score == "exclude_title", ]$n  #1886
  # number of documents after title screening
  length(levels(as.factor(overview$ID))) - score_title[score_title$title_score == "duplicate", ]$n -
                                            score_title[score_title$title_score == "exclude_title", ]$n  #2061

# number of documents excluded at the abstract stage
score_abstract <- overview %>%
            group_by(abstract_score) %>% 
            count()

score_abstract[score_abstract$abstract_score == "exclude_abstract", ]$n  #1417
  # number of documents after abstract screening
  length(levels(as.factor(overview$ID))) - score_title[score_title$title_score == "duplicate", ]$n -
                                           score_title[score_title$title_score == "exclude_title", ]$n -
                                           score_abstract[score_abstract$abstract_score == "exclude_abstract", ]$n  #644

# number of unretrievable full texts            
score_fulltext <- overview %>%
            group_by(full_text_score) %>% 
            count()

score_fulltext[score_fulltext$full_text_score == "exclude_no_full_text", ]$n  #44
  # number of documents after removing non-retrievable documents
  length(levels(as.factor(overview$ID))) - score_title[score_title$title_score == "duplicate", ]$n -
                                           score_title[score_title$title_score == "exclude_title", ]$n -
                                           score_abstract[score_abstract$abstract_score == "exclude_abstract", ]$n -
                                           score_fulltext[score_fulltext$full_text_score == "exclude_no_full_text", ]$n  #600

# there are several categories that refer to unsuitable document types:
# corrections to published articles (indicated in column title_score = correction/reply), raw datasets, maps and supplementary materials
score_fulltext[score_fulltext$full_text_score == "exclude_other", ]$n  #37

# reasons for excluding at full text stage
score_fulltext[score_fulltext$full_text_score == "exclude_population_coordinates", ]$n  # 109
score_fulltext[score_fulltext$full_text_score == "exclude_population_habitat", ]$n      # 18
score_fulltext[score_fulltext$full_text_score == "exclude_exposure_herbivore", ]$n      # 13
score_fulltext[score_fulltext$full_text_score == "exclude_comparator", ]$n              # 146
score_fulltext[score_fulltext$full_text_score == "exclude_study_design", ]$n            # 42
  sum(score_fulltext[score_fulltext$full_text_score %in% 
                     c("exclude_outcome_effect","exclude_outcome"), ]$n)                # 20
score_fulltext[score_fulltext$full_text_score == "exclude_redundant", ]$n               # 14

# total number of studies excluded at full text stage
sum(score_fulltext[score_fulltext$full_text_score %in% 
                     c("exclude_population_coordinates","exclude_population_habitat",
                       "exclude_exposure_herbivore", "exclude_comparator", "exclude_study_design",
                       "exclude_outcome_effect", "exclude_outcome", "exclude_redundant",
                       "exclude_other"), ]$n)  #399

# number of documents after full text screening
length(levels(as.factor(overview$ID))) - 
  score_title[score_title$title_score == "duplicate", ]$n -
  score_title[score_title$title_score == "exclude_title", ]$n -
  score_abstract[score_abstract$abstract_score == "exclude_abstract", ]$n -
  score_fulltext[score_fulltext$full_text_score == "exclude_no_full_text", ]$n -
    sum(score_fulltext[score_fulltext$full_text_score %in% 
                     c("exclude_population_coordinates","exclude_population_habitat",
                       "exclude_exposure_herbivore", "exclude_comparator", "exclude_study_design",
                       "exclude_outcome_effect", "exclude_outcome", "exclude_other", 
                       "exclude_redundant"), ]$n) 
# this should be the same as 
score_fulltext[score_fulltext$full_text_score == "include", ]$n
  
  
  
# 2. repeatability of screening process ----------------------------------------
# here we want to calculate agreement between observers in the different stages
# of article screening

## 2.1 title screening ---------------------------------------------------------
# how many observers scored articles by title?
overview %>% 
  # remove duplicates
  filter(!title_score %in% c("duplicate")) %>% 
  select(ID, person_title_1, person_title_2) %>% 
  pivot_longer(cols=c(person_title_1, person_title_2), names_to = "person",
                      values_to = "name") %>% 
    drop_na() %>%  #drop NAs
    group_by(name) %>% count() %>% arrange(-n)

# how many titles were there without duplicates?  
sum(score_title[score_title$title_score != "duplicate", ]$n) # 3947 titles

title.scores <- overview %>% 
                  # remove duplicates and correction/reply
                  filter(search_source %in% c("Scopus", "WoS"),
                         !title_score %in% c("duplicate", "correction/reply")) %>% 
                  # select the columns with the title scores (either include or exclude_title)
                  select(exclusion_score_title_1, exclusion_score_title_2) %>% 
                  drop_na()
                  
(3520/3947)*100 #89% of titles were assessed by two reviewers

# percent agreement: number of agreements divided by the total number of screened articles
agree(title.scores)

# Cohen's kappa is a measure of agreement that takes into account
# that agreement may arise by chance. Cohen’s kappa also deals with situations 
# where observers use some of categories more than others (which affects the 
# calculation of how likely it is that observers agree by chance).
# Kappa results are interpreted as follows: 
#   - values ≤ 0 indicate no agreement 
#   - 0.01–0.20 none to slight agreement
#   - 0.21–0.40 fair agreement
#   - 0.41– 0.60 moderate agreement 
#   - 0.61–0.80 substantial agreement 
#   - 0.81–1.00 almost perfect agreement
kappa2(title.scores)
    #Kappa is 0.531 --> moderate agreement


## 2.2 abstract screening ------------------------------------------------------
# how many observers scored articles by abstract?
overview %>% 
  # include only articles included at title stage
  filter(title_score == "include") %>% 
  select(ID, person_abstract_1, person_abstract_2) %>% 
  pivot_longer(cols=c(person_abstract_1, person_abstract_2), names_to = "person",
                      values_to = "name") %>% 
    drop_na() %>%  #drop NAs  
  group_by(name) %>% count() %>%  arrange (-n)

abstract.scores <- overview %>% 
                      # include only articles included at title stage
                      filter(title_score == "include") %>% 
                      # select the columns with the abstract scores (either include or exclude_abstract)
                      select(exclusion_score_abstract_1, exclusion_score_abstract_2) %>% 
                      # drop abstracts not coded by two observers
                      drop_na()

# how many abstracts were screened by two observers? 
length(overview[overview$title_score == "include", ]$ID) #2050 articles included based on title
length(abstract.scores$exclusion_score_abstract_1) #204 abstracts coded by two observers
(204/2050)*100 #9.95 % of abstracts were screened by two observers

agree(abstract.scores)

kappa2(abstract.scores)
    #Kappa is 0.804 --> substantial agreement


## 2.3 full text screening ------------------------------------------------------
#how many observers scored articles by full text?
overview %>% 
  # include only articles included at abstract stage
  filter(abstract_score == "include") %>% 
  select(ID, person_full_text_1, person_full_text_2) %>% 
  pivot_longer(cols=c(person_full_text_1, person_full_text_2), names_to = "person",
                      values_to = "name") %>% 
    drop_na() %>%  #drop NAs
  group_by(name) %>% count() %>% arrange(-n)

full.text.scores <- overview %>% 
                      # include only articles included at abstract stage
                      filter(abstract_score == "include") %>% 
                      # select the columns with the abstract scores (either include or exclude_abstract)
                      select(exclusion_score_full_text_1, exclusion_score_full_text_2) %>% 
                      # drop full texts not coded by two observers
                      drop_na()

# how many abstracts were screened by two observers? 
length(na.omit(overview[overview$abstract_score == "include", ]$ID)) #633 articles included based on full text score
length(full.text.scores$exclusion_score_full_text_1) #28 full texts were coded by two observers
(28/633)*100 #4.4 % of full texts were screened by two observers

agree(full.text.scores) # 96.4 %

kappa2(full.text.scores)
    #Kappa is 0.92 --> near perfect agreement :)


# how many observers coded articles?
overview %>% 
  # include only articles included at title stage
  filter(full_text_score == "include") %>% 
  select(ID, coded_person) %>% 
  group_by(coded_person) %>% count() %>% arrange(-n)

# how many observers did snowballing?
overview %>% 
  # remove duplicates and correction/reply
  filter(snowballing == "yes") %>% 
    select(ID, snowballing_person) %>% 
      group_by(snowballing_person) %>% count()


# the last step from the upper part of the ROSES diagram is to report the number 
# of studies included in the systematic review dataset. We have this information
# in the other dataset (coded_data):
coded_data %>% summarize(n = length(study_ID)) # 3713 studies

# check that the number of included articles is the same in both databases :)
coded_data %>% summarize(n = length(unique(article_ID))) # 201 articles

# how many studies per article?
coded_data %>% group_by(article_ID) %>% 
                summarize(nr.studies = n()) %>% 
               ungroup() %>% 
                summarize(min = min(nr.studies),
                          max = max(nr.studies),
                          mean = mean(nr.studies),
                          SE = sd(nr.studies) / sqrt(length(nr.studies))) 



# 3. study validity assessment -------------------------------------------------
# here we calculate the agreement between observers during study validity assessment

# study validity was assessed for each study, so we need to check first 
# if all studies in an article have the same assessment

diff_studies_list <- coded_data %>%
  # filter the original dataset to keep only the rows where article_ID 
  # is in articles_with_diff_studies
  #filter(article_ID %in% articles_with_diff_studies) %>%
  # group the dataset by article_ID and study_ID to check the scores for risk of bias
  group_by(article_ID, study_ID) %>%
  # create a variable with a combination of values for all criteria
  # for each article_ID and study_ID
  mutate(combinations = paste(bias_risk_criterion1,
                              bias_risk_criterion2,
                              bias_risk_criterion3or4,
                              bias_risk_criterion5,
                              bias_risk_criterion6,
                              bias_risk_criterion7,
                              overall_criterion, sep = "; ")) %>% 
  # just to simplify, select only relevant columns :)
  select(article_ID, study_ID, combinations) %>% 
  ungroup() %>% group_by(article_ID) %>%
  # find how many distinct combinations there are per article
  summarize(num_combinations = n_distinct(combinations)) %>% 
  # find the rows where the number of unique combinations is greater than 1
  filter(num_combinations > 1) %>%
  select(article_ID) %>%
  arrange(article_ID)

# print the list of articles with different bias risk criteria
if (nrow(diff_studies_list) == 0) {
  cat("All articles have the same score for all bias risk criteria.\n")
} else {
  cat("The following articles have different scores for bias risk criteria:\n")
  print(diff_studies_list)
}

# 8 articles have two combinations of scores of risk of bias, so we have 
# to report results per study, not pooled per article
study_validity <- coded_data %>% 
  select(article_ID, study_ID, 
         # scores given by the first reviewer
         bias_risk_criterion1, bias_risk_criterion2,
         bias_risk_criterion3or4, bias_risk_criterion5,
         bias_risk_criterion6, bias_risk_criterion7, overall_criterion,
         # scores given by the second reviewer
         bias_risk_criterion1_2, bias_risk_criterion2_2, 
         bias_risk_criterion3or4_2, bias_risk_criterion5_2,
         bias_risk_criterion6_2, bias_risk_criterion7_2, overall_criterion_2) %>% 
  # keep only studies assessed by both reviewers
  drop_na()

study_validity %>% summarize(nr.articles = length(unique(article_ID)),
                             nr.studies = length(study_ID))
(856/3713)*100 # 856 studies (23.1% of the database)
(21/201)*100 # 21 articles (10.5% of the database)

## 3.1 risk of bias 1 ----------------------------------------------------------
study_validity_1 <- study_validity %>% 
                      select(bias_risk_criterion1, bias_risk_criterion1_2)
agree(study_validity_1) # 76.2% agreement
kappa2(study_validity_1) # Cohen's kappa = 0.46 --> moderate agreement

## 3.2 risk of bias 2 ----------------------------------------------------------
study_validity_2 <- study_validity %>% 
                      select(bias_risk_criterion2, bias_risk_criterion2_2)
agree(study_validity_2)  # 76.9% agreement
kappa2(study_validity_2) # Cohen's kappa = 0.55 --> moderate agreement

## 3.3 risk of bias 3 or 4 -----------------------------------------------------
study_validity_3or4 <- study_validity %>% 
                         select(bias_risk_criterion3or4, bias_risk_criterion3or4_2)
agree(study_validity_3or4)  # 87.9% agreement
kappa2(study_validity_3or4) # Cohen's kappa = 0.67 --> substantial agreement

## 3.4 risk of bias 5 ----------------------------------------------------------
study_validity_5 <- study_validity %>% 
                      select(bias_risk_criterion5, bias_risk_criterion5_2)
agree(study_validity_5)  # 75.6% agreement
kappa2(study_validity_5) # Cohen's kappa = 0.49 --> moderate agreement

## 3.5 risk of bias 6 ----------------------------------------------------------
study_validity_6 <- study_validity %>% 
                      select(bias_risk_criterion6, bias_risk_criterion6_2)
agree(study_validity_6)  # 76.2% agreement
kappa2(study_validity_6) # Cohen's kappa = 0.48--> moderate agreement

## 3.6 risk of bias 7 ----------------------------------------------------------
study_validity_7 <- study_validity %>% 
                      select(bias_risk_criterion7, bias_risk_criterion7_2)
agree(study_validity_7)  # 83.1% agreement
kappa2(study_validity_7) # Cohen's kappa = 0.44 --> moderate agreement

## 3.7 overall risk of bias ----------------------------------------------------
study_validity_overall <- study_validity %>% 
                            select(overall_criterion, overall_criterion_2)
agree(study_validity_overall)  # 89.4% agreement
kappa2(study_validity_overall) # Cohen's kappa = 0.67 --> substantial agreement


## 3.8 risk of bias plot -------------------------------------------------------
# this is Figure 4 in the manuscript
# preparing the dataset
study_validity_plot <- coded_data %>% 
                         select(article_ID, study_ID, 
                         # scores for risk of bias
                         bias_risk_criterion1, bias_risk_criterion2,
                         bias_risk_criterion3or4, bias_risk_criterion5,
                         bias_risk_criterion6, bias_risk_criterion7, overall_criterion) %>%
                         # change the name of the different categories of overall criterion
                       mutate(overall_criterion = case_when(
                            overall_criterion == "Overall low risk of bias" ~ "low_risk",
                            overall_criterion == "Overall medium risk of bias" ~ "medium_risk",
                            overall_criterion == "Overall high risk of bias" ~ "high_risk")) %>%  
                      # select relevant columns (do not need article_ID or study_ID)
                      select(bias_risk_criterion1, bias_risk_criterion2,
                             bias_risk_criterion3or4, bias_risk_criterion5,
                             bias_risk_criterion6, bias_risk_criterion7,
                             overall_criterion) %>% 
                      # count how many studies had a certain score of risk of bias for each criterion
                      mutate(across(everything(), factor, levels = c("high_risk", "medium_risk", "low_risk"))) %>% 
                      pivot_longer(everything(), names_to = "criterion", values_to = "level") %>%
                        count(criterion, level)

# some descriptive stuff
study_validity_plot %>% filter(criterion == "overall_criterion") %>% 
    mutate(percentage = (n/3713)*100)

study_validity_plot %>% filter(level == "high_risk")

# preparing the plot
# define the colors for each level
level_colors <- c("low_risk" = "#006699FF", "medium_risk" = "#98C2D9FF", "high_risk" = "#B35900FF")

# define the labels for each criterion
new_names <- c("bias_risk_criterion1" = "Risk of confounding bias (criterion 1)",
               "bias_risk_criterion2" = "Risk of post-intervention sampling bias (criterion 2)",
               "bias_risk_criterion3or4" = "Risk of misclassified comparison (criterion 3)
               or performance bias (criterion 4)",
               "bias_risk_criterion5" = "Risk of measurement bias (criterion 5)",
               "bias_risk_criterion6" = "Risk of outcome reporting bias (criterion 6)",
               "bias_risk_criterion7" = "Risk of outcome assessment biases (criterion 7)",
               "overall_criterion" = "Overall risk")

# change the order of the criteria because we are going to flip the coordinates in the plot
desired_order <- c("overall_criterion", "bias_risk_criterion7", 
                   "bias_risk_criterion6", "bias_risk_criterion5",
                   "bias_risk_criterion3or4", "bias_risk_criterion2",
                   "bias_risk_criterion1")
study_validity_plot$criterion <- factor(study_validity_plot$criterion, levels = desired_order)

# plot
p <- ggplot(study_validity_plot, aes(x = criterion, y = n, fill = level)) +
     geom_col(position = "stack") +
     scale_fill_manual(values = level_colors,
                       labels = c("low_risk" = "low risk", 
                               "medium_risk" = "medium risk", 
                               "high_risk" = "high risk")) +
     guides(fill = guide_legend(title = NULL, reverse = TRUE)) +
        # reverse = TRUE, reverses the fill legend (colours) so low risk is the first one
     scale_x_discrete(labels = new_names) +
     labs(x = "", y = "Number of studies", fill = "Risk Level") +
     coord_flip() +
     theme_systrev() +
     theme(legend.position = "bottom",
           axis.text = element_text(size = 12, color="black"))
p

ggsave(p, file = "./figures/risk_of_bias.png", dpi = 600)

