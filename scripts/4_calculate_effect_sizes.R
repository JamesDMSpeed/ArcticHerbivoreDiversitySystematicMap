#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       Systematic review on herbivory diversity in Arctic tundra
#                       Calculating effect sizes
#     Jonas Trepel (jonas.trepel@gmail.com / jonas.trepel@bio.au.dk)
#                 15-April-2023, UPDATED: 19-Sept-2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# script to calculate effect sizes from studies retrieved in the systematic review 

# libraries --------------------------------------------------------------------
library(data.table)
library(metafor)
library(estmeansd)
library(tidyverse)
library(tidyr)
library(esc)      # to calculate effect sizes (SMD) from one-way ANOVA
library(metagear) # to impute missing SD


# load data --------------------------------------------------------------------

rm(list=ls())

# in this script we will use the file coded_data.csv (from script 1), which
# contains coded raw data extracted from studies

# make sure we load the most recent data set
# path name here is a path to an existing folder in your working directory
file <- list.files("main_data", pattern = '.csv',
                   full.names = T)

if(length(file) == 1){
  dt <- fread(file)
  print("YAY")
}else{
  print("CLEAN FILE FOLDER")
}

names(dt)

# set the working directory
# ICB: setwd("C:/Users/isabel/OneDrive - Menntaský/TUNDRAsalad systematic review/R/data")
# LBP: setwd("C:/Users/laura_jndxf5s/OneDrive - Menntaský/systematic review/R/data")
# dt <- fread("coded_data.csv")


## build scratch file ----------------------------------------
scratch <- dt  %>% 
  # select only the columns of interest   
  select("article_ID", "study_ID", 
         "measured_response_variable", "MA.value", "outcome_type", 
         "value_higher_diversity", "variability_higher_diversity",
         "value_lower_diversity", "variability_lower_diversity", "value_type", 
         "variability_type", "sample_size_higher_diversity", "sample_size_lower_diversity",
         "effect_size", "effect_size_type", "effect_size_variability", "effect_size_variability_type",
         "effect_size_direction", "effect_size_sample_size", "effect_size_comments", 
         "statistical_test", "statistical_test_value", "statistical_test_df", "statistical_test_p",
         "outcome_comments", "more_data") %>% 
  # for studies reporting IQR, range and some CIs we have a range of values for variability
  # keep those separate from studies reporting one value of variability
  mutate(variability_higher_diversity = ifelse(variability_higher_diversity == "NA", NA, variability_higher_diversity),
         HIGH_variability_range = ifelse(is.na(suppressWarnings(as.numeric(variability_higher_diversity))), 
                                    ifelse(is.na(variability_higher_diversity), NA, variability_higher_diversity), NA),
            # make sure that variability is numerical
            variability_higher_diversity = suppressWarnings(as.numeric(variability_higher_diversity)),
         variability_lower_diversity = ifelse(variability_lower_diversity == "NA", NA, variability_lower_diversity),
         LOW_variability_range = ifelse(is.na(suppressWarnings(as.numeric(variability_lower_diversity))), 
                                    ifelse(is.na(variability_lower_diversity), NA, variability_lower_diversity), NA),
            # make sure that variability is numerical
            variability_lower_diversity = suppressWarnings(as.numeric(variability_lower_diversity)),
         # make sure sample sizes are numeric
         sample_size_higher_diversity = suppressWarnings(as.numeric(sample_size_higher_diversity)),
         sample_size_lower_diversity = suppressWarnings(as.numeric(sample_size_lower_diversity)),
         # make sure values of high and low are numeric
         value_higher_diversity = suppressWarnings(as.numeric(value_higher_diversity)),
         value_lower_diversity = suppressWarnings(as.numeric(value_lower_diversity))) %>% 
  # need to split ranges into upper and lower limits
  separate(HIGH_variability_range, c("H_lower", "H_upper"), "; ") %>% 
  separate(LOW_variability_range, c("L_lower", "L_upper"), "; ") %>% 
  mutate(H_lower = as.numeric(H_lower),
         H_upper = as.numeric(H_upper),
         L_lower = as.numeric(L_lower),
         L_upper = as.numeric(L_upper)) 
summary(scratch) # scratch contains all observations

scratch %>% distinct(article_ID)

# some checks: do we have studies with SD or SE with a sample size of 1?
scratch %>% filter(variability_type == "SD" | variability_type == "SE") %>% 
            filter(!is.na(variability_higher_diversity) & !is.na(variability_lower_diversity)) %>% 
            filter(sample_size_higher_diversity == 1 | sample_size_lower_diversity == 1)
# there is only one study (51_e) -- checked it and they reported SE with a n=1 :( 

# LBP: scratch <- as.data.table(scratch)


## split by cases ---------------------------------------------------------------------------------
# split by cases depending on what measure of variability the studies report
# we will remove studies with NAs in the key variables (mean, variability, sample size)
# but we need to do that after splitting because variability for CIs and IQRs
# is stored in different columns :)
scratch %>% distinct(variability_type)

scratch <- scratch %>% 
            # create an additional variable for cases when we need to estimate variability
            mutate(error_type = case_when(variability_type %in% c("SE", "SD") ~ "reported",
                                          TRUE ~ variability_type))

  simple <- scratch[variability_type %in% c("SE", "SD")]           # 2725 studies
  CIs <- scratch[variability_type == "CI", ]                       # 178 studies
  IQRs <- scratch[variability_type == "IQR", ]                     # 255 studies
  range <- scratch[variability_type == "range", ]                  # 11 studies
  no.variability <- scratch %>% filter(is.na(variability_type))    # 544 studies
      no.variability %>% group_by(outcome_type) %>% summarize(n = n())
      # from these there are 177 reporting effect sizes and 17 reporting outcomes of statistical tests
      # and 350 reporting raw values -- let's check these
            no.variability %>% filter(outcome_type == "raw_value") %>% 
              group_by(value_type) %>% summarize(n = n())
            # 288 studies report a mean without variability
            # 62 studies report only one value ("total") for high and one for low diversity
            # with no associated variability measure
      # create separate subsets for using later
      no.and.mean <- no.variability[value_type %in% c("mean", "total")] # 350 studies
      no.no.mean <- no.variability[!value_type %in% c("mean", "total")] # 194 studies (these are the ones reporting results from 
                                                                        # statistical tests or effect sizes)

      # add columns to merge later
      no.and.mean[, error_type := "estimated"] # call this "estimated" if doing imputations (see section 2)
        no.and.mean[, `:=` (Mean_High_Diversity = value_higher_diversity,
                            Mean_Low_Diversity = value_lower_diversity)]
        no.and.mean[, `:=` (SD_High_Diversity = variability_higher_diversity,
                            SD_Low_Diversity = variability_lower_diversity)]
      summary(no.and.mean)

# just a quick check -- this should be zero :)
length(scratch$study_ID) - (length(simple$study_ID) + length(CIs$study_ID) +
                              length(IQRs$study_ID) + length(range$study_ID) + 
                              length(no.and.mean$study_ID) +
                              length(no.no.mean$study_ID))


# 1. convert to mean and SD ----------------------------------------------------
### 1.1 simple cases -----------------------------------------------------------
## let's start with the simple ones  :)
# variability of 0?
simple %>% filter(variability_higher_diversity == 0 | variability_lower_diversity == 0) %>% 
    distinct(study_ID) 
    # 82 studies report variability of zero (these were double-checked and actually reported 0)

# for each study we need mean and SD; for studies that report these values directly, 
# we store them in new columns that we will use in the meta-analyses
simple[value_type == "mean", 
       `:=` (Mean_High_Diversity = value_higher_diversity,
             Mean_Low_Diversity = value_lower_diversity) ]
simple[variability_type == "SD", 
       `:=` (SD_High_Diversity = variability_higher_diversity,
             SD_Low_Diversity = variability_lower_diversity)]

# for studies reporting mean and SE we need to convert SE to SD
simple %>% group_by(variability_type) %>% summarize(n = n()) # 2109 studies

simple[variability_type == "SE", 
       `:=` (SD_High_Diversity = variability_higher_diversity * sqrt(sample_size_higher_diversity),
             SD_Low_Diversity = variability_lower_diversity * sqrt(sample_size_lower_diversity))]

# there are some missing values for variability (98 studies)
simple %>% filter(is.na(simple$SD_Low_Diversity) | is.na(simple$SD_High_Diversity)) %>% distinct(study_ID)
# we mark these studies as "estimated" for error_type if we are doing imputations
simple <- simple %>% mutate(error_type = ifelse(is.na(simple$SD_Low_Diversity) | is.na(simple$SD_High_Diversity),
                                                "estimated", error_type))


### 1.2 studies reporting IQRs -------------------------------------------------
# we will use the package estmeansd to estimate mean and SD from medians and IQRs
# this is controversial because it only works for ca normally distributed data
# so we will check later if this transformation affects the results
out <- c()
    for(i in 1:nrow(IQRs)){
      out <- qe.mean.sd(q1.val = IQRs$H_lower[i],
                        med.val = IQRs$value_higher_diversity[i],
                        q3.val = IQRs$H_upper[i],
                        n = IQRs$sample_size_higher_diversity[i])
      IQRs[i, `:=` (Mean_High_Diversity = out$est.mean,
                    SD_High_Diversity = out$est.sd)]
      out <- qe.mean.sd(q1.val = IQRs$L_lower[i],
                        med.val = IQRs$value_lower_diversity[i],
                        q3.val = IQRs$L_upper[i],
                        n = IQRs$sample_size_lower_diversity[i])
      IQRs[i, `:=` (Mean_Low_Diversity = out$est.mean,
                    SD_Low_Diversity = out$est.sd)]
      print(i)
    }
# the function stores the estimated mean and SD in the IQRs dataset
summary(IQRs)


### 1.3 studies reporting 95% CIs ----------------------------------------------
# most studies report one value for CI --> in those cases we assume CIs are
# centered around the mean; for studies reporting two values, we use lower and upper CIs 

# we will use the following formula to calculate SD from CI 
# sqrt(n) * (upper limit - lower limit) / qt(1-0.05/2, df = n-1) 
# https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm

CIs <- CIs %>%
  mutate(HIGH_lower_CI = ifelse(is.na(variability_higher_diversity), H_lower,
                                 value_higher_diversity - variability_higher_diversity),
         HIGH_upper_CI = ifelse(is.na(variability_higher_diversity), H_upper,
                                 value_higher_diversity + variability_higher_diversity),
         LOW_lower_CI = ifelse(is.na(variability_lower_diversity), L_lower,
                                 value_lower_diversity - variability_lower_diversity),
         LOW_upper_CI = ifelse(is.na(variability_lower_diversity), L_upper,
                                 value_lower_diversity + variability_lower_diversity)) %>% 
  mutate(Mean_High_Diversity = case_when(value_type == "mean" ~ value_higher_diversity),
         Mean_Low_Diversity = case_when(value_type == "mean" ~ value_lower_diversity),
         SD_High_Diversity = sqrt(sample_size_higher_diversity) * 
           (HIGH_upper_CI - HIGH_lower_CI)/qt(1-0.05/2, df = (sample_size_higher_diversity-1)),
         SD_Low_Diversity = sqrt(sample_size_lower_diversity) * 
           (LOW_upper_CI - LOW_lower_CI)/qt(1-0.05/2, df = (sample_size_lower_diversity-1))) %>% 
  # remove extra columns for merging later
  select(-HIGH_lower_CI, -HIGH_upper_CI, -LOW_lower_CI, -LOW_upper_CI)

summary(CIs)


### 1.4 studies reporting ranges -----------------------------------------------
# SD can be estimated by quartering the range ( 1/4 * (Max − Min) )
# see: https://environmentalevidencejournal.biomedcentral.com/articles/10.1186/s13750-023-00299-x
range <- range %>%
  mutate(Mean_High_Diversity = case_when(value_type == "mean" ~ value_higher_diversity),
         Mean_Low_Diversity = case_when(value_type == "mean" ~ value_lower_diversity),
         SD_High_Diversity = 1/4*(H_upper - H_lower), 
         SD_Low_Diversity = 1/4*(L_upper - L_lower)) 
summary(range)



# 2. recombine all measures ----------------------------------------------------
# to deal with missing values for SD we can impute those values
# following Nakagawa et al (2023) via the impute_SD() function in metagear.

# this method yields some outliers when calculating effect sizes
# that we will remove (beyond some arbitrary threshold)
# but we will also run analyses with and without the imputed values

# for the imputations we take the values of simple and CI 
# (as those are directly reported by the studies; exclude IQR and ranges), 
# and impute the missing values in simple and no.and.mean

responses.imp <- rbind(simple, CIs, no.and.mean, fill = T)
      
responses.imp %>%  filter(is.na(SD_High_Diversity) | is.na(SD_Low_Diversity)) %>% 
  summarize(n = n()) # 448 rows need imputing missing values for SD (this are the 
                     # 350 from no.and.mean and the 98 missing one or other)

# all of them are recorded as error_type = "estimated"
responses.imp %>%  filter(is.na(SD_High_Diversity) | is.na(SD_Low_Diversity)) %>% distinct(error_type)

# the imputations should be done for each response variable separately
# we will do that with a loop. We create first a list of data, split 
# by response variable, to run the loop on each element of that list
variable_list <- split(responses.imp, responses.imp$MA.value)

# create a dataframe to store output
responses.impSD <- responses.imp[1,] 
  responses.impSD[1,] <- NA

set.seed(161)

for(i in 1:length(variable_list)){
  # tryCatch allows running the loop to the end even if we have errors
  # (here because we cannot do the imputations for all response variables
  # because there are some that have no SD values at all)
  result <- tryCatch({ 
        responses.tmp <- impute_SD(variable_list[[i]], 
                               columnSDnames = c("SD_High_Diversity", "SD_Low_Diversity"),    # columns with missing SDs
                               columnXnames = c("Mean_High_Diversity", "Mean_Low_Diversity"), # columns with mean values for each SD
                               method = "HotDeck")
        responses.impSD <- rbind(responses.impSD, responses.tmp)
        }, error = function(e) {
        return(NULL)
      })
      if (is.null(result)) {
        next
      }
}
  
# the loop drops some variables for which it was not possible to calculate SDs (5 studies)
# we will add them again (with missing values for SD), so that we do not miss those studies
missing.rows <- responses.imp %>% filter(!study_ID %in% responses.impSD$study_ID)

imputed.dataset <- rbind(responses.impSD, missing.rows)[-1,] #remove first row (all NAs)

# quick check -- the 448 imputed SD values are recorded as error_type = "estimated"
imputed.dataset %>% group_by(error_type) %>% summarize(n = n())

# paste everything together
responses <- rbind(imputed.dataset, IQRs, range, fill = T)


scratch[!study_ID %in% responses$study_ID] # 194 studies --> these are the ones 
                                           # reporting effect sizes or statistical tests
                                           # (the ones in no.no.mean)

summary(responses) # there should be only 5 NAs for SDs (the studies for which 
                   # we could not impute SD values)



# 3. calculate Hedges g (SMD) --------------------------------------------------
### 3.1 from studies for which we have mean, sd and n for both groups----
responses.out <- escalc(measure = "SMD", # standardized mean difference (Hedges g)
                        m1i = Mean_High_Diversity,
                        m2i = Mean_Low_Diversity,
                        sd1i = SD_High_Diversity,
                        sd2i = SD_Low_Diversity,
                        n1i = sample_size_higher_diversity,
                        n2i = sample_size_lower_diversity,
                        var.names=c("yi_smd","vi_smd"),
                        data = responses)
responses.out %>% filter(is.na(yi_smd) | is.na(vi_smd)) %>% summarize(n = n())
  # effect size cannot be calculated for 115 studies

# now we will calculate effect sizes from other sources: t-tests, one-way
# ANOVAs and reported effect sizes


### 3.2 from studies reporting t-test----
no.no.mean %>% group_by(outcome_type, statistical_test) %>% 
  summarize(n = n()) # there is one study reporting paired t-test that we cannot use

t_test <- no.no.mean %>% 
  filter(outcome_type == "statistical_test" & statistical_test == "t-test") %>%
  # statistical_test_value is read as character, change to numeric
  mutate(statistical_test_value = as.numeric(statistical_test_value))

responses.out_t <- escalc(measure   = "SMD", # standardized mean difference (Hedges g)
                          ti        = statistical_test_value,
                          n1i       = sample_size_higher_diversity,
                          n2i       = sample_size_lower_diversity,
                          var.names = c("yi_smd","vi_smd"),
                          data      = t_test) %>% 
                  # add the missing columns before we merge
                  mutate(error_type = "t",
                  Mean_High_Diversity = NA,
                  Mean_Low_Diversity = NA,
                  SD_High_Diversity = NA,
                  SD_Low_Diversity = NA)

### 3.3 from studies reporting one-way ANOVA----
responses.out_a <- no.no.mean %>% 
  filter(outcome_type == "statistical_test" & statistical_test == "one-way ANOVA") %>%
  # statistical_test_value is read as character, change to numeric
  mutate(statistical_test_value = as.numeric(statistical_test_value)) %>% 
  # add the missing columns before we merge
  mutate(error_type = "ANOVA",
         Mean_High_Diversity = NA,
         Mean_Low_Diversity = NA,
         SD_High_Diversity = NA,
         SD_Low_Diversity = NA)

# to calculate effect sizes we have to use a loop for each row in the dataframe
for (i in 1:nrow(responses.out_a)) {
  f <- responses.out_a$statistical_test_value[i] # get values from the current row
    grp1n <- responses.out_a$sample_size_higher_diversity[i]
    grp2n <- responses.out_a$sample_size_lower_diversity[i]
    
  # call the esc_f function to calculate effect size
  es_output <- esc_f(f = f,         # F value of the one-way anova
                     grp1n = grp1n, # sample size of group 1 
                     grp2n = grp2n, # sample size of group 2
                     es.type = "d") # SMD
  
  # extract the effect size and variance values
  vi_smd <- es_output$es
  yi_smd <- es_output$var
  
  # update the anovas dataframe with the new values
  responses.out_a$yi_smd[i] <- yi_smd
  responses.out_a$vi_smd[i] <- vi_smd
}             


### 3.4 from studies reporting effect sizes----
no.no.mean %>% filter(outcome_type == "effect_size") %>% 
  group_by(effect_size_type) %>%  summarize(n = n())
  # from these we can only transform to Hedges' g the study reporting change as
  # difference in means (= Cohen's d) and pooled variability; 
  # the others (171 studies) we have to leave out :(
scratch %>% filter(outcome_type == "effect_size") %>%
  filter(!effect_size_type == "change (grubbed vs control)") %>% summarize(n = n())

# which articles are these?
no.no.mean %>% filter(outcome_type == "effect_size") %>% 
  distinct(article_ID)

effect_sizes <- no.no.mean %>% filter(outcome_type == "effect_size") %>%
                            filter(effect_size_type == "change (grubbed vs control)") %>% 
                            mutate(ES = effect_size,
                                   SD_ES = as.numeric(effect_size_variability))

responses.out_es <- escalc(measure   = "SMD", # standardized mean difference (Hedges g)
                           di        = ES / SD_ES, # Cohen's d: change/ pooled SD
                           n1i       = sample_size_higher_diversity,
                           n2i       = sample_size_lower_diversity,
                           var.names = c("yi_smd", "vi_smd"),
                           data      = effect_sizes) %>% 
                  # add the missing columns before we merge
                  mutate(error_type = "ES",
                         Mean_High_Diversity = NA,
                         Mean_Low_Diversity = NA,
                         SD_High_Diversity = NA,
                         SD_Low_Diversity = NA) %>% 
                  select(-"ES", -"SD_ES")

         
# 4. add rows to the dataset----------------------------------------------------
responses.out <- rbind(responses.out, responses.out_t, responses.out_a, responses.out_es)

    # nr of studies with no effect size
    responses.out %>% filter(is.na(yi_smd)) %>% summarize(n = n()) # 115 studies

# yi is the effect size and vi the corresponding variance 
# (which will be used to weigh the data points in the meta-analyses)

setDT(responses.out)
summary(responses.out)
  quantile(responses.out$yi_smd, c(.005, .995), na.rm = TRUE)
  quantile(responses.out$vi_smd, c(.005, .995), na.rm = TRUE)


# some effect sizes could not be calculated
no.eff.size <- responses.out %>% filter(is.na(vi_smd) | is.na(yi_smd)) %>%  # 115 observations
  select(study_ID, Mean_High_Diversity, Mean_Low_Diversity, SD_High_Diversity, SD_Low_Diversity,
         sample_size_higher_diversity, sample_size_lower_diversity, error_type)
# this happens when studies report variability = 0 or sample size = 1
# for any of the treatments (high vs low), or for those where variability is missing

no.eff.size %>% filter(sample_size_higher_diversity == 1 | sample_size_lower_diversity == 1) # 70 studies
no.eff.size %>% filter(SD_High_Diversity == 0 | SD_Low_Diversity == 0) # 50 studies
no.eff.size %>% filter(is.na(SD_High_Diversity) | is.na(SD_Low_Diversity)) # 5 studies
no.eff.size %>% filter(sample_size_higher_diversity == 1 | sample_size_lower_diversity == 1 | 
                         SD_High_Diversity == 0 | SD_Low_Diversity == 0 | is.na(SD_High_Diversity) | is.na(SD_Low_Diversity)) # 115 studies

## there are some extreme outliers, most of them estimated
responses.out[yi_smd < quantile(responses.out$yi_smd, .005, na.rm = TRUE), ] # 18 outliers
responses.out[yi_smd > quantile(responses.out$yi_smd, .995, na.rm = TRUE), ] # 18 outliers

responses.out %>% filter(yi_smd < quantile(responses.out$yi_smd, .005, na.rm = TRUE) | 
                           yi_smd > quantile(responses.out$yi_smd, .995, na.rm = TRUE)) %>% 
                              summarize(n = n()) # 36 studies

range(responses.out$yi_smd, na.rm = TRUE) # ranges from -6893.7 to 11841.9

# write out as a side-file which can be added to the main dataset when needed
effect.sizes <- responses.out %>% 
                  # remove NAs
                  filter(!is.na(yi_smd)) %>%   
                  select(study_ID, yi_smd, vi_smd, error_type)
fwrite(effect.sizes, "side_files/effect_sizes.csv")

dt.smd <- left_join(dt, effect.sizes, by = "study_ID") %>% 
  # make sure that variables contributing "negatively" have a negative ES
  mutate(yi_smd = ifelse(measured_response_variable_new %in% c(
    "extinction", "mortality", "species_lost", "cellulose", "lignin",
    "ash_content", "shoot_death", "shoot_dieback", "leaf_death", "dead_biomass"), -yi_smd, yi_smd))
  summary(dt.smd)

fwrite(dt.smd, "effect_sizes/data_with_effect_sizes.csv")



# 5. extra code ---------------------------------------------------------------
# alternatively, for imputing missing SD values 
# we could estimate missing SDs as 1/10 mean
       ## when studies report only a mean but no error measurements it is possible
       ## to assign 1/10 of the mean as SD (Luo et al 2006, Ecology)
       ## https://doi.org/10.1890/04-1724 

no.and.mean[, error_type := "estimated"]
  no.and.mean[, `:=` (Mean_High_Diversity = value_higher_diversity,
                      Mean_Low_Diversity = value_lower_diversity)]
  no.and.mean[, `:=` (SD_High_Diversity = abs(Mean_High_Diversity/10),
                      SD_Low_Diversity = abs(Mean_Low_Diversity/10))]
summary(no.and.mean)

# do the same for studies that are missing values for variability in one of the treatments
# and record also their error_type as estimated
simple <- simple %>%  mutate(error_type = case_when(is.na(SD_High_Diversity) ~ "estimated",
                                                    is.na(SD_Low_Diversity) ~ "estimated",
                                                    TRUE ~ error_type),
                             SD_High_Diversity = ifelse(error_type == "estimated", abs(Mean_High_Diversity/10),
                                                        SD_High_Diversity),
                             SD_Low_Diversity = ifelse(error_type == "estimated", abs(Mean_Low_Diversity/10),
                                                        SD_Low_Diversity))
nrow(simple[is.na(simple$SD_High_Diversity)]) 
nrow(simple[is.na(simple$SD_Low_Diversity)]) 
# no NAs now! :)




