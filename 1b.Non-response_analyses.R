###############################################
# 1b. Non-response analyses GenR rsfMRI static
###############################################

#####
# set environment
#####

rm(list=ls())

# load env 

source("PUTRELEVANTPATHWHERESOURCEFILEIS/0.Source_file.R")

# load data 

tab <- "/PATHWHEREYOUWANTTHETABLES/"

dd <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds")) # sub-sample

comp <- readRDS(paste0(indata, "GenR_attrition_analyses_sample.rds")) # sample for comparison  - 8548 children invited at @9


# check missingness 

missings <- lapply(dd[,c("puberty", "sex", "ethn", "mat_edu")], function(x) length(which(is.na(x))))

missings2 <- lapply(comp[,c("puberty", "sex", "ethn", "mat_edu")], function(x) length(which(is.na(x))))

cbind(missings, missings2) 


comp_clean <- comp[!is.na(comp$ethn) & !is.na(comp$puberty) & !is.na(comp$mat_edu) & !is.na(comp$sex), ] # 3601 ppl
dd_clean <- dd[!is.na(dd$ethn) & !is.na(dd$puberty) & !is.na(dd$mat_edu) & !is.na(dd$sex), c("id", "ethn", "sex", "mat_edu", "puberty")] # 902 ppl 


# select relevant vars and compare

covs <- c("puberty", "sex", "ethn", "mat_edu")


tab1_subsample <- CreateTableOne(c(covs), data = dd_clean) 
tab1_fullsample <- CreateTableOne(c(covs), data = comp_clean)
comparison <- compare_them(dd_clean, comp_clean, covs) 


# save

comparison_table <- print(comparison, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE) # statistical comparison bw the two samples
sumstats_subsample <- print(tab1_subsample, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE) # summary stats table for subsample
sumstats_fullsample <- print(tab1_fullsample, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE) # summary stats table for full sample


## Save to a CSV file

write.csv2(comparison_table, paste0(tab, "Attritionanalyses_proportioncomparison.csv"))
write.csv2(sumstats_subsample, paste0(tab, "Attritionanalyses_proportions_subsample.csv"))
write.csv2(sumstats_fullsample, paste0(tab, "Attritionanalyses_proportions_fullsample.csv"))


