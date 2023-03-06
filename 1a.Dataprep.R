###################
# 1. Dataprep GenR 
###################
# project: Exploring the longitudinal associations of functional network 
# connectivity and psychiatric symptom changes in youth
# cohort: GenR (focus at @9 and @13)
# author script: Lorenza Dall'Aglio (l.dallaglio@erasmusmc.nl) 
# license CC 4.0 international
# NB before loading the datasets, ensure missingness is dealt with in SPSS
# (all missings and invalid values to missing systems)


####
# environment prep
####

rm(list=ls())

source("PATHTORELEVANTDIRECTORYWHERESOURCEFILEIS/0.Source_file.R")


#####
# load & merge the data - general data & behavioral data
#####
# NB the following code won't work for long format df!

### load and merge the spss data ###

data <- list.files(pattern = ".sav") 
# loads every file in the working directory that ends with the .sav pattern

data 

# 5 elements in the list: general data, CBCL 9, core MRI, puberty, CBCL 13
# the general data should be first for merging purposes


# read in all files into a list

data2 <- list()

for (i in data){
  data2[[i]] <- read.spss(i, to.data.frame = T)
}


# names to lower case 

data3 <- lapply(data2, function(x) {names(x) <- tolower(names(x)); x}) 
data4 <- lapply(data3, as.data.table)


# order df by idc (id var)

lapply(data4, function (x) {x[order(x$idc), ]})


# merge

dd <- Reduce(function(d1, d2) merge(d1, d2, by = intersect(names(d1), names(d2)), all.x = TRUE), data4)  
#reduce = member of apply family enabling to take into account result from previous iteration for the next

dd2 <- as.data.frame(dd)



#####
# load & merge the data - MRI data 
#####

### core MRI data ### 

mri_core <- readRDS("genr_mri_core_data_20220311.rds") 

mri_core <- mri_core[order(mri_core$idc), ] # order by idc 

# this data is already in wide format and can be merged with the rest directly 


### static functional connectivity ###

static <- read.csv("rsfmri_graphtheory_static.csv") # rsfMRI data - graph theory metric

# this df is in long format so it needs to be restructured into wide to merge 
# this needs some preparation first 

static <- static[order(static$idc), ]

static <- as_tibble(static)

# change the str of the var in the df

static$X <- NULL
static$ses <- factor(static$ses)
static$idc <- as.numeric(as.character(static$idc))


## pivot to wider to go from long to wide format 
# first turn the var ses into t - so it can be t1 and t2 for the timepoints 

static$ses <- gsub("ses-" , "", static$ses)
static$ses <- gsub("F09", "t1", static$ses)
static$ses <- gsub("F13", "t2", static$ses)

static$t <- factor(static$ses)
static$acq <- NULL


# now pivot (from long to wide format)

static2 <- static %>% 
  pivot_wider(
    names_from = t, # the variable you want to use for the name of the vars
    id_cols = idc, # variable identifying the ids  
    values_from = c(cpl, ge, cc, mod), # vars that we need to move to wide format 
    names_sep = "_" # separator for the name of the vars - it's gonna be e.g. cpl_t1 
  )



### dynamic functional connectivity ###

dynamic <- load("dfnc_metrics_subjects_clean.Rdata") # dynamic rsfMRI data - df within it loaded as df_id2

dynamic <- df_id2 # rename


## from long to wide df for rsfMRI data

dynamic <- dynamic[order(dynamic$idc), ]
dynamic <- as_tibble(dynamic)
dynamic$t <- recode(dynamic$ses, F09 = "t1", F13 = "t2")  

dynamic2 <- dynamic %>% 
  pivot_wider(
    names_from = t,
    id_cols = idc, 
    values_from = c(mdt1, mdt2, mdt3, mdt4, mdt5, fr1, fr2, fr3, fr4, fr5, nt), 
    names_sep = "_" 
  )



#### merge all data together ####

# merge general data with core mri data 

mri_core <- mri_core[order(mri_core$idc), ] 
dd2 <- dd2[order(dd2$idc), ] 

dd3 <- merge(dd2, mri_core, by = "idc", all = T) 


# merge this data with rsfmri data 

dynamic2 <- dynamic2[order(dynamic2$idc), ]
dd3 <- dd3[order(dd3$idc), ]

dd4 <- merge(dd3, dynamic2, by = 'idc', all = T)  

static2 <- static2[order(static2$idc), ]
dd4 <- dd4[order(dd4$idc), ]

dd5 <- merge(dd4, static2, by = "idc", all = T) 




#######
# Select the variables we need
#######

vars <- c("idc", "idm", "mother", # ids
          "sum_ext_9m", "sum_int_9m",  "sum_int_14", "sum_ext_14", # main vars
          "mdt1_t1", "mdt1_t2", "mdt2_t1", "mdt2_t2", "mdt3_t1", "mdt3_t2", "mdt4_t1", "mdt4_t2", "mdt5_t1", "mdt5_t2",
          "nt_t1", "nt_t2",  # dynamic data
          "cc_t1", "cc_t2", "cpl_t1", "cpl_t2", "mod_t1", "mod_t2", # static data
          "ethninfv3", "gender", "educm5", "agechild_cbcl9m", "agechild_gr1093", "age_child_mri_f09", "age_child_mri_f13", 
          "puberty13", # covs
          "has_braces_mri_f09", "has_braces_mri_f13", "mri_consent_f09", "mri_consent_f13", "rsfmri_has_nii_f09", "rsfmri_has_nii_f13",
          "num_vols_bold_f09", "num_vols_bold_f13", "mean_bold_rms_f09", "mean_bold_rms_f13",
          "exclude_bold_f09", "exclude_bold_f13", "exclude_reg_prob_bold_f09", "exclude_reg_prob_bold_f13", 
          "exclude_incidental_f09", "exclude_incidental_f13",
          "mr750_softwareversionshort_dicom", "startfase3_9", # inc/exc vars
          "sum_anx_9m", "sum_wit_9m", "sum_som_9m", "sum_sop_9m","sum_tho_9m", "sum_att_9m", "sum_rul_9m", "sum_agg_9m", # syndrome scales
          "sum_anx_14", "sum_wit_14","sum_som_14", "sum_sop_14", "sum_tho_14", "sum_att_14", "sum_rul_14", "sum_agg_14",
          "internalizing_problems_tscore_9m", "externalizing_problems_tscore_9m", "internalizing_problems_tscore_14", "externalizing_problems_tscore_14" # the t scores for peer review sens. analyses
) 



# filter the df for the variables we need (specified in vars)

dd6 <- dd5[ , names(dd5) %in% vars] 


####
# Variable recoding & naming
####

dd7 <- dd6 %>% 		mutate(educm_2 = recode(educm5, "no education finished" = 'low', 
                                         "primary" = "low",  "secondary, phase 1" = "int", 
                                         "secondary, phase 2" = "int", "higher, phase 1"= "high",
                                         "higher, phase 2" = "high"),
                        idc = as.character(idc),
                        idm = as.character(idm), 
                        mother = as.character(mother),
                        ethninfv3 = recode(ethninfv3, "Dutch" = "dutch", "Moroccan" = "other", 
                                           "Indonesian" = "other", "Cape Verdian" = "other", "Dutch Antilles" = "other", 
                                           "Surinamese" = "other", "Turkish" = "european", "Surinamese-Creole" = "other",
                                           "Surinamese-Hindustani" = "european", "Surinamese-Unspecified" = "other", 
                                           "African" = "other", "American,western" = "other", "American, non western" = "other", "Asian, western" = "other",
                                           "Asian, non western" = "other", "European" = "european", "Oceanie" = "other")
) %>% 
  
  # rename vars 
  rename(mat_edu = educm_2, ethn = ethninfv3, sex = gender, agecbcl_t1 = agechild_cbcl9m, 
         int_t1 = sum_int_9m, int_t2 = sum_int_14, ext_t1 = sum_ext_9m, ext_t2 = sum_ext_14, 
         agecbcl_t2 = agechild_gr1093, agemri_t1 = age_child_mri_f09, agemri_t2 = age_child_mri_f13,
         puberty = puberty13) 



#####
# Other data prep. 
#####

# create age diff variables for behavior-rsfMRI
# this is necessary because the behavioral and MRI data were collected at slightly different time-points
# i.e., within the same assessment wave, children visited the center multiple times

dd7$agediff_t1 <- dd7$agemri_t1 - dd7$agecbcl_t1 
dd7$agediff_t2 <- dd7$agemri_t2 - dd7$agecbcl_t2

summary(dd7)

dd8 <- dd7


#### get comparison sample for non-response analyses ####

summary(dd7$startfase3_9)  # 8,548 ppl participated at this wave 
# (but they might not have participated in the focus centre for MRI measurements)

dd7b <- dd7[dd7$startfase3_9 == "Yes" & !is.na(dd7$startfase3_9), ]

vars <- c("idc", "ethn", "sex", "mat_edu", "puberty") 

dd7c <- dd7b[ , names(dd7b) %in% vars]

saveRDS(dd7c, "GenR_attrition_analyses_sample.rds")


#### rename/recode vars  ###

dd9 <- dd8 %>% 
  mutate(mat_edu = ordered(mat_edu)) %>%
  rename(id = idc,
         anxdep_t1 = sum_anx_9m,
         anxdep_t2 = sum_anx_14, 
         withdep_t1 = sum_wit_9m,
         withdep_t2 = sum_wit_14, 
         somatic_t1 = sum_som_9m, 
         somatic_t2 = sum_som_14, 
         social_t1 = sum_sop_9m, 
         social_t2 = sum_sop_14, 
         thought_t1 = sum_tho_9m, 
         thought_t2 = sum_tho_14, 
         att_t1 = sum_att_9m, 
         att_t2 = sum_att_14, 
         rulebr_t1 = sum_rul_9m, 
         rulebr_t2 = sum_rul_14, 
         agg_t1 = sum_agg_9m,
         agg_t2 = sum_agg_14)



### transform non-normally distributed behavioral variables ####

# square root transformation of psychiatric problems

cols <- c("int_t1", "int_t2", 
          "ext_t1", "ext_t2",
          "anxdep_t1", "anxdep_t2", 
          "withdep_t1", "withdep_t2", 
          "somatic_t1", "somatic_t2",
          "social_t1", "social_t2", 
          "thought_t1", "thought_t2",
          "att_t1", "att_t2", 
          "rulebr_t1", "rulebr_t2", 
          "agg_t1", "agg_t2")

dd9[cols] <- lapply(dd9[cols], sqrt)



######
# Data inclusion & exclusion (main) 
######
# inclusion: data on beh. and brain at either @9 or @13
# NB for the rsfMRI data, both static and dynamic are obtained in the same way, so either is fine for inclusion
# exclusion: failed QC, substantial incidental findings
# for remaining twins / siblings, only one randomly included


### Apply exclusion criteria to those at @13 with data on MRI ###

# which kids pass the QC at @13?

df <- dd9


df$data_13 <-  ifelse((!is.na(df$mri_consent_f13) & df$mri_consent_f13 == "yes" & 
                         !is.na(df$num_vols_bold_f13) & df$num_vols_bold_f13 == 200 & 
                         !is.na(df$mean_bold_rms_f13) & df$mean_bold_rms_f13 <= 0.25 &
                         is.na(df$exclude_reg_prob_bold_f13) & 
                         !is.na(df$exclude_bold_f13) & df$exclude_bold_f13 == "include" &
                         !is.na(df$exclude_incidental_f13) & df$exclude_incidental_f13 == "include"),
                         "pass", "fail") # pass if they meet all conditions, fail if they don't 
                          # these are all QC conditions and IF 


df$data_13 <- factor(df$data_13)
summary(df$data_13)
# 2212 children pass the QC at 13


# and which kids pass at @9?

df$data_9 <- ifelse((!is.na(df$mri_consent_f09) & df$mri_consent_f09 == "yes" & 
                    !is.na(df$num_vols_bold_f09) & df$num_vols_bold_f09 == 200 &
                    !is.na(df$mean_bold_rms_f09) & df$mean_bold_rms_f09 <= 0.25 &
                    is.na(df$exclude_reg_prob_bold_f09) & 
                    !is.na(df$exclude_bold_f09) & df$exclude_bold_f09 == "include" &
                    !is.na(df$exclude_incidental_f09) & df$exclude_incidental_f09 == "include"), 
                    "pass", "fail") 

df$data_9 <- factor(df$data_9)
summary(df$data_9)
# 2592 children pass


## select children with passing data at @9 or @13
# NB before set to NA the MRI data where the QC wasn't passed
# so you're sure that children with useable data at one time point
# but not at the other are included correctly (i.e. only their useable data is included)

df2 <- df

for(i in 1:nrow(df2)){
  if(!is.na(df2$data_9[i]) & df2$data_9[i] == "fail"){ # if the child failed the QC at 9 or has no data
    df2$mod_t1[i] <- NA # set all values for rsfMRI at that time point to NA
    df2$cc_t1[i] <- NA
    df2$cpl_t1[i] <- NA
    df2$nt_t1[i] <- NA
    df2$mdt1_t1[i] <- NA
    df2$mdt2_t1[i] <- NA
    df2$mdt3_t1[i] <- NA
    df2$mdt4_t1[i] <- NA
    df2$mdt5_t1[i] <- NA
  }
}




for(i in 1:nrow(df2)){
  if(!is.na(df2$data_13[i]) & df2$data_13[i] == "fail"){ # same as above but for @13
    df2$mod_t2[i] <- NA
    df2$cc_t2[i] <- NA
    df2$cpl_t2[i] <- NA
    df2$nt_t2[i] <- NA
    df2$mdt1_t2[i] <- NA
    df2$mdt2_t2[i] <- NA
    df2$mdt3_t2[i] <- NA
    df2$mdt4_t2[i] <- NA
    df2$mdt5_t2[i] <- NA
  }
}


# subset for children with a pass at 9 or 13

df3 <- subset(df2, (data_9 == "pass" | data_13 == "pass")) # 3767 
# 3767


# subset for children with all rsfMRI data at either 9 or 13
# and with psych problems (inclusion criteria)

df3$has_both_F09 <- NA
df3$has_both_F14 <- NA

# has data at @9

df3$has_both_F09[!is.na(df3$mdt1_t1) & !is.na(df3$mdt2_t1) & !is.na(df3$mdt3_t1) & !is.na(df3$mdt4_t1) & !is.na(df3$mdt5_t1) & !is.na(df3$nt_t1)  
                & !is.na(df3$int_t1) & !is.na(df3$ext_t1)]<- 1





# has data at @13

df3$has_both_F14[!is.na(df3$mdt1_t2) & !is.na(df3$mdt2_t2) & !is.na(df3$mdt3_t2) & !is.na(df3$mdt4_t2) & !is.na(df3$mdt5_t2) & !is.na(df3$nt_t2)
                & !is.na(df3$int_t2) & !is.na(df3$ext_t2)]<- 1


table(df3$has_both_F09) # 2208 with data at 9 that passed the QC
table(df3$has_both_F14) # 1939 with data at 13 that passed the QC


# select for children with rsfMRI and psych data at either time point

df3 <- subset(df3, df3$has_both_F09==1 | df3$has_both_F14==1) # 3296 


# randomly include one sibling/twin for each set

sum(duplicated(df3$mother)) 

set.seed(2022) # set a seed for randomisation

df3 <- df3[sample(nrow(df3)), ] # randomly order data

final_dd <- df3[!duplicated(df3$mother), ] # 3131 children in the final df


# subset for sensitivity analysis where both data at 9 and 13 is needed 

temp <- subset(final_dd, final_dd$has_both_F09==1 & final_dd$has_both_F14==1)  
# 806 ppl

saveRDS(temp, paste0(indata, "completecases_GenR_rsfmri_psych.rds"))


######
# Clean up the dataframe
######

names(final_dd)

exc <- c("idm", "agecbcl_t1", "agecbcl_t2", "mother", "startfase3_9", "educm5", 
         "mr750_softwareversionshort_dicom", "mri_consent_f09", "mri_consent_f13", 
         "rsfmri_has_nii_f09", "rsfmri_has_nii_f13", "has_braces_mri_f09", 
         "has_braces_mri_f13", "exclude_incidental_f09", "exclude_incidental_f13", 
         "num_vols_bold_f09", "num_vols_bold_f13", "mean_bold_rms_f09", "mean_bold_rms_f13",
         "exclude_bold_f09", "exclude_bold_f13", "exclude_reg_prob_bold_f09", 
         "exclude_reg_prob_bold_f13", "data_13",
         "has_both_F09", "has_both_F14", "data_9") # set the variables which we do not need

final_dd <- final_dd[ , names(final_dd) %notin% exc]


# Checks:
# are there enough people in each category for categorical variable? 
# are the min and max values reasonable for the continuous variables? 
# the MDT can be very high for some ppl (almost 200, but that is actually normal as based
# on Lopez-Vincente et al. 2021, Frontiers in Systems Neuroscience)
# all fine


######
# save 
######

saveRDS(final_dd, paste0(indata, "genr_data_main_rsfMRI_psych.rds"))

# 3131 children




