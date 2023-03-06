##################
# 0b. Preparation connectivity data - graph theory files
##################
# working in the server 
# so from bash load R

###
# Set environment
###

library(stringr)
library(tidyverse)
library(dplyr)

indata <- "/PUT/RELEVANT/PATH/WHERE/DATA/IS/"
setwd(indata)

dd <- read.csv("combined.csv", header = F) # read in the file where rsfmri graph theory data was saved


####
# Format data
####

# give colnames to the df 

names(dd) <- c("inputFile", "subj", "ses", "acq", "conmatType", "cpl",
               "ge", "ge2", "cc", "mod")


# change id format so it becomes like IDC 

dd$subj[1:5] # this is an example of what ids look like rn 
# sub-9_ses-F13_task-rest_acq-0006_run-0_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii 
# you just need the number after sub-

dd$id <- NA


# remove all info till the subject id 

dd <- tibble(dd)

dd2 <- dd %>% mutate(
  id =str_remove(subj, pattern = "^.*?-"),
  # remove till the - symbol 
)


dd2$id2 <- NA

for(i in 1:nrow(dd2)){
  dd2$idc[i] <- (str_split_fixed(dd2$id, "_", 2)[i])
}

# select cols 

dd3 <- dd2[ , c("ses", "acq", "cpl", "cc", "mod", "idc")]


# save

write.csv(dd3, paste0(indata, "rsfmri_graphtheory_static.csv"))
