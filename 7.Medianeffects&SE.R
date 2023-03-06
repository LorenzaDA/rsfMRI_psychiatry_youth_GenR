##############
# Check median effect sizes and std errors for everything
##############


indata_stat <- "/PATH/WHERE/STATIC/DATA/WAS/SAVED/"
indata_dyn <- "/PATH/WHERE/DYNAMIC/DATA/WAS/SAVED/"

stat2 <-  read.csv(paste0(indata_stat, "results_syndromes.csv"))
stat1 <- read.csv(paste0(indata_stat, "results_main.csv"))

dyn1 <- read.csv(paste0(indata_dyn, "results_main.csv"))
dyn2 <- read.csv(paste0(indata_dyn, "results_syndromes.csv"))

stat <- rbind(stat2, stat1)
dyn <- rbind(dyn1,dyn2)

median(abs(stat$est.std)) # 0.014
median(abs(dyn$est.std)) # 0.021

median(abs(stat$se)) # 0.023
median(abs(dyn$se)) # 0.022


######
# check n children with clinically relevant symptoms
######

source("/PATH/WHERE/SOURCE/FILE/IS/0.Source_file.R")

rs <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds"))

quantile(rs$int_t1, 0.93, na.rm = T)
# 3.61

quantile(rs$int_t2, 0.93, na.rm = T) 
# 4 

quantile(rs$ext_t1, 0.93, na.rm = T)
# 3.32

quantile(rs$ext_t2, 0.93, na.rm = T)
# 3.61

temp <- rs[!is.na(rs$int_t1) & (rs$int_t1 > 3.61 | rs$int_t1 == 3.61), ]
# 188

temp <- rs[!is.na(rs$int_t2) & (rs$int_t2 > 4 | rs$int_t2 == 4), ]
# 210

temp <- rs[!is.na(rs$ext_t1) & (rs$ext_t1 > 3.32 | rs$ext_t1 == 3.32), ]
# 198

temp <- rs[!is.na(rs$ext_t2) & (rs$ext_t2 > 3.61 | rs$ext_t2 == 3.61), ]
# 158

#  188.5 = average n kids with clinically relevant problems at either time point

