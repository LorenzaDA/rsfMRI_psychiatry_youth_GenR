##############################################
# 12. Sensitivity analyses for MDT 5 analyses 
###############################################
# issue: some associations with MDT5 seemed to be driven by 
# a few observations on the far end (ext, agg, att)
# winsorize those and see if associations stayed.

rm(list=ls())

# load env 

source("/PUT/RELEVANT/PATHS/0.Source_file.R")

# load data 

dd <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds"))

tab <- "/PUT/RELEVANT/PATH/FOR/OUTPUT/"


####
# factor scaling 
####
# this is needed so that the variables are more or less on the same order of magnitude
# FA/MD values are generally 100 times smaller than behavioural ones. --> use scale factor of 100 

scale_factor = 100

summary(dd$mdt5_t1) # 0 to 170
summary(dd$mdt5_t2) # 0 to 81

# winsorize extreme vals - generally it's 95th quantile
quantile(dd$mdt5_t1, 0.95, na.rm = T) # 43.65 

# winsorize with max val the 95th quantile 
dd$mdt5_t1_try <- DescTools::Winsorize(dd$mdt5_t1, minval = NULL, maxval = 43.65, na.rm = T)


# scale
dd$mdt5_t1 <- dd$mdt5_t1_try/scale_factor
dd$mdt5_t2 <- dd$mdt5_t2/scale_factor

####
# Recoding of vars for lavaan
####
# Before running lavaan, variables need to be coded in a certain way if they are not continuous
# this varies depending on whether the var is exogeneous (ind) or endogeneous (dep)
# more info at https://lavaan.ugent.be/tutorial/cat.html


### Make ordered variables, continuous/numeric ###
# maternal education only here

dd$mat_edu_num <- as.numeric(dd$mat_edu)



### dummy coding for binary/nominal variables ### 
# to make dummy cols from factors or character vars you can use dummy_cols()
# NB check which vars have factors or chr 
# id is chr - change it just temporarily so it is not transformed into a dummy 
# only nominal vars here are ethnicity and sex

rownames(dd) <- dd$id
dd$id <- NULL

dd2 <- fastDummies::dummy_cols(dd)   
# NB the function dichotomises everything w/o setting ref cateogory so you need to take care of that


# delete reference vars 
# and dummies with NAs 

dd2$ethn_dutch <- NULL
dd2$ethn_NA <- NULL
dd2$sex_girl <- NULL


# model for MDT5 with att, agg, ext

m_general  <- '

# lagged and stability 
phen_t2 + mdt5_t2 ~ phen_t1 + mdt5_t1

# covariances
phen_t1 ~~ mdt5_t1
phen_t2 ~~ mdt5_t2

# autocorr
phen_t1 ~~ phen_t1
mdt5_t1 ~~ mdt5_t1
phen_t2 ~~ phen_t2
mdt5_t2 ~~ mdt5_t2

# covs
phen_t2 ~ sex_boy + ethn_other + ethn_european + mat_edu_num

mdt5_t2 ~ puberty + agediff_t2 + agemri_t2

mdt5_t1 ~ agediff_t1 + agemri_t1

'

# create empty lists where to put the model specifications for each measure
# From the general model specifications, we create specific model specifications for each
# syndrome scale. The "var" within the specification will be substituted with each rsfMRI measure

specifications_mdt5 <- list()


# create vector for the substitution of the var specifications into syndrome scale names 

names <- c("ext_", "att_", "agg_")


# loop to create the various model specifications per syndrome scale 

for (i in names){
  specifications_mdt5[i] <- gsub("phen_", i, m_general)
}


fits_all <- sapply(specifications_mdt5, function(x) sem(x, data = dd2, missing = "fiml", fixed.x = F))

table_fits <- sapply(fits_all, function(x) fitMeasures(x, c("cfi", "tli", "rmsea", "srmr")))

table_fits


# get standardised estimates 

output_std <- lapply(fits_all, function(x) standardizedSolution(x)) # get std estimates, ci and p vals info 

output_std2 <- as.data.frame(output_std) # change format for saving 

output_std3 <- output_std2[2:3, ] # get the lagged coefficients

cols <- grep(".pvalue", names(output_std3), value = T)

ps <- output_std3[1, names(output_std3) %in% cols]

write.csv(output_std3, paste0(tab, "mdt5_winsorized_sensanalysis.csv"))
