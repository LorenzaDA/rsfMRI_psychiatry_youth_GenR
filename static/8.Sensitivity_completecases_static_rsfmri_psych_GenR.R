#####################
# 9. CLPMs - with complete cases
#####################
# sensitivity analysis
# asked in peer review

####
# Data prep
####

rm(list=ls())

# load env 

source("/RELEVANT/PATH/0.Source_file.R")

tab <- "/RELEVANT/PATH/"

# load data 

dd <- readRDS(paste0(indata, "completecases_GenR_rsfmri_psych.rds"))


####
# factor scaling 
####
# this is needed so that the variables are more or less on the same order of magnitude

scale_factor = 10

rsfmri <- c("cpl_t1", "cpl_t2", "cc_t1", "cc_t2", "mod_t1", "mod_t2")

dd[rsfmri] <- dd[rsfmri] * scale_factor


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
# to make dummy cols from factors or character vars you can use dummy_cols(), which will 
# automatically turn factors and characters into dummy vars
# NB id is chr - change it just temporarily so it is not transformed into a dummy 
# the only nominal vars here are ethnicity and sex

rownames(dd) <- dd$id
dd$id <- NULL

dd2 <- fastDummies::dummy_cols(dd)   
# NB the function dichotomises everything w/o setting ref cateogory so you need to take care of that


# delete reference vars 
# and dummies with NAs 

dd2$ethn_dutch <- NULL
dd2$ethn_NA <- NULL
dd2$sex_girl <- NULL


#####
# Model specifications 
#####
# specify the models you want to test
# of note, these specifications are kept general in part (var_t1 & var_t1)
# this is so that different vars can be tested with the same models (i.e., CPL, CC, etc..)

m_int  <- '

# lagged and stability 
int_t2 + var_t2 ~ int_t1 + var_t1

# covariances
int_t1 ~~ var_t1
int_t2 ~~ var_t2

# autocorr
int_t1 ~~ int_t1
var_t1 ~~ var_t1
int_t2 ~~ int_t2
var_t2 ~~ var_t2

# covs
int_t2 ~ sex_boy + ethn_other + ethn_european + mat_edu_num

var_t2 ~ puberty + agediff_t2 + agemri_t2

var_t1 ~ agediff_t1 + agemri_t1

'



m_ext <- '
# lagged and stability 
ext_t2 + var_t2 ~ ext_t1 + var_t1

# covariances
ext_t1 ~~ var_t1
ext_t2 ~~ var_t2

# autocorr
ext_t1 ~~ ext_t1
var_t1 ~~ var_t1
ext_t2 ~~ ext_t2
var_t2 ~~ var_t2

# covs
ext_t2 ~ sex_boy + ethn_other + ethn_european + mat_edu_num

var_t2 ~ puberty + agediff_t2 + agemri_t2

var_t1 ~ agediff_t1 + agemri_t1

'


# create empty lists where to put the model specifications for each measure
# From the general model specifications, we create specific model specifications for each
# syndrome scale. The "var" within the specification will be substituted with each rsfMRI measure

specifications_int <- list()
specifications_ext <- list()

# create vector for the substitution of the var specifications into syndrome scale names 

names <- c("cpl_", "cc_", "mod_")


# loop to create the various model specifications per syndrome scale 

for (i in names){
  specifications_int[i] <- gsub("var_", i, m_int) # gsub changes the "var_" string into every element of names in the m_fa specification
  specifications_ext[i] <- gsub("var_", i, m_ext) # as above, but for md specifications
}


names(specifications_int) <- paste0(names(specifications_int), "int")
names(specifications_ext) <- paste0(names(specifications_ext), "ext")


# put all model specifications together

specifications <- append(specifications_int, specifications_ext)


#####
# Model evaluation
#####

# fit all models with sapply

fits_all <- sapply(specifications, function(x) sem(x, data = dd2, missing = "fiml", fixed.x = F))


# get all model fit measures that you are interested in 

table_fits <- sapply(fits_all, function(x) fitMeasures(x, c("cfi", "tli", "rmsea", "srmr")))

table_fits <- as.data.frame(t(table_fits)) # for saving, transform into dataframe format after transposing (so fit indices are the cols and models the rows)

table_fits$model <- rownames(table_fits) # create new col with the model name 

table_fits <- table_fits[, c(5, 3, 4, 1, 2)] # change order of cols

table_fits$model <- NULL # for rounding we need to have all vars as numeric

table_fits <- round(table_fits, digits = 3) # round to 3 decimals

rownames(table_fits) <- gsub("cc", "Modularity", # NB this is because the data was inverted 
                             gsub("mod", "Clustering coefficient", # as above
                                  gsub("cpl", "Characteristic path length",
                                       gsub("int", "Internalizing problems",
                                            gsub("ext", "Externalizing problems", 
                                                 gsub("_", " - ", 
                                                      rownames(table_fits)))))))




########
# Get estimates 
########

# get standardised estimates for all tested models

output_std <- lapply(fits_all, function(x) standardizedSolution(x)) 

output_std2 <- as.data.frame(output_std) # change format for saving - from list to df


# select only lagged paths (i.e. the longitudinal assoc. between brain and behavior)

output_std3 <- output_std2[2:3, ] 

#### create output table for in-text ####

output_std3$FDR <- as.numeric("NA")

# create an empty df where the output will go 

output_table <- data.frame(model = 1:12, 
                           lhs = "NA",
                           op = "NA", 
                           rhs = "NA", 
                           est.std = as.numeric("NA"),
                           se = as.numeric("NA"), 
                           pvalue = as.numeric("NA"), 
                           fdr = as.numeric("NA"))


# get the cols that will be needed to fill in the output df

cols <- paste0("output_std3$", names(output_std3))


## get p values

cols_p <- grep(".pvalue", cols, value = T) # select only those for the p value

cols_p <- paste(cols_p, collapse = ", ") # put together the elements of cols_p, and separate them by a comma
# copy paste and put in the vector below

output_table$pvalue <- c(output_std3$cpl_int.pvalue, 
                         output_std3$cc_int.pvalue, 
                         output_std3$mod_int.pvalue, 
                         output_std3$cpl_ext.pvalue,
                         output_std3$cc_ext.pvalue, 
                         output_std3$mod_ext.pvalue)
# put p values into the output table 


## get FDR p

pvals_adj <- p.adjust(output_table$pvalue, method = "fdr") # FDR correction to p vals

output_table$fdr <- pvals_adj # fill in the table col for fdr with the obtained values


## get the left part of the lagged path regression (i.e. outcome)

cols_lhs <- grep(".lhs", cols, value = T)

cols_lhs <- paste(cols_lhs, collapse = ", ")

output_table$lhs <- c(output_std3$cpl_int.lhs, 
                      output_std3$cc_int.lhs, 
                      output_std3$mod_int.lhs, 
                      output_std3$cpl_ext.lhs, 
                      output_std3$cc_ext.lhs, 
                      output_std3$mod_ext.lhs)


## get the right part of the lagged path regression (i.e. predictor)

cols_rhs <- grep(".rhs", cols, value = T)

cols_rhs <- paste(cols_rhs, collapse = ", ")

output_table$rhs <- c(output_std3$cpl_int.rhs, 
                      output_std3$cc_int.rhs, 
                      output_std3$mod_int.rhs, 
                      output_std3$cpl_ext.rhs, 
                      output_std3$cc_ext.rhs, 
                      output_std3$mod_ext.rhs)

output_table$op <- "~"


## estimates

cols_est.std <- grep(".est.std", cols, value = T)

cols_est.std <- paste(cols_est.std, collapse = ", ")

output_table$est.std <- c(output_std3$cpl_int.est.std, 
                          output_std3$cc_int.est.std, 
                          output_std3$mod_int.est.std, 
                          output_std3$cpl_ext.est.std, 
                          output_std3$cc_ext.est.std, 
                          output_std3$mod_ext.est.std)


## SE

cols_se <- grep(".se", cols, value = T)

cols_se <- paste(cols_se, collapse = ", ")

output_table$se <-  c(output_std3$cpl_int.se,
                      output_std3$cc_int.se, 
                      output_std3$mod_int.se, 
                      output_std3$cpl_ext.se, 
                      output_std3$cc_ext.se,
                      output_std3$mod_ext.se)


## model

output_table$model <- paste0(output_table$lhs, " ", output_table$op, " ", output_table$rhs)


## clean up the table 

# rename 

output_table$model2 <- gsub("cc", "Modularity", # NB this is because the data was inverted 
                            gsub("mod", "Clustering coefficient", # as above
                                 gsub("cpl", "Characteristic path length", 
                                      gsub("int", "Internalizing problems",
                                           gsub("ext", "Externalizing problems", 
                                                gsub("_t2", " T2", 
                                                     gsub("_t1", " T1",                                                                                                                      
                                                          output_table$model)))))))


# separate into temporal direction (brain --> beh; beh --> brain)

BrtoBeh <- grep("cpl_t1|mod_t1|cc_t1", output_table$model, value = T)

out_BrtoBeh <- output_table[output_table$model %in% BrtoBeh, ]
out_BehtoBr <- output_table[output_table$model %notin% BrtoBeh, ]


## clean up 

# model as rownames

rownames(out_BrtoBeh) <- out_BrtoBeh$model2
rownames(out_BehtoBr) <- out_BehtoBr$model2

# select only relevant cols

out_BrtoBeh <- out_BrtoBeh[, c(5, 6, 7, 8)]
out_BehtoBr <- out_BehtoBr[, c(5, 6, 7, 8)]

# round

out_BrtoBeh <- round(out_BrtoBeh, digits = 3)
out_BehtoBr <- round(out_BehtoBr, digits = 3)

## save 

write.csv(out_BrtoBeh, paste0(tab, "results_sens_CC_BrtoBeh.csv"))
write.csv(out_BehtoBr, paste0(tab, "results_sens_CC_BehtoBr.csv"))

