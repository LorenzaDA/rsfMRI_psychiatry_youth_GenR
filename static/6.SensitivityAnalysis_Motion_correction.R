##############
# 6. sensitivity analysis  with more stringent motion correction - static
##############

####
# Data prep
####

rm(list=ls())

# load env 

source("/RELEVANT/PATH/0.Source_file.R")

tab <- "/RELEVANT/PATH/"

# load data 

dd <- readRDS(paste0(indata, "genr_data_sensitivity_rsfMRI_psych.rds"))

summary(dd)


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


#####
# Model specifications 
#####

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

# fit all models

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







write.csv(table_fits, paste0(tab, "model_fits_sensitivity.csv"))



########
# Get estimates 
########

# get standardised estimates 

output_std <- lapply(fits_all, function(x) standardizedSolution(x)) # get std estimates, ci and p vals info 

output_std2 <- as.data.frame(output_std) # change format for saving 

output_std3 <- output_std2[2:3, ] # get the lagged coefficients
output_std3


#### create output table for in-text ####

output_std3$FDR <- as.numeric("NA")

names(output_std3)

output_table <- data.frame(model = 1:12, 
                           lhs = "NA",
                           op = "NA", 
                           rhs = "NA", 
                           est.std = as.numeric("NA"),
                           se = as.numeric("NA"), 
                           z = as.numeric("NA"), 
                           pvalue = as.numeric("NA"), 
                           ci.lower = as.numeric("NA"), 
                           ci.upper = as.numeric("NA"),
                           fdr = as.numeric("NA"))



cols <- paste0("output_std3$", names(output_std3))

cols_p <- grep(".pvalue", cols, value = T)

cols_p <- paste(cols_p, collapse = ", ")

output_table$pvalue <- c(output_std3$cpl_int.pvalue, output_std3$cc_int.pvalue, output_std3$mod_int.pvalue, output_std3$cpl_ext.pvalue, output_std3$cc_ext.pvalue, output_std3$mod_ext.pvalue)
pvals_adj <- p.adjust(output_table$pvalue, method = "fdr") 

output_table$fdr <- pvals_adj



cols_lhs <- grep(".lhs", cols, value = T)

cols_lhs <- paste(cols_lhs, collapse = ", ")


output_table$lhs <- c(output_std3$cpl_int.lhs, output_std3$cc_int.lhs, output_std3$mod_int.lhs, output_std3$cpl_ext.lhs, output_std3$cc_ext.lhs, output_std3$mod_ext.lhs)

cols_rhs <- grep(".rhs", cols, value = T)

cols_rhs <- paste(cols_rhs, collapse = ", ")

output_table$rhs <- c(output_std3$cpl_int.rhs, output_std3$cc_int.rhs, output_std3$mod_int.rhs, output_std3$cpl_ext.rhs, output_std3$cc_ext.rhs, output_std3$mod_ext.rhs)
output_table$op <- "~"


cols_est.std <- grep(".est.std", cols, value = T)

cols_est.std <- paste(cols_est.std, collapse = ", ")



output_table$est.std <- c(output_std3$cpl_int.est.std, output_std3$cc_int.est.std, output_std3$mod_int.est.std, output_std3$cpl_ext.est.std, output_std3$cc_ext.est.std, output_std3$mod_ext.est.std)

cols_se <- grep(".se", cols, value = T)

cols_se <- paste(cols_se, collapse = ", ")

output_table$se <-  c(output_std3$cpl_int.se, output_std3$cc_int.se, output_std3$mod_int.se, output_std3$cpl_ext.se, output_std3$cc_ext.se, output_std3$mod_ext.se)

cols_z <- grep(".z", cols, value = T)

cols_z <- paste(cols_z, collapse = ", ")

output_table$z <-  c(output_std3$cpl_int.z, output_std3$cc_int.z, output_std3$mod_int.z, output_std3$cpl_ext.z, output_std3$cc_ext.z, output_std3$mod_ext.z)



cols_ci.lower <- grep(".ci.lower", cols, value = T)

cols_ci.lower <- paste(cols_ci.lower, collapse = ", ")


output_table$ci.lower <- c(output_std3$cpl_int.ci.lower, output_std3$cc_int.ci.lower, output_std3$mod_int.ci.lower, output_std3$cpl_ext.ci.lower, output_std3$cc_ext.ci.lower, output_std3$mod_ext.ci.lower)

cols_ci.upper <- gsub(".lower", ".upper", cols_ci.lower)

output_table$ci.upper <- c(output_std3$cpl_int.ci.upper, output_std3$cc_int.ci.upper, output_std3$mod_int.ci.upper, output_std3$cpl_ext.ci.upper, output_std3$cc_ext.ci.upper, output_std3$mod_ext.ci.upper)


output_table$model <- paste0(output_table$lhs, " ", output_table$op, " ", output_table$rhs)

out <- output_table[ , c(1, 5, 6, 9, 10, 8, 11)] 
out

rownames(out) <- out$model

out$model <- NULL

out2 <- round(out, digits = 3)

out2


out2$model2 <- gsub("cc", "Modularity", # NB this is because the data was inverted 
                            gsub("mod", "Clustering coefficient", # as above
                                 gsub("cpl", "Characteristic path length", 
                                           gsub("int", "Internalizing problems",
                                                gsub("ext", "Externalizing problems", 
                                                     gsub("_t2", " T2", 
                                                          gsub("_t1", " T1",                                                                                                                      
                                                              rownames(out2))))))))



# separate into temporal direction (brain --> beh; beh --> brain)

BrtoBeh <- grep("cpl_t1|mod_t1|ge_t1|cc_t1", rownames(out2), value = T)

out_BrtoBeh <- out2[rownames(out2) %in% BrtoBeh, ]
out_BehtoBr <- out2[rownames(out2) %notin% BrtoBeh, ]


# model as rownames

rownames(out_BrtoBeh) <- out_BrtoBeh$model2
rownames(out_BehtoBr) <- out_BehtoBr$model2

out_BrtoBeh$model2 <- NULL
out_BehtoBr$model2 <- NULL

out <- rbind(out_BrtoBeh, out_BehtoBr)
write.csv(out, paste0(tab, "results_sensitivity.csv"))
