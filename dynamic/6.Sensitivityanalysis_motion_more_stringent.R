#############
# sensitivity analysis  with more stringent motion correction - dynamic
#############

####
# data prep
####

rm(list=ls())

# load env 

source("/PUT/RELEVANT/PATH/0.Source_file.R")

# load data 

dd <- readRDS(paste0(indata, "genr_data_sensitivity_rsfMRI_psych.rds"))

tab <- "/RELEVANT/OUTPUT/PATH"


####
# factor scaling 
####
# this is needed so that the variables are more or less on the same order of magnitude

scale_factor = 100


mdts <- grep("mdt", names(dd), value = T)

dd[mdts] <- dd[mdts] / scale_factor



####
# Recoding of vars for lavaan
####

dd$mat_edu_num <- as.numeric(dd$mat_edu)

dd$income <- as.numeric(dd$income)

rownames(dd) <- dd$id
dd$id <- NULL

dd2 <- fastDummies::dummy_cols(dd)   

dd2$ethn_dutch <- NULL
dd2$ethn_NA <- NULL
dd2$sex_girl <- NULL


####
# Model specifications
####


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

names <- c("mdt1_", "mdt2_", "mdt3_", "mdt4_", "mdt5_", "nt_")


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
# model evaluations
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

write.csv(table_fits, paste0(tab, "model_fits_sensitivity_rounded.csv"))


#####
# model estimates
#####

# get standardised estimates 

output_std <- lapply(fits_all, function(x) standardizedSolution(x)) # get std estimates, ci and p vals info 

output_std2 <- as.data.frame(output_std) # change format for saving 

output_std3 <- output_std2[2:3, ] # get the lagged coefficients


#### create output table for in-text ####

output_std3$FDR <- as.numeric("NA")

names(output_std3)

output_table <- data.frame(model = 1:24, 
                           lhs = "NA",
                           op = "NA", 
                           rhs = "NA", 
                           est.std = as.numeric("NA"),
                           se = as.numeric("NA"), 
                           pvalue = as.numeric("NA"), 
                           fdr = as.numeric("NA"))



cols <- paste0("output_std3$", names(output_std3))

cols_p <- grep(".pvalue", cols, value = T)

cols_p <- paste(cols_p, collapse = ", ")

output_table$pvalue <- c(output_std3$mdt1_int.pvalue, output_std3$mdt2_int.pvalue, output_std3$mdt3_int.pvalue, output_std3$mdt4_int.pvalue, output_std3$mdt5_int.pvalue, output_std3$nt_int.pvalue, output_std3$mdt1_ext.pvalue, output_std3$mdt2_ext.pvalue, output_std3$mdt3_ext.pvalue, output_std3$mdt4_ext.pvalue, output_std3$mdt5_ext.pvalue, output_std3$nt_ext.pvalue)

pvals_adj <- p.adjust(output_table$pvalue, method = "fdr") 

output_table$fdr <- pvals_adj



cols_lhs <- grep(".lhs", cols, value = T)

cols_lhs <- paste(cols_lhs, collapse = ", ")


output_table$lhs <- c(output_std3$mdt1_int.lhs, output_std3$mdt2_int.lhs, output_std3$mdt3_int.lhs, output_std3$mdt4_int.lhs, output_std3$mdt5_int.lhs, output_std3$nt_int.lhs, output_std3$mdt1_ext.lhs, output_std3$mdt2_ext.lhs, output_std3$mdt3_ext.lhs, output_std3$mdt4_ext.lhs, output_std3$mdt5_ext.lhs, output_std3$nt_ext.lhs)



cols_rhs <- grep(".rhs", cols, value = T)

cols_rhs <- paste(cols_rhs, collapse = ", ")

output_table$rhs <- c(output_std3$mdt1_int.rhs, output_std3$mdt2_int.rhs, output_std3$mdt3_int.rhs, output_std3$mdt4_int.rhs, output_std3$mdt5_int.rhs, output_std3$nt_int.rhs, output_std3$mdt1_ext.rhs, output_std3$mdt2_ext.rhs, output_std3$mdt3_ext.rhs, output_std3$mdt4_ext.rhs, output_std3$mdt5_ext.rhs, output_std3$nt_ext.rhs)

output_table$op <- "~"


cols_est.std <- grep(".est.std", cols, value = T)

cols_est.std <- paste(cols_est.std, collapse = ", ")



output_table$est.std <- c(output_std3$mdt1_int.est.std, output_std3$mdt2_int.est.std, output_std3$mdt3_int.est.std, output_std3$mdt4_int.est.std, output_std3$mdt5_int.est.std, output_std3$nt_int.est.std, output_std3$mdt1_ext.est.std, output_std3$mdt2_ext.est.std, output_std3$mdt3_ext.est.std, output_std3$mdt4_ext.est.std, output_std3$mdt5_ext.est.std, output_std3$nt_ext.est.std)


cols_se <- grep(".se", cols, value = T)

cols_se <- paste(cols_se, collapse = ", ")

output_table$se <-  c(output_std3$mdt1_int.se, output_std3$mdt2_int.se, output_std3$mdt3_int.se, output_std3$mdt4_int.se, output_std3$mdt5_int.se, output_std3$nt_int.se, output_std3$mdt1_ext.se, output_std3$mdt2_ext.se, output_std3$mdt3_ext.se, output_std3$mdt4_ext.se, output_std3$mdt5_ext.se, output_std3$nt_ext.se)


output_table$model <- paste0(output_table$lhs, " ", output_table$op, " ", output_table$rhs)

out <- output_table[ , c(1, 5, 6, 7, 8)] 
out

rownames(out) <- out$model
out$model <- NULL

out2 <- round(out, digits = 3)

# separate into temporal direction (brain --> beh; beh --> brain)

BrtoBeh <- grep("mdt1_t1|mdt2_t1|mdt3_t1|mdt4_t1|mdt5_t1| nt_t1", rownames(out2), value = T)

out_BrtoBeh <- out2[rownames(out2) %in% BrtoBeh, ]
out_BehtoBr <- out2[rownames(out2) %notin% BrtoBeh, ]


out2 <- rbind(out_BrtoBeh, out_BehtoBr)

rownames(out2) <- gsub("mdt1", "MDT state 1", 
                   gsub("mdt2", "MDT state 2",
                        gsub("mdt3", "MDT state 3", 
                             gsub("mdt4", "MDT state 4",
                                  gsub("mdt5", "MDT state 5",
                                       gsub("nt ", "NT ",
                                            gsub("int", "Internalizing problems",
                                                 gsub("ext", "Externalizing problems", 
                                                      gsub("_t1", " T1",
                                                           gsub("_t2", " T2",
                                                                rownames(out2)))))))))))


write.csv(out2, paste0(tab, "results_sensitivity_rounded.csv"))


