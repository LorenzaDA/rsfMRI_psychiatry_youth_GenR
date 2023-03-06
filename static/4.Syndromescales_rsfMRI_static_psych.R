#####################################
# 4. Syndrome scales - rsfMRI static
#####################################
# project: Exploring the longitudinal associations of functional network connectivity
# and psychiatric symptom changes in youth
# cohort: GenR (focus at @9 and @13)
# author script: Lorenza Dall'Aglio (l.dallaglio@erasmusmc.nl) 
# license CC 4.0 international


#####
# Set environment
#####

source("/PATH/SOURCE/FILE/0.Source_file.R")

# load data 

dd <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds"))

tab <- "/PATH/OUTPUT/"


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
# more info at https://lavaan.ugent.be/tutorial/cat.html

### Make ordered variables, continuous/numeric ###
# maternal education only here

dd$mat_edu_num <- as.numeric(dd$mat_edu)

### dummy coding for binary/nominal variables ### 

rownames(dd) <- dd$id
dd$id <- NULL

dd2 <- fastDummies::dummy_cols(dd)   

# delete reference vars 
# and dummies with NAs 

summary(dd2$ethn) # largest category = dutch 
dd2$ethn_dutch <- NULL
dd2$ethn_NA <- NULL

summary(dd2$sex) # largest cat = girls
summary(dd2$sex)
dd2$sex_girl <- NULL


#####
# Model specifications 
#####

# general model which you will later fill in for each syndrome scale and static measure

m_general  <- '

# lagged and stability 
phen_t2 + state_t2 ~ phen_t1 + state_t1

# covariances
phen_t1 ~~ state_t1
phen_t2 ~~ state_t2

# autocorr
phen_t1 ~~ phen_t1
state_t1 ~~ state_t1
phen_t2 ~~ phen_t2
state_t2 ~~ state_t2

# covs
phen_t2 ~ sex_boy + ethn_other + ethn_european + mat_edu_num

state_t2 ~ puberty + agediff_t2 + agemri_t2

state_t1 ~ agediff_t1 + agemri_t1

'

# substitute the models "state" witht the specific measures

m_cc <- gsub("state", "cc", m_general)
m_mod <- gsub("state", "mod", m_general)
m_cpl <- gsub("state", "cpl", m_general)



### Specific model specifications 

# create empty lists for each rsfmri measure
# within each list you need to have models for each syndrome scale

specifications_cc <- list()
specifications_mod <- list()
specifications_cpl <- list()


# create vector for the substitution of the var specifications into syndrome scale names 

names <- c("anxdep_", "withdep_", "somatic_", "social_", "thought_", "att_", "agg_", "rulebr_")


# loop to create the various model specifications per syndrome scale within a list 
# of rsfmri data

for (i in names){
  specifications_cc[i] <- gsub("phen_", i, m_cc) # gsub changes the "var_" string into every element of names in the m_fa specification
  specifications_mod[i] <- gsub("phen_", i, m_mod)
  specifications_cpl[i] <- gsub("phen_", i, m_cpl)
}


names(specifications_cc) <- paste0(names(specifications_cc), "cc")
names(specifications_mod) <- paste0(names(specifications_mod), "mod")
names(specifications_cpl) <- paste0(names(specifications_cpl), "cpl")


# put all specifications together 

specifications <- c(specifications_cc, specifications_mod, specifications_cpl) 


######
# Model evaluation
######

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
                                            gsub("anxdep", "Anxious/Depressed",
                                                 gsub("withdep", "Withdrawn/Depressed", 
                                                      gsub("somatic", "Somatic complaints",
                                                           gsub("social", "Social problems",
                                                                gsub("thought", "Thought problems",
                                                                     gsub("att", "Attention problems",
                                                                          gsub("agg", "Aggressive behaviors", 
                                                                               gsub("rulebr", "Rule-breaking behaviors",
                                                      gsub("_", " - ", 
                                                           rownames(table_fits)))))))))))))





write.csv(table_fits, paste0(tab, "model_fits_syndromes.csv"))




########
# Get estimates 
########

# get standardised estimates 

output_std <- lapply(fits_all, function(x) standardizedSolution(x)) # get std estimates, ci and p vals info 

output_std2 <- as.data.frame(output_std) # change format for saving 

output_std3 <- output_std2[2:3, ] # get the lagged coefficients


#### create output table for in-text ####

output_std3$FDR <- as.numeric("NA")

names(output_std3)

output_table <- data.frame(model = 1:48, 
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

output_table$pvalue <- c(output_std3$anxdep_cc.pvalue, output_std3$withdep_cc.pvalue, output_std3$somatic_cc.pvalue, output_std3$social_cc.pvalue, output_std3$thought_cc.pvalue, output_std3$att_cc.pvalue, output_std3$agg_cc.pvalue, output_std3$rulebr_cc.pvalue, output_std3$anxdep_mod.pvalue, output_std3$withdep_mod.pvalue, output_std3$somatic_mod.pvalue, output_std3$social_mod.pvalue, output_std3$thought_mod.pvalue, output_std3$att_mod.pvalue, output_std3$agg_mod.pvalue, output_std3$rulebr_mod.pvalue, output_std3$anxdep_cpl.pvalue, output_std3$withdep_cpl.pvalue, output_std3$somatic_cpl.pvalue, output_std3$social_cpl.pvalue, output_std3$thought_cpl.pvalue, output_std3$att_cpl.pvalue, output_std3$agg_cpl.pvalue, output_std3$rulebr_cpl.pvalue)
pvals_adj <- p.adjust(output_table$pvalue, method = "fdr") 

output_table$fdr <- pvals_adj



cols_lhs <- gsub(".pvalue", ".lhs", cols_p)


output_table$lhs <- c(output_std3$anxdep_cc.lhs, output_std3$withdep_cc.lhs, output_std3$somatic_cc.lhs, output_std3$social_cc.lhs, output_std3$thought_cc.lhs, output_std3$att_cc.lhs, output_std3$agg_cc.lhs, output_std3$rulebr_cc.lhs, output_std3$anxdep_mod.lhs, output_std3$withdep_mod.lhs, output_std3$somatic_mod.lhs, output_std3$social_mod.lhs, output_std3$thought_mod.lhs, output_std3$att_mod.lhs, output_std3$agg_mod.lhs, output_std3$rulebr_mod.lhs, output_std3$anxdep_cpl.lhs, output_std3$withdep_cpl.lhs, output_std3$somatic_cpl.lhs, output_std3$social_cpl.lhs, output_std3$thought_cpl.lhs, output_std3$att_cpl.lhs, output_std3$agg_cpl.lhs, output_std3$rulebr_cpl.lhs)

cols_rhs <- gsub(".lhs", ".rhs", cols_lhs)

output_table$rhs <- c(output_std3$anxdep_cc.rhs, output_std3$withdep_cc.rhs, output_std3$somatic_cc.rhs, output_std3$social_cc.rhs, output_std3$thought_cc.rhs, output_std3$att_cc.rhs, output_std3$agg_cc.rhs, output_std3$rulebr_cc.rhs, output_std3$anxdep_mod.rhs, output_std3$withdep_mod.rhs, output_std3$somatic_mod.rhs, output_std3$social_mod.rhs, output_std3$thought_mod.rhs, output_std3$att_mod.rhs, output_std3$agg_mod.rhs, output_std3$rulebr_mod.rhs, output_std3$anxdep_cpl.rhs, output_std3$withdep_cpl.rhs, output_std3$somatic_cpl.rhs, output_std3$social_cpl.rhs, output_std3$thought_cpl.rhs, output_std3$att_cpl.rhs, output_std3$agg_cpl.rhs, output_std3$rulebr_cpl.rhs)

output_table$op <- "~"


cols_est.std <- gsub(".rhs", ".est.std", cols_rhs)


output_table$est.std <- c(output_std3$anxdep_cc.est.std, output_std3$withdep_cc.est.std, output_std3$somatic_cc.est.std, output_std3$social_cc.est.std, output_std3$thought_cc.est.std, output_std3$att_cc.est.std, output_std3$agg_cc.est.std, output_std3$rulebr_cc.est.std, output_std3$anxdep_mod.est.std, output_std3$withdep_mod.est.std, output_std3$somatic_mod.est.std, output_std3$social_mod.est.std, output_std3$thought_mod.est.std, output_std3$att_mod.est.std, output_std3$agg_mod.est.std, output_std3$rulebr_mod.est.std, output_std3$anxdep_cpl.est.std, output_std3$withdep_cpl.est.std, output_std3$somatic_cpl.est.std, output_std3$social_cpl.est.std, output_std3$thought_cpl.est.std, output_std3$att_cpl.est.std, output_std3$agg_cpl.est.std, output_std3$rulebr_cpl.est.std)

cols_se <- gsub(".est.std", ".se", cols_est.std)


output_table$se <-  c(output_std3$anxdep_cc.se, output_std3$withdep_cc.se, output_std3$somatic_cc.se, output_std3$social_cc.se, output_std3$thought_cc.se, output_std3$att_cc.se, output_std3$agg_cc.se, output_std3$rulebr_cc.se, output_std3$anxdep_mod.se, output_std3$withdep_mod.se, output_std3$somatic_mod.se, output_std3$social_mod.se, output_std3$thought_mod.se, output_std3$att_mod.se, output_std3$agg_mod.se, output_std3$rulebr_mod.se, output_std3$anxdep_cpl.se, output_std3$withdep_cpl.se, output_std3$somatic_cpl.se, output_std3$social_cpl.se, output_std3$thought_cpl.se, output_std3$att_cpl.se, output_std3$agg_cpl.se, output_std3$rulebr_cpl.se)

output_table$model <- paste0(output_table$lhs, " ", output_table$op, " ", output_table$rhs)

out <- output_table[ , c(1, 5, 6, 7, 8)] 

out2 <- out
rownames(out2) <- out2$model 
out2$model <- NULL

out2 <- round(out2, digits = 3)

# separate into temporal direction (brain --> beh; beh --> brain)

BrtoBeh <- grep("cpl_t1|mod_t1|cc_t1", rownames(out2), value = T)

out_BrtoBeh <- out2[rownames(out2) %in% BrtoBeh, ]
out_BehtoBr <- out2[rownames(out2) %notin% BrtoBeh, ]


out2 <- rbind(out_BrtoBeh, out_BehtoBr)


rownames(out2) <- gsub("cc", "Modularity", # NB this is because the data was inverted 
                             gsub("mod", "Clustering coefficient", # as above
                                  gsub("cpl", "Characteristic path length", 
                                            gsub("anxdep", "Anxious/Depressed",
                                                 gsub("withdep", "Withdrawn/Depressed", 
                                                      gsub("somatic", "Somatic complaints",
                                                           gsub("social", "Social problems",
                                                                gsub("thought", "Thought problems",
                                                                     gsub("att", "Attention problems",
                                                                          gsub("agg", "Aggressive behaviors", 
                                                                               gsub("rulebr", "Rule-breaking behaviors",
                                                                                    gsub("_t2", " T2", 
                                                                                         gsub("_t1", " T1",
                                                                                         rownames(out2))))))))))))))








write.csv(out2, paste0(tab, "results_syndromes.csv"))

out2
