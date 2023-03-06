####################################
# 3. Analyses by sex - static rsfMRI 
####################################
# project: Exploring the longitudinal associations of functional network connectivity
# and psychiatric symptom changes in youth
# cohort: GenR (focus at @9 and @13)
# author script: Lorenza Dall'Aglio (l.dallaglio@erasmusmc.nl) 
# license CC 4.0 international


#####
# Set environment
####

rm(list=ls())

source("/PATH/WHERE/SOURCE/FILE/IS/0.Source_file.R")

tab <- "/OUTPUT/PATH/"

# load data 

dd <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds"))


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
# based on https://lavaan.ugent.be/tutorial/cat.html

### Make ordered variables, continuous/numeric ###
# maternal education only here

dd$mat_edu_num <- as.numeric(dd$mat_edu)


### dummy coding for binary/nominal variables ### 

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
# general model specifications which will be adapted at a later stage

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
int_t2 ~ ethn_other + ethn_european + mat_edu_num

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
ext_t2 ~ ethn_other + ethn_european + mat_edu_num

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



######
# Model identification and results
######

# group by sex

fits_all <- sapply(specifications, function(x) sem(x, data = dd2, missing = "fiml", fixed.x = F, group = "sex"))

# group by sex setting the coefficients as equal to test for sex differences

fits_all_reg <- sapply(specifications, function(x) sem(x, data = dd2, missing = "fiml", fixed.x = F, group = "sex", group.equal = "regressions"))


comp1 <- lavTestLRT(fits_all$cpl_int, fits_all_reg$cpl_int)
comp2 <- lavTestLRT(fits_all$cc_int, fits_all_reg$cc_int)
comp3 <- lavTestLRT(fits_all$mod_int, fits_all_reg$mod_int)

comp4 <- lavTestLRT(fits_all$cpl_ext, fits_all_reg$cpl_ext)
comp5 <- lavTestLRT(fits_all$cc_ext, fits_all_reg$cc_ext)
comp6 <- lavTestLRT(fits_all$mod_ext, fits_all_reg$mod_ext)


comparison_tab <- rbind(comp1, comp2, comp3, comp4, comp5, comp6)

comparison_tab
# no differences for any


## clean up the table
# change the rownames so they can be directly used in the supplement
rownames(comparison_tab) <- gsub("cc_", "Modularity - ", # NB this is because the data was inverted 
                            gsub("mod_", "Clustering coefficient - ", # as above
                                 gsub("cpl_", "Characteristic path length - ", 
                                           gsub("int", "Internalizing problems",
                                                gsub("ext", "Externalizing problems",
                                                               gsub("fits_all", "fit ", 
                                                                   gsub("\\$", " ", 
                                                                    gsub("_reg", "(equal estimates) ",
                                                               rownames(comparison_tab)))))))))


# round to 3 decimals 

comparison_tab <- round(comparison_tab, digits = 3)


# save table

write.csv(comparison_tab, paste0(tab, "comparison_across_sexes.csv"))


