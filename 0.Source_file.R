########
# 0. Source file
########
# Project: Exploring the longitudinal associations of functional network 
# connectivity and psychiatric symptom changes in youth
# from the general population
# Cohort: Generation R Study
# Author script: Lorenza Dall'Aglio (l.dallaglio@erasmusmc.nl) 
# License CC 4.0 international


### set environment ###

indata <- "SETWHEREYOURDATAIS"
fig <- "SETWHEREYOUWANTOUTPUTFIGURES"
tab <- "SETWHEREYOUWANTOUTPUTTABS"


setwd(indata)

x <- c("dplyr", "data.table", "ggplot2", "cowplot", "ggpubr", "beanplot", 
       "tidyverse", "naniar", "gtsummary", "flextable", "RColorBrewer","foreign",
       "sjPlot", "lavaan", "bestNormalize", "MASS", "lavaan", "fastDummies", "semTable",
       "tableone", "lsr")

lapply(x, require, character.only = T) 


### useful functions ###

`%notin%` <- Negate(`%in%`) # negative %in%


# to compare for non-response analysis (by Rosa Mulder, PhD)

compare_them <- function(data1,data2,vars) {
  
  names(vars) <- vars
  
  stat <- lapply(vars, function(x) { 
    if(is.numeric(data1[,x])){
      result <- as.numeric(unlist(t.test(data1[,x],data2[,x]))[1:3])
    }else if(is.factor(data1[,x])){
      d1 <- cbind(data1[,x],"d1")
      d2 <- cbind(data2[,x],"d2")
      d  <- as.data.frame(rbind(d1,d2))
      colnames(d) <- c("var","set")
      result <- as.numeric(unlist(chisq.test(d$var, d$set))[1:3])
    }else{
      result <- c(NA,NA,NA)                   #please check if NAs occur, should not happen
    }
  })
  names  <- vars
  save   <- as.data.frame(stat[1:length(stat)], c("t_or_X_stat","df","p"))
  save   <- tFrame(save)
  
}

