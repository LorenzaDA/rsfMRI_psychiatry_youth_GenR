#######################
# 11. regression plots
#######################
# as asked by reviewers, we're showing regression plots
# for nominally sign. associations. 
# of note,this is only for complete cases
# and with lm instead of clpms
# so they don't necessarily map exactly to clpms results

# import data

source("/RELEVANT/PATH/0.Source_file.R")

dd <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds"))

# subset to complete cases for this to work 

dd <- dd[complete.cases(dd), ]


### mdt5 --> ext

fit.2 <- lm(ext_t2 ~ sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + ext_t1 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

# plot 
a <- ggplot(dd, aes(x = mdt5_t1, y = residual2))+
  geom_point(aes(colour = "#009E73", fill = "#009E73"))+
  geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
  xlab("MDT state 5 T1") + 
  ylab("Externalizing pr. T2")+ 
  theme_classic()+
  theme(legend.position = "none")
  
  
### int --> mdt3

fit.2 <- lm(mdt3_t2 ~ mdt3_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(b <- ggplot(dd, aes(x = int_t1, y = residual2))+
  geom_point(aes(colour = "#009E73", fill = "#009E73"))+
  geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm") + 
  xlab("Internalizing pr. T1") + 
  ylab("MDT state 3 T2")+ 
  theme_classic()+
  theme(legend.position = "none"))


### MDT5 --> attention

fit.2 <- lm(att_t2 ~ att_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(c <- ggplot(dd, aes(x = mdt5_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
    xlab("MDT state 5 T1") + 
    ylab("Attention pr. T2")+ 
    theme_classic()+
    theme(legend.position = "none"))



### MDT5 --> aggression 

fit.2 <- lm(agg_t2 ~ agg_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(d <- ggplot(dd, aes(x = mdt5_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
    xlab("MDT state 5 T1") + 
    ylab("Aggression T2")+ 
    theme_classic()+
    theme(legend.position = "none"))


### Rule breaking → MDT2 

fit.2 <- lm(mdt2_t2 ~ mdt2_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(e <- ggplot(dd, aes(x = rulebr_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
    xlab("Rule-breaking T1") + 
    ylab("MDT state 2 T2")+ 
    theme_classic()+
    theme(legend.position = "none"))


### somatic --> MDT3

fit.2 <- lm(mdt3_t2 ~ mdt3_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(f <- ggplot(dd, aes(x = somatic_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
    xlab("Somatic pr. T1") + 
    ylab("MDT state 3 T2")+ 
    theme_classic()+
    theme(legend.position = "none"))



### thought problems → MDT3

fit.2 <- lm(mdt3_t2 ~ mdt3_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(g <- ggplot(dd, aes(x = thought_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
    xlab("Thought pr. T1") + 
    ylab("MDT state 3 T2")+ 
    theme_classic()+
    theme(legend.position = "none"))


### attention problems → MDT4 

fit.2 <- lm(mdt4_t2 ~ mdt4_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(h <- ggplot(dd, aes(x = att_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm") + 
    xlab("Attention pr. T1") + 
    ylab("MDT state 4 T2")+ 
    theme_classic()+
    theme(legend.position = "none"))


### Thought problems → NT

fit.2 <- lm(nt_t2 ~ nt_t1 + sex + mat_edu + agediff_t1 + agediff_t2 + agemri_t1 + agemri_t2 + puberty + ethn, data = dd)

dd$residual2 <- residuals(fit.2)

(i <- ggplot(dd, aes(x = thought_t1, y = residual2))+
    geom_point(aes(colour = "#009E73", fill = "#009E73"))+
    geom_smooth(aes(colour = "#E69F00", fill = "E69F00"), method = "lm")+ 
    xlab("Thought pr. T1") + 
    ylab("NT T2")+ 
    theme_classic()+
    theme(legend.position = "none"))


# brain to beh: a, d, c
# beh to brain:b, e, f, g, h, i

all_brtobeh <- ggarrange(a, d, c,
                  ncol = 3, nrow = 1)

ggsave("Rel_BrtoBeh_rsfmri_psych.png", path = fig)


all_behtobr <- ggarrange(b, e, f, g, h, i,
                         ncol = 3, nrow = 2)


ggsave("Rel_BehtoBr_rsfmri_psych.png", path = fig)
