##################
# 1c. Figures & tables
##################

rm(list=ls())

source("/PUT/PATH/WHERE/SOURCE/FILE/IS/0.Source_file.R")

rs <- readRDS(paste0(indata, "genr_data_main_rsfMRI_psych.rds"))

fig <- "/PATH/WHERE/YOU/WANT/YOUR/FIGURES"


#####
# Tables
#####


# get general summary 

summary_descriptives <- summary(rs)
sd(rs$puberty, na.rm = T) # 0.752
sd(rs$cpl_t1, na.rm = T) # 0.320
sd(rs$cpl_t2, na.rm = T) # 0.357
sd(rs$cc_t1, na.rm = T) # 0.028
sd(rs$cc_t2, na.rm = T) # 0.030
sd(rs$nt_t1, na.rm = T) # 2.87
sd(rs$nt_t2, na.rm = T) # 3.32
sd(rs$agemri_t1, na.rm = T) # 0.60
sd(rs$agemri_t2, na.rm = T) # 0.57


write.csv(summary_descriptives, paste0(fig, "descriptives_sample.csv"))


# select relevant vars

cols <- grep("mdt", names(rs), value = T) 

all_cols <- c("int_t1", "int_t2", 
              "ext_t1", "ext_t2",
              "anxdep_t1", "anxdep_t2", 
              "withdep_t1", "withdep_t2", 
              "somatic_t1", "somatic_t2",
              "social_t1", "social_t2", 
              "thought_t1", "thought_t2",
              "att_t1", "att_t2", 
              "rulebr_t1", "rulebr_t2", 
              "agg_t1", "agg_t2", "sex", "ethn",  "mat_edu", "agemri_t1", 
              "agemri_t2", "agediff_t1", "agediff_t2", "nt_t1", "nt_t2", "cpl_t1", 
              "cpl_t2", "cc_t1", "cc_t2", "mod_t1", "mod_t2",
              cols)

all_cols

temp <- rs[ , names(rs) %in% all_cols]

table1 <- tbl_summary(temp) %>% 
  modify_header(label = "**Variable**") %>% 
  italicize_levels()



# save table as word doc

table1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(table1, 
                          path = paste0(fig, "table1_static&dynamic.docx"))




######
# Data prep
#####


# select cols which have values at both t1 and t2

names(rs)

cols2 <- c("int_t1", "int_t2", 
           "ext_t1", "ext_t2","agemri_t1", 
           "agemri_t2", "agediff_t1", "agediff_t2", "nt_t1", "nt_t2", "cpl_t1", 
           "cpl_t2", "mod_t1", "mod_t2", "cc_t1", "cc_t2", cols)

cols_all <- c(cols2, "sex", "ethn", "mat_edu", "puberty", "id")

rs <- rs[ , names(rs) %in% cols_all]


# to long format 

dd_long <- rs %>% pivot_longer(
  cols = all_of(cols2), 
  names_to = c(".value", "t"),
  names_sep = "_") 

dd_long



# set list of coloblind friendly palettes

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# yellow = F0E442
# orange = E69F00 (light); D55E00 (dark)
# blue = 56B4E9 (light),  0072B2 (dark)
# green = 009E73
# pink = CC79A7
# white = FFFFFF
# black = 000000


######
# Violin plots
######


a <- ggplot(dd_long, aes(x = t, y = mdt1)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#CC79A7", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#CC79A7"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#CC79A7"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic() + 
  labs(y = "MDT1") + 
  theme(axis.title.x = element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.line.x=element_blank())



b <- ggplot(dd_long, aes(x = t, y = mdt2)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#56B4E9", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#56B4E9"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#56B4E9"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic() + 
  labs(y = "MDT2") + 
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x =element_blank())


c <- ggplot(dd_long, aes(x = t, y = mdt3)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#0072B2", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#0072B2"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#0072B2"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic() + 
  labs(y = "MDT3") + 
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank())


d <- ggplot(dd_long, aes(x = t, y = mdt4)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#009E73", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#009E73"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#009E73"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic() + 
  labs(y = "MDT4") + 
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank())



e <- ggplot(dd_long, aes(x = t, y = mdt5)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#F0E442", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#F0E442"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#F0E442"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic() + 
  labs(y = "MDT5", x = "time point") 




f <- ggplot(dd_long, aes(x = t, y = nt)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#E69F00", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#E69F00"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#E69F00"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic() + 
  labs(y = "NT", x = "time point")




g <- ggplot(dd_long, aes(x = t, y = int)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#E69F00", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#E69F00"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#E69F00"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic()




h <- ggplot(dd_long, aes(x = t, y = ext)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#D55E00", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#D55E00"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#D55E00"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic()


i <- ggplot(dd_long, aes(x = t, y = cpl)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#CC79A7", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#CC79A7"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#CC79A7"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic()






k <- ggplot(dd_long, aes(x = t, y = cc)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#009E73", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#009E73"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#009E73"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic()




l <- ggplot(dd_long, aes(x = t, y = mod)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2,
    fill = "#F0E442", 
    slab_colour = "#000000"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    colour = "#000000",
    fill = "#F0E442"
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3,
    colour = "#F0E442"
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_classic()




all1 <- ggarrange(a, b, c, d, e, f,
                  ncol = 2, nrow = 3, 
                  labels = "A")


ggsave("Distributions_mainvars_dynamic.png", path = fig)  


all2 <- ggarrange(i, k, l,
                  ncol = 3, nrow = 1)

ggsave("Distributions_mainvars_static.png", path = fig)



all3 <- ggarrange(g, h,
                  ncol = 2, nrow = 1)

ggsave("Distributions_mainvars_exposure.png", path = fig)


##########
# Results for syndrome scales 
##########

## beh --> brain 
## for dynamic connectivity 

out <- "/PATH/WHERE/YOU/SAVED/RESULTS/"

syn <- read.csv(paste0(out, "dynamic/", "results_syndromes.csv"))


syn_behtobr <- syn[49:96, ] # second 48 rows are about the behavior to brain 


### plot ####

# get CIs first

syn_behtobr$cilower <- syn_behtobr$est.std - (syn_behtobr$se * 1.96)
syn_behtobr$ciupper <- syn_behtobr$est.std + (syn_behtobr$se * 1.96)


# separate the rsfmri and syndrome info 

syn_behtobr$syndrome <- sub(" ~ ", "__", syn_behtobr$X)
temp <- strsplit(syn_behtobr$syndrome, "__", fixed = T)
part1    <- unlist(temp)[2*(1:length(syn_behtobr$syndrome))-1]
part2    <- unlist(temp)[2*(1:length(syn_behtobr$syndrome))  ]

syn_behtobr$syndrome <- part2
syn_behtobr$rsfmri <- part1


# order alphabetically the rsfmri data so you can get it in the right order
# for the y axis later 
syn_behtobr$rsfmri2<- factor((syn_behtobr$rsfmri), levels=rev(sort(unique(syn_behtobr$rsfmri), decreasing = F)))

# double check it worked 
cbind(syn_behtobr$rsfmri, syn_behtobr$rsfmri2)


# plot 
syn_behtobr_plots <- ggplot(data = syn_behtobr, aes(x = est.std, y = rsfmri2, xmin = cilower, xmax = ciupper, color = rsfmri2))+
  geom_vline(xintercept = 0, linetype ="dotted", color = "black", size=0.8) + # line in the middle that is dotted - put it before the points so it is in the background
  geom_pointrange(group = syn_behtobr$rsfmri2) + 
  geom_errorbar(width = 0.2, cex=0.2) + # whiskers on the range
  # add a vertical line at the 0 
  facet_wrap(~syndrome, ncol = 4) + # this makes us divide the plot by model (one square per model)
  scale_color_manual(values = c("MDT state 1 T2" = "#CC79A7",
                                "MDT state 2 T2" = "#56B4E9",
                                "MDT state 3 T2" = "#0072B2",
                                "MDT state 4 T2" = "#009E73",
                                "MDT state 5 T2" = "#F0E442",
                                "NT T2" = "#E69F00"
                               )) + # set the colours you want for the studies
  theme(legend.position = "bottom", strip.background = element_rect(colour = "black", fill = "white")) + 
  xlab("Standardised estimate (95% CI)") + # Label on the x axis 
  ylab("") + # Label on the y axis
  theme(axis.line.x = element_line(colour = "black", size = 0.2), 
        axis.line.y = element_line(colour = "black", size = 0.2),
        panel.background = element_blank(), 
        axis.title.x = element_text(colour = "Black", size = 9),
        axis.text = element_text(size = 7, color = "Black"), 
        text = element_text(size = 8)) + 
  theme(legend.position = "none")



effects <- ggarrange(syn_behtobr_plots, 
                     labels = "B")

all_disteffects <- ggarrange(all1, effects,
                  ncol = 1, nrow = 2)


all_disteffects


###########
# Distributions of effect sizes in the exploratory analyses
###########

median(syn$est.std) # 0.009
median(syn$se) # 0.0225


a <- ggplot(syn) + 
  geom_histogram(mapping = aes(x = abs(est.std)), color =  "#D55E00", fill = "#D55E00", bins = 10) + 
  theme(legend.position = "bottom", strip.background = element_rect(colour = "black", fill = "white")) + 
  xlab("Estimates") + # Label on the x axis 
  ylab("Count") + # Label on the y axis
  theme(axis.line.x = element_line(colour = "black", size = 0.2), 
        axis.line.y = element_line(colour = "black", size = 0.2),
        panel.background = element_blank(),
       axis.title.x = element_text(colour = "Black", size = 6),
      axis.text = element_text(size = 5, color = "Black"), 
       text = element_text(size = 6))

b <- ggplot(syn) + 
  geom_histogram(mapping = aes(x = abs(se)), color =  "#D55E00", fill = "#D55E00", bins = 10) + 
  theme(legend.position = "bottom", strip.background = element_rect(colour = "black", fill = "white")) + 
  xlab("Standard errors") + # Label on the x axis 
  ylab("Count") + # Label on the y axis
  theme(axis.line.x = element_line(colour = "black", size = 0.2), 
        axis.line.y = element_line(colour = "black", size = 0.2),
        panel.background = element_blank(), 
        axis.title.x = element_text(colour = "Black", size = 6),
        axis.text = element_text(size = 5, color = "Black"), 
        text = element_text(size = 6))


all_hists <- ggarrange(a, b,
          labels = "C")


all_fin <- ggarrange(all1, effects, all_hists,  
                             ncol = 2, nrow = 2)

library(gridExtra)

gt <- arrangeGrob(all1,
                  effects, all_hists,
                  ncol = 3, nrow = 2, 
                  layout_matrix = cbind(c(1,2), c(1, 2), c(3, 2)))

p <- ggplotify::as.ggplot(gt) 

p

ggexport(p, filename = paste0(out, "dynamic/", "Figure1_associations_exploratory.pdf"))


