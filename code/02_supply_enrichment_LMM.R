# -----------------------------------------------------------
# Quantifying fish-derived nutrient hotspots across reefscapes
# By Noelle K. Helder

# Updated August 27 2024

# Manuscript DOI: submitted to Ecol Apps  
# -----------------------------------------------------------
# Script objectives
#  This script fits LMMs to macroalgal nutrient enrichment (tissue %N and %P) as a function 
#     of estimated N/P nutrient supply from fishes  

# Outputs: 
#   - Model Table: Table 2: %N ~ N supply - "output/Table2.doc"
#   - Model Table: Table 3: %P ~ P supply - "output/Table3.doc"
#   - Visual: Fig 5: %N/%P ~ N/P supply - "output/Fig-5.png"
#   - Visual: Fig S2: %N/%P ~ N/P supply by site - "output/Fig-S2.png"
# --------------------------------------

# Load Packages: 
library(MuMIn)
library(tidyverse)
library(ggeffects)
library(lme4)
library(cowplot)
# install.packages("svglite")
library(svglite)

# --------------------------------------
# Load  data
# --------------------------------------
# 1) Intra-reef structural complexity data from photogrammetry DEMs
predictors <- read.csv(
  # "~/UofA/Thesis/r_projects/thesis_chapters/fish_habitat_models/data/predictors.sub.csv")
  "data/reef_structure.csv")

# 2) Estimated nutrient supply (mg/m2/day) from all fishes per plot combined
mean_supply <- read.csv("data/predicted_nutrient_supply.csv") %>% 
  dplyr::select(-'X') %>%  
  ungroup() %>% 
  group_by(site, grid) %>% 
  # calculate average and SD supply (Nitrogen + phosphorous) across visits for each plot 
  dplyr::summarize(mean.N=mean(sum.N.m2, na.rm=TRUE),
                   sd.N=sd(sum.N.m2, na.rm=TRUE),
                   mean.P=mean(sum.P.m2, na.rm=TRUE),
                   sd.P=sd(sum.P.m2, na.rm=TRUE)) %>% 
  left_join(predictors) %>% 
  ungroup () %>% 
  # PC and PR are 2 reefscapes that are not a part of this analysis. 
  dplyr::filter(site != "PR") %>% 
  dplyr::filter(site != "PC") %>% 
  drop_na()


# 3) Macroalgal tissue content data measured from samples at BASL
#     Data has been cleaned up from raw BASL format to match up with each corresponding plot. 
NDR.fp <- read.csv("data/input/NDR.fp.csv") # 151 observations of 12 vars for NDR
CR.fp <- read.csv("data/input/CR.fp.csv") # 121 obs of 12 vars for CR


# 4) Metadata - corresponding information about each plot survey
metadata <- read.csv("data/input/metadata.csv")


# Prep and clean data for analysis
# Create separate Nitrogen and Phosphorous datasets from macroalgae sample data
# 96 N samples at CR
CR.N <- CR.fp %>% dplyr::select(-'P', -'date') %>% 
  drop_na() 
# 25 P samples at CR
CR.P <- CR.fp %>% dplyr::select(site, fp, transect, grid, distance, P) %>% 
  drop_na()
# 110 samples at NDR
NDR.N <- NDR.fp %>% dplyr::select(-'P', -'date')  %>% 
  drop_na()
# 41 P samples at NDR
NDR.P <- NDR.fp %>% dplyr::select(site, fp, transect, grid, distance, P) %>% 
  drop_na() 

# Combine for a total N dataframe: 
96+110 # total N observations: 206
dat.N <- NDR.N %>% rbind(CR.N) %>% left_join(mean_supply)

# Combine for a total P dataframe: 
25+41 # total P observations: 66
dat.P <- NDR.P %>% rbind(CR.P) %>% left_join(mean_supply)

# Clean up workspace: 
rm(CR.fp, NDR.fp)

# --------------------------------------
# Modeling the relationship between predicted fish-derived nitrogen and macroalgal enrichment

# Linear Mixed Effects Model including a quadratic term to test 
# for non-linear relationships. 
# Plot (focal point, or 'fp') nested within reefscape as random effects to control for sampling design. 

# Units: 
# N/P are %Nitrogen and %Phosphorous
# mean.N/P are mg/m2/day

# Data includes all samples (N=206, P=66; Reefscapes = Carysfort and North Dry Rocks)
# --------------------------------------
# Remove missing data
dat.N <- dat.N %>% drop_na()

# Total model (all sampling sites combined)
N.mod <- lmer(N ~ poly(mean.N,2) + (1|fp/site), data=dat.N)
N.mod2 <- lmer(N ~ mean.N + (1|fp/site), data=dat.N)

AIC(N.mod, N.mod2)
# Explore random effects structure
# N.mod2 <- lmer(N ~ mean.N + (1|fp/site), data=dat.N)
# N.mod3 <- lmer(N ~ mean.N + (1|grid/fp/site), data=dat.N)
# N.mod4 <- lmer(N ~ mean.N + (1|grid/site), data=dat.N)
# N.mod5 <- lmer(N ~ poly(mean.N, 2) + (1|grid/site), data=dat.N)
# N.mod6 <- lmer(N ~ poly(mean.N,2)+ (1|grid/site), data=dat.N)
# N.mod7 <- lmer(N ~ poly(mean.N,2) + (1|grid/site) + (1|fp), data=dat.N)

# Compare AIC scores
# AIC(N.mod, N.mod2, N.mod3, N.mod4, N.mod5, N.mod6, N.mod7)

# Select lowest AIC scores and then choose simplest model. 5, 6, and 7 are all basically equal, 
# but mod6 is the simplest. 
# summary(N.mod6)
# r.squaredGLMM(N.mod6)

# Check for assumptions of best model
# plot(N.mod)
# abline(0,0)

# Final model results: 
# summary(N.mod)

# Clean up
# rm(N.mod2, N.mod3, N.mod4, N.mod5, N.mod6, N.mod7, pred.mm, N.combined)

# Table 2
stargazer(N.mod, 
          type="html",
          ci= TRUE, ci.separator = " - ",
          column.labels = c("Macroalgal %N"),
          digits = 2,
          star.cutoffs = c(0.05),
          digit.separator = "",
          covariate.labels=c("N mg m(2)day(-1)", "N(2) mg m(2)day(-1"),
          omit.table.layout = "sn",
          out="output/Ecol-Apps/Table2.doc")
# -----------------------------------------------------------
# Modeling the relationship between predicted fish-derived phosphorous and macroalgal enrichment
# -----------------------------------------------------------
dat.P <- dat.P %>% drop_na()

# Base model
P.mod <- lmer(P ~ mean.P + (1|grid), data=dat.P)

# Check for non-linearity
P.poly <- lmer(P ~ poly(mean.P,2) + (1|grid), data=dat.P)

# Polynomial model is a better fit (lower AIC)
AIC(P.mod, P.poly)
r.squaredGLMM(P.poly)

# Table 3
stargazer(P.poly, 
          type="html",
          ci= TRUE, ci.separator = " - ",
          column.labels = c("Macroalgal %P"),
          digits = 2,
          star.cutoffs = c(0.05),
          digit.separator = "",
          covariate.labels=c("P mg m(2)day(-1)", "P(2) mg m(2)day(-1"),
          omit.table.layout = "sn",
          out="output/Ecol-Apps/Table3.doc")
# -----------------------------------------------------------
# Create Figure 5
# Visualize relationships between fish nutrient supply and benthic enrichment 
# -----------------------------------------------------------

# Define a base theme for consistency
base_theme <- theme_classic() +
  theme(
    text = element_text(size = 8),   # Adjust text size to 8 points
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )

# Visualize %N as a function of N supply
# get model predictions
pred.mm <- ggpredict(N.mod, terms = c("mean.N [all]")) 
N.combined <- ggplot(pred.mm) + 
  geom_line(aes(x = x, y = (predicted))) +          
  geom_ribbon(aes(x = x, 
                  ymin = (predicted) - std.error, 
                  ymax = (predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  
  geom_jitter(data=dat.N, aes(x=mean.N, 
                              y=N,
                              color=site,
                              shape=site), 
              size=2, width=0.1,height=0.1,alpha=0.8) +
  labs(x=expression(N~Supply~(mg~m^-2~day^-1)), 
       y="Macroalgae %N") +
  base_theme +
  scale_color_grey() +
  guides(shape=guide_legend(title="Reefscape"),
         color=guide_legend(title="Reefscape"))

N.combined

# Phosphorous
pred.P <- ggpredict(P.poly, terms = c("mean.P [all]"))

P.combined <- ggplot(pred.P) + 
  geom_line(aes(x = x, y = (predicted))) +    
  geom_ribbon(aes(x = x, ymin = (predicted) - std.error, 
                  ymax = (predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +
  geom_jitter(data=dat.P, aes(x=mean.P, y=P,shape=site, color=site), 
              size=2, 
              width=0.05,
              height=0.05,
              alpha=0.8) + 
  labs(x=expression(P~Supply~(mg~m^-2~day^-1)), 
       y="Macroalgae %P") +
  base_theme + 
  scale_color_grey()
P.combined

# Combine plots:
fig5 <- cowplot::plot_grid(
  N.combined + theme(legend.position="none"),
  P.combined + theme(legend.position="none"),
  align='vh',
  labels=c('(a)', '(b)'),
  label_size = 8,
  nrow=1,
  rel_widths = c(0.5, 0.5))

# Note on rel_widths: Allocate widths to plots and legend. 
# For example, c(0.45, 0.45, 0.10) means each plot takes 45% of the width, 
# and the legend takes 10%. Adjust these ratios to fit your layout.
# Save the combined plot
ggsave(
  filename = "output/Ecol-Apps/Fig-5.jpg",  # Filename for the saved plot
  plot = fig5,                    # The plot to save
  width = 18,                             # Width in cm
  height = 8,                             # Height in cm (adjust as needed)
  units = "cm",                           # Units for width and height
  dpi = 500                                # Resolution in dots per inch
)


# -----------------------------------------------------------
# Supplementary Info - Parsing out hotspots by site
# -----------------------------------------------------------
# Phosphorous 
# One model for NDR + CR separately: separate dataframes
NDR.P <- NDR.P %>% left_join(mean_supply) %>% drop_na()
CR.P <- CR.P %>% left_join(mean_supply) %>%  drop_na()


# NDR Phosphorous model: poly model is better than linear: 
P.mod <- lmer(P ~ mean.P + (1|grid), data=NDR.P)
P.poly <- lmer(P ~ poly(mean.P,2) + (1|grid), data=NDR.P)
summary(P.mod)
summary(P.poly)

AIC(P.mod, P.poly)
r.squaredGLMM(P.poly)

# Visualize NDR phosphorous trends
pred.P <- ggpredict(P.poly, terms = c("mean.P [all]")) 
NDR.P.plot <- ggplot(pred.P) + 
  geom_line(aes(x = x, y = (predicted))) +    
  geom_ribbon(aes(x = x, ymin = (predicted) - std.error, 
                  ymax = (predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +
  
  geom_point(data = NDR.P, shape=17, size=1,               
             aes(x = mean.P, y = P)) +
  labs(x=expression(P~Supply~(mg~m^-2~day^-1)), 
       y="Macroalgae %P") +
  base_theme +
  scale_color_grey() +
  guides(shape=guide_legend(title="Site")) +
  theme(legend.position = "none")

NDR.P.plot


# Carysfort (CR1)
P.cr <- lm(P ~ mean.P, data=CR.P)
P.poly <- lm(P ~ poly(mean.P,2), data=CR.P)

AIC(P.cr, P.poly)

# Carysfort phosphorous trends 
pred.P.cr <- ggpredict(P.cr, terms = c("mean.P")) 
P.cr.plot <- ggplot(pred.P.cr) + 
  geom_line(aes(x = x, y = (predicted))) +    
  geom_ribbon(aes(x = x, ymin = (predicted) - std.error, 
                  ymax = (predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +
  geom_jitter(data=CR.P, aes(x=mean.P, y=P,shape=site), 
              size=1, width=0.1,height=0.1,alpha=0.8) +
 
  labs(x=expression(P~Supply~(mg~m^-2~day^-1)), 
       y="Macroalgae %P") +
  base_theme +
  scale_color_grey() +
  guides(shape=guide_legend(title="Site")) +
  theme(legend.position = "none")
P.cr.plot


# Isolate NDR fp1
dat.N.subset.ndr <- dat.N %>% dplyr::filter(site == "NDR")
dat.N.subset.ndr1 <- dat.N.subset.ndr %>% dplyr::filter(fp == "1") # df with just NDR 1
dat.N.subset.ndr <- dat.N.subset.ndr %>% dplyr::filter(fp == "2")
dat.N.subset <- dat.N %>% dplyr::filter(site != "NDR") %>% 
  rbind(dat.N.subset.ndr)


# 3 fps without NDR fp 1 (2 at CR, 1 at NDR)
N.mod.sub <- lmer(N~ mean.N + (1|grid/site), data=dat.N.subset)

# Visualize hotspot nitrogen trends (excluding NDR hotspot)
pred.sub <- ggpredict(N.mod.sub, terms = c("mean.N"))  

# Plot the predictions
subset <- ggplot(pred.sub) + 
  geom_line(aes(x = x, y = (predicted))) +          # slope
  geom_ribbon(aes(x = x, ymin = (predicted) - std.error, ymax = (predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_jitter(data=dat.N.subset, aes(x=mean.N, y=N,shape=site), 
              size=1, width=0.1,height=0.1,alpha=0.8) +
  labs(x=expression(N~Supply~(mg~m^-2~day^-1)), 
       y="Macroalgae %N") +
  base_theme +
  scale_color_grey() +
  guides(shape=guide_legend(title="Reefscape"),
         color=guide_legend(title="Reefscape")) +
  theme(legend.position="none")
subset  



# Just NDR hotspot model 
ndr1.mod <- lmer(N ~ poly(mean.N,2) + (1|grid), data=dat.N.subset.ndr1)


# Get model predictions for just NDR 1 model: use [all] to grab polynomial terms
pred.ndrfp1 <- ggpredict(ndr1.mod, terms = c("mean.N [all]")) 

# Visualize 
ndrfp1 <- ggplot(pred.ndrfp1) + 
  geom_line(aes(x = x, y = (predicted))) +          # slope
  geom_ribbon(aes(x = x, ymin = (predicted) - std.error, ymax = (predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_jitter(data=dat.N.subset.ndr1, aes(x=mean.N, y=N,
                                          color=as.factor(site), 
                                          shape=site), 
              size=1, width=0.4,height=0.1,alpha=0.8, shape=17) +
  labs(x=expression(N~Supply~(mg~m^-2~day^-1)), 
       y="Macroalgae %N") +
  base_theme +
  scale_color_grey() +
  theme(legend.position="none") 
ndrfp1

# ----------------------
# Create Figure S2
# ----------------------
# Extract the legend and set its background to transparent
legend <- get_legend(subset + theme(
  legend.position = "left",
  legend.background = element_rect(fill = "white", color = "white"), # White background for legend
  legend.box.background = element_rect(fill = "white", color = "white") # White background for legend box
))


figs2 <- cowplot::plot_grid(ndrfp1, subset, legend,
  NDR.P.plot+ theme(legend.position="none"),
  P.cr.plot + theme(legend.position="none"), NULL,
  align='vh',
  labels=c('(a)', '(b)', '', '(c)', '(d)'),
  label_size = 8,
  nrow=2,
  rel_widths = c(0.35, 0.35, 0.2)) + 
  theme(plot.background = element_rect(fill = "white", colour = NA))
figs2
?cowplot::plot_grid
?ggsave
# Save the combined plot
ggsave(
  filename = "output/Ecol-Apps/Fig-S2.png",  # Filename for the saved plot
  plot = figs2,                    # The plot to save
  width = 12,                             # Width in cm
  height = 10,                             # Height in cm (adjust as needed)
  units = "cm",                           # Units for width and height
  dpi = 300                                # Resolution in dots per inch
)




