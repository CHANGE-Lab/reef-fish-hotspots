# -----------------------------------------------------------
# Quantifying fish-derived nutrient hotspots across reefscapes
# By Noelle K. Helder
# Updated August 23, 2024

# Manuscript DOI: TBD 
# -----------------------------------------------------------
# Script objectives
#  This script fits LMMs to estimated fish nutrient supply (N & P) as a function of 
#     habitat complexity (relief + VRM + VRM dev)
# Outputs: 
#   - Model Table: Table 1: N & P ~ complexity - "output/supply.complexity.lmm.table.doc"
#   - Visual: Fig 3: N & P ~ complexity - "output/Fig-3.png"
#   - Visual: LMM interaction terms - "output/Fig-4.png")
# --------------------------------------


# Packages: 
library(tidyverse) # general
library(sjPlot) # plot effects 
library(stargazer) # model table outputs
library(lme4) # lme
library(DHARMa) # residuals 
library(ggeffects) # visualize
library(cowplot) # visualize 
library(svglite)


# --------------------------------------
# Load data
# Habitat complexity metrics derived from DEMs
complexity <- read.csv("data/reef_structure.csv")

# Nutrient supply
mean_supply <- read.csv("data/predicted_nutrient_supply.csv") %>% 
  dplyr::select(-'X') %>%  
  ungroup() %>% 
  group_by(site, grid) %>% 
  # calculate mean supply per plot
  dplyr::summarize(mean.N=mean(sum.N.m2, na.rm=TRUE),
            sd.N=sd(sum.N.m2, na.rm=TRUE),
            mean.P=mean(sum.P.m2, na.rm=TRUE),
            sd.P=sd(sum.P.m2, na.rm=TRUE)) %>% 
  left_join(complexity) %>% 
  drop_na()

rm(complexity)

# --------------------------------------
# Prep data for modelling: 
# --------------------------------------
# Remove two sites that are not included in this analysis
mean_supply <- mean_supply %>% ungroup () %>% 
  dplyr::filter(site != "PR") %>% 
  dplyr::filter(site != "PC") 

# Fix names
df <- mean_supply
df$relief <- df$height_diff



# Scale and center continuous complexity: 
df$vrm_1cm_sc <- scale(df$vrm_1cm,center = T, scale = T)
df$vrm_dev_sc <- scale(df$vrm_dev,center = T, scale = T)
df$relief_sc <- scale(df$relief,center = T, scale = T)

# --------------------------------------
# Use Linear Mixed Effects Model to quantify relationship between 
#     estimated nitrogen from fish and intra-reef structural complexity (relief, vrm, vrm_dev)
#     Site (reefscape) is a random effect. 

# We test for interactions between relief and both VRM metrics (2 total interaction terms), 
# and include a polynomial term to test for non-linearity. 

# We use backwards stepwise selection to identify a 'top' model based on delta AIC and 
# the principle of parsimony when required. 

# Model selection tables are shown in scripts 01b_P_complexity_LMM_TableS2.R and
# 01c_N_complexity_LMM_TableS3.R. Only top models are visualized here. 
# --------------------------------------

N.best <- lmer(log(mean.N) ~ poly(relief,2)*vrm_dev_sc + 
                   (vrm_1cm_sc) +
                   (1|site), REML=FALSE, data=df)



P.best <- lmer(log(mean.P) ~ poly(relief,2)*vrm_dev_sc + vrm_1cm_sc +
                   (1|site), REML=FALSE, data=df)

# --------------------------------------


# --------------------------------------
# Table 1 (and Table S2, Table S3) 
# Model table for supply (Nitrogen, Phosphorous) and structural complexity
# --------------------------------------
# summary(lmm.best)
# summary(P.best)

# Table 1
stargazer(N.best, P.best, 
          type="html",
          ci= TRUE, ci.separator = " - ",
          column.labels = c("N", "P"),
          digits = 2,
          star.cutoffs = c(0.05),
          digit.separator = "",
          covariate.labels=c("Relief", "Relief(2)","VRM Dev", "VRM (1cm)", "Relief*VRM Dev",
                             "Relief(2)*VRM Dev"),
          omit.table.layout = "sn",
          out="output/supply.complexity.lmm.table.doc")


# --------------------------------------
# Create Fig 3 - Relationships between supply and complexity 
# Goal: Visualize  model predictions over raw data values with ggpredict
# --------------------------------------

# Define a base theme for consistency
base_theme <- theme_classic() +
  theme(
    text = element_text(size = 8),   # Adjust text size to 8 points
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )

# Nitrogen supply-complexity relationships
pred.N.rel <- ggpredict(N.best, terms = c("relief [all]"), ci.lvl=0.95)  


# Nitrogen-relief relationship
N.p1 <- ggplot(pred.N.rel) +
  # slope 
  geom_line(aes(x = x, 
                y = log(predicted))) +    
  # standard errors 
  geom_ribbon(aes(x = x, ymin = log(predicted) - std.error, 
                         ymax = log(predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +
  # add raw data (scaled)
  geom_point(data = df,                      
          aes(x = relief, y = log(mean.N)), alpha = 0.5) + 
  # labels
  labs(x = "Relief (m)", 
       y=expression(ln~N~Supply~(mg~m^-2~day^-1))) + 
  base_theme

N.p1


# Nitrogen-VRM (1cm) relationship 
pred.N.vrm <- ggpredict(N.best, terms = c("vrm_1cm_sc [all]"), ci.lvl=0.95)  # this gives overall predictions for the model

N.p2 <- ggplot(pred.N.vrm) + 
  geom_line(aes(x = x, y = log(predicted))) +          # slope
  geom_ribbon(aes(x = x, ymin = log(predicted) - std.error, 
                  ymax = log(predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = df,                      # adding the raw data (scaled values)
            aes(x = vrm_1cm_sc, y = log(mean.N)), alpha = 0.5) + 
  labs(x = "VRM (1cm)", 
       y=expression(ln~N~Supply~(mg~m^-2~day^-1))) + 
  base_theme

N.p2

# Nitrogen-VRM deviation relationship
pred.N.vrmdev <- ggpredict(N.best, terms = c("vrm_dev_sc [all]"), ci.lvl=0.95)  

N.p3 <- ggplot(pred.N.vrmdev) + 
  geom_line(aes(x = x, y = log(predicted))) +         
  geom_ribbon(aes(x = x, ymin = log(predicted) - std.error, 
                  ymax = log(predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = df,                      # adding the raw data (scaled values)
             aes(x = vrm_1cm_sc, y = log(mean.N)), alpha = 0.5) + 
  labs(x = "VRM Deviation", 
       y=expression(ln~N~Supply~(mg~m^-2~day^-1))) + 
  base_theme
N.p3


# Phosphorous-supply-complexity relationships
pred.P.rel <- ggpredict(P.best, terms = c("relief [all]"), ci.lvl=0.95)  


# Phosphorous-relief relationship
P.p1 <- ggplot(pred.P.rel) +
  # slope 
  geom_line(aes(x = x, 
                y = log(predicted))) +    
  # standard errors 
  geom_ribbon(aes(x = x, ymin = log(predicted) - std.error, 
                  ymax = log(predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +
  # add raw data (scaled)
  geom_point(data = df,                      
             aes(x = relief, y = log(mean.P)), alpha = 0.5) + 
  # labels
  labs(x = "Relief (m)", 
       y=expression(ln~P~Supply~(mg~m^-2~day^-1))) + 
  base_theme

P.p1


# Phosphorous-VRM (1cm) relationship 
pred.P.vrm <- ggpredict(P.best, terms = c("vrm_1cm_sc [all]"), ci.lvl=0.95)  # this gives overall predictions for the model

P.p2 <- ggplot(pred.P.vrm) + 
  geom_line(aes(x = x, y = log(predicted))) +          # slope
  geom_ribbon(aes(x = x, ymin = log(predicted) - std.error, 
                  ymax = log(predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = df,                      # adding the raw data (scaled values)
             aes(x = vrm_1cm_sc, y = log(mean.P)), alpha = 0.5) + 
  labs(x = "VRM (1cm)", 
       y=expression(ln~P~Supply~(mg~m^-2~day^-1))) + 
  base_theme

P.p2

# Phosphorous-VRM deviation relationship
pred.P.vrmdev <- ggpredict(P.best, terms = c("vrm_dev_sc [all]"), ci.lvl=0.95)  # this gives overall predictions for the model

P.p3 <- ggplot(pred.P.vrmdev) + 
  geom_line(aes(x = x, y = log(predicted))) +          # slope
  geom_ribbon(aes(x = x, ymin = log(predicted) - std.error, 
                  ymax = log(predicted) + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = df,                      # adding the raw data (scaled values)
             aes(x = vrm_1cm_sc, y = log(mean.P)), alpha = 0.5) + 
  labs(x = "VRM Deviation", 
       y=expression(ln~P~Supply~(mg~m^-2~day^-1))) + 
  base_theme

P.p3

# Combine plots for supply-structure model - Ecol. Apps submission 
# Size and Proportion
# • Images should be of adequate resolution (300-600 dpi).
# • When possible, submit figures in the size intended for publication in the typeset PDF of your
#   manuscript, using single-column width sizing of 8.5 cm for figures when possible. The maximum
#   allowable size is 18 cm wide by 24 cm tall (to accommodate a caption underneath).
# • All text must be sized between 6 and 10 points when the image is sized for publication.
fig3 <- cowplot::plot_grid(N.p1, N.p2, N.p3, 
                                    P.p1, P.p2, P.p3, 
                                    nrow=2, 
                                    labels=c('(a)', '(b)', '(c)',
                                             '(d)', '(e)', '(f)'),
                                    label_size = 8) # Adjust label_size as needed)

# Save the combined plot
ggsave(
  filename = "output/Fig-3.jpg",  # Filename for the saved plot
  plot = fig3,                    # The plot to save
  width = 18,                             # Width in cm
  height = 16,                             # Height in cm (adjust as needed)
  units = "cm",                           # Units for width and height
  dpi = 500                                # Resolution in dots per inch
)
# --------------------------------------
# Create Fig 4:
# Visualize model interaction terms 
# --------------------------------------
# Interaction between VRM and Relief: Positive effect of VRM  in low-relief habitat

P.int <- plot_model(P.best, type = "pred", 
                    terms = c("vrm_dev_sc","relief")) +
   scale_color_sjplot("viridis", labels=c("Low", "Medium", "High")) +
   scale_fill_sjplot("viridis") +
   ggplot2::labs(colour = "Relief") +
   labs(x="VRM Deviation", y=expression(ln~Mean~P~Supply~(mg~m^-2~day^-1))) +
   labs(title = "") +
   base_theme
P.int


N.int <- plot_model(N.best, type = "pred", 
                     terms = c("vrm_dev_sc","relief")) +
    scale_color_sjplot("viridis", labels=c("Low", "Medium", "High")) +
    scale_fill_sjplot("viridis") +
    ggplot2::labs(colour = "Relief") +
    labs(x="VRM Deviation", y=expression(ln~Mean~N~Supply~(mg~m^-2~day^-1))) +
    labs(title = "") +
    base_theme

# Extract the legend from one of the plots (e.g., N.int)
legend <- get_legend(N.int + theme(legend.position = "left"))

# Combined figure for manuscript (Figure 4)
fig4 <- cowplot::plot_grid(
  N.int+ theme(legend.position="none"),
  P.int + theme(legend.position="none"),
  legend,
  align='vh',
  labels=c('(a)', '(b)'),
  label_size = 8,
  nrow=1,
  rel_widths = c(0.45, 0.45, 0.10))
fig4

# Save the combined plot
ggsave(
  filename = "output/Fig-4.jpg",  # Filename for the saved plot
  plot = fig4,                    # The plot to save
  width = 18,                             # Width in cm
  height = 8,                             # Height in cm (adjust as needed)
  units = "cm",                           # Units for width and height
  dpi = 300                                # Resolution in dots per inch
)
