# -----------------------------------------------------------
# Quantifying fish-derived nutrient hotspots across reefscapes
# By Noelle K. Helder
# Updated March 31 2023

# Manuscript DOI: TBD Ecol Apps submitted
# -----------------------------------------------------------
# This script create manuscript Figure 6 exploring hotspot spatial
# extents. 

# Methodology notes:
#     Plots show vertical relief, estimated Nitrogen supply from fish, macroalgal tissue 
#     %Nitrogen, and (for NDR1 only) A. cervicornis growth rates.
#     Plots are for focal plots (fp) at 2 reefscapes (NDR1 and CR1). Each reefscape contained
#     2 fp. Each fp consisted of 3x15 m transects starting at the fp (0m). Transects therefore 
#     extended through 4 total plots each at 0, 5, 10, and 15m. 
#
#
# Outputs:
#   NDR1 hotspot panel plot (high res)
#   CR1 hotspot panel plot (high res)
# -----------------------------------------------------------
# install.packages("PNWColors")

# Load packages
library(tidyverse)
library(PNWColors)


# Load input data
# Relief SD was calculated in the original R script from Burns et al. Here I just filtered
# the spreadsheet to just Relief and its SD, and just for NDR and CR. 
relief <- read.csv("data/input/hotspot.relief.df.csv") %>% 
  dplyr::filter(fac == "1") %>% # keep only 1 cm res estimates
  dplyr::select(-'fac') %>% 
  dplyr::select(-'height_diff') # drop extra col

# Macroalgal nutrient content per plot at NDR and CR focal points (fp), includes complexity data
N.fp <- read.csv("data/processed/N.fp.test.csv")

# Estimated Nutrient supply from fishes across all surveys (n=4-6 replicates per plot per reefscape)
supply.load <- read.csv("data/processed/pooled.predicted.nutrients.v2.csv") %>% 
  dplyr::select(-'X')

# Prep data for plotting 
# refactor for joins
N.fp$site <- as.factor(N.fp$site)
relief$site <- as.factor(relief$site)
N.fp$grid <- as.factor(N.fp$grid)
relief$grid <- as.factor(relief$grid)

# add vertical relief standard deviation (SD) to the main dataframe
N.fp <- N.fp %>% left_join(relief)

# Hotspot dataframes used in this visual: 
# North Dry Rocks (NDR) 
ndr.hotspot <- N.fp %>% filter(site == "NDR") %>% 
  filter(fp == "1") # hypothesized hotspot sampling location 

# Carysfort Reef (CR)
cr.hotspot <- N.fp %>% filter(site == "CR") %>% 
  filter(fp == "1") # hypothesized hotspot sampling location

ndr.hotspot$var <- (ndr.hotspot$sd_relief)^2
cr.hotspot$var <- (cr.hotspot$sd_relief)^2


hotspot.pal <- pnw_palette(6, name = "Bay", type = "continuous")
# get color codes
# view(hotspot.pal)
# -----------------------------------------------------------
# NDR1 Hotspot (focal point 1)
# -----------------------------------------------------------

# Transect 1: Vertical relief by distance from hotspot 
ndr.1.relief <- ggplot(subset(ndr.hotspot, trans.cat == "1"),
       aes(x=dist.cat, y=height_diff)) +
  geom_point(color='#00496F') +
  geom_errorbar(aes(ymin=height_diff-var/2, ymax=height_diff+var/2), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~trans.cat) +
  ylim(1, 3.5) +
  ylab("Relief (m)") +
  xlab("") +
  base_theme +
  geom_line()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
# ndr.1.relief

# Transect 1: %N (macroalgal samples) by distance from hotspot 
ndr.1.percN <- ggplot(subset(ndr.hotspot, trans.cat == "1"),
                       aes(x=factor(dist.cat), 
                           y=N)) +
  geom_boxplot(fill='#00496F') +
  base_theme + 
  labs(x="",
       y="Macroalgal %N") +
  scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))

ndr.1.percN

# Transect 2: Vertical relief by distance from hotspot 
ndr.2.relief <- ggplot(subset(ndr.hotspot, trans.cat == "2"),
                       aes(x=dist.cat, y=height_diff)) +
  geom_point(color='#0C7996') +
  geom_errorbar(aes(ymin=height_diff-var/2, ymax=height_diff+var/2), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~trans.cat) +
  geom_line() +
  ylim(1, 3.5) +
  base_theme +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )  +
  ylab("Relief (m)") +
  xlab("") 
# ndr.2.relief

# Transect 2:  %N (macroalgal samples) by distance from hotspot
ndr.2.percN <- ggplot(subset(ndr.hotspot, trans.cat == "2"),
                      aes(x=factor(dist.cat), 
                          y=N)) +
  geom_boxplot(fill='#0C7996') +
  #geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
  #             position=position_dodge(.9))  +
  base_theme +
  labs(x="",
       y="Macroalgal %N") +
  scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
# ndr.2.percN


# Transect 3: Vertical relief by distance from hotspot 
ndr.3.relief <- ggplot(subset(ndr.hotspot, trans.cat == "3"),
                       aes(x=dist.cat, y=height_diff)) +
  geom_point(color='#94B669') +
  geom_errorbar(aes(ymin=height_diff-var/2, ymax=height_diff+var/2), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~trans.cat) +
  base_theme +
  geom_line() +
  ylim(1, 3.5)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  ylab("Relief (m)") +
  xlab("")
ndr.3.relief

# fix(ndr.hotspot)
# Transect 3: %N (macroalgal samples) by distance from hotspot 
ndr.3.percN <- ggplot(subset(ndr.hotspot, trans.cat == "3"),
                      aes(x=factor(dist.cat), 
                          y=N)) +
  geom_boxplot(fill='#94B669') +
  #geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
  #             position=position_dodge(.9))  +
  base_theme +
  labs(x="",
       y="Macroalgal %N") +
  scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
ndr.3.percN

# Prep estimated nutrient supply from fish df
# Note: this dataframe includes the individual survey supply rates.
supply <- supply.load

NDR.1 <- supply %>% filter(site == "NDR")

# Create sub data for plots within each transect
F1.1.df <- NDR.1 %>% filter(grid == "63" | grid == "64" | grid == "65" | grid == "71") %>% 
  dplyr::select(site, grid, visit, sum.N.m2)
F1.2.df <- NDR.1 %>% filter(grid == "63" |grid == "68" | grid == "67" | grid == "73")%>% 
  dplyr::select(site, grid, visit, sum.N.m2)
F1.3.df <- NDR.1 %>% filter(grid == "63" |grid == "57" | grid == "51" | grid == "45")%>% 
  dplyr::select(site, grid, visit, sum.N.m2)


# NDR.1b <- F1.1.df %>% rbind(F1.2.df) %>% rbind(F1.3.df)
# specify factor level order
F1.3.df$grid <- factor(F1.3.df$grid, 
                       levels=c("63", "57", "51", "45"))

# define plotting labels 
supply.labs <- c("0", "5", "10", "15")

# Transect 1: estimated nitrogen supply from fish 
ndr.1.supply <- ggplot(F1.1.df,
                     aes(x=factor(grid), 
                         y=sum.N.m2)) +
  geom_boxplot(fill='#00496F') +
  base_theme +
  labs(x="",
       y=expression(atop("N Supply", paste("(mg ",m^{-2}~day^{-1},")")))) + 
  scale_x_discrete(labels=supply.labs)
  # scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
ndr.1.supply


# Transect 2: estimated nitrogen supply from fish 
ndr.2.supply <- ggplot(F1.2.df,
                       aes(x=factor(grid), 
                           y=sum.N.m2)) +
  geom_boxplot(fill='#0C7996') +
  base_theme +
  labs(x="",
       y=expression(atop("N Supply", paste("(mg ",m^{-2}~day^{-1},")")))) +
  scale_x_discrete(labels=supply.labs)
# scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
ndr.2.supply


# Transect 3: estimated nitrogen supply from fish 
ndr.3.supply <- ggplot(F1.3.df,
                       aes(x=factor(grid), 
                           y=sum.N.m2)) +
  geom_boxplot(fill='#94B669') +
  base_theme +
  labs(x="",
       y=expression(atop("N Supply", paste("(mg ",m^{-2}~day^{-1},")")))) +
      #  y=expression(N~Supply~(mg~m^-2~day^-1))) +
  scale_x_discrete(labels=supply.labs)
# scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
ndr.3.supply

#
# Figure 6 Export: this includes ndr1.cerv and ndr3.cerv created from the script 
# 'manuscript_coral_growth_processing.R'
# png("/output/NDR_panel_plot_highres.png", width = 26, height = 14, units = 'cm', res = 300)
fig6 <- cowplot::plot_grid(ndr.1.relief, ndr.1.supply, ndr.1.percN, ndr1.cerv,
                   # ndr.2.relief, ndr.2.supply, ndr.2.percN,
                   ndr.3.relief, ndr.3.supply, ndr.3.percN, ndr3.cerv,
                   align = "v",
                   labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'),
                   label_size = 8,
                   ncol=4)
fig6

# Save the combined plot
ggsave(
  filename = "output/Ecol-Apps/Fig-6.jpg",  # Filename for the saved plot
  plot = fig6,                    # The plot to save
  width = 18,                             # Width in cm
  height = 10,                             # Height in cm (adjust as needed)
  units = "cm",                           # Units for width and height
  dpi = 500                                # Resolution in dots per inch
)

# -----------------------------------------------------------
# Appendix S1: Figure S4 - CR1 Hotspot (focal point 3)
# -----------------------------------------------------------
# CR focal points: FP 1 === F3 (grid 60) and FP 2 == F4 (grid 42)

CR.1 <- supply %>% filter(site == "CR")

# For each CR transect (A, C, and E)
F3.B.df <- CR.1 %>% filter(grid == "60" | grid == "50" | grid == "49" | grid == "48")
F3.C.df <- CR.1 %>% filter(grid == "60" |grid == "70" | grid == "80" | grid == "90")
F3.E.df <- CR.1 %>% filter(grid == "60" |grid == "69" | grid == "68" | grid == "67")

# For plotting, I need to specify the order of the grid cells to be increasing distances. Here, 
# I only had to do it for E, because the rest were in numerical order already, which is the
# standard plotting order.
F3.B.df$grid <- factor(F3.B.df$grid, levels=c("60", "49", "48", "38")) # specify factor level order
F3.C.df$grid <- factor(F3.C.df$grid, levels=c("60", "70", "80", "90"))
F3.E.df$grid <- factor(F3.E.df$grid, levels=c("60", "69", "68", "67"))


# Transect 1: estimated nutrient supply: 
cr.1.supply <- ggplot(F3.B.df,
                       aes(x=factor(grid), 
                           y=sum.N.m2)) +
  geom_boxplot(fill='#EDB829') +
  base_theme +
  labs(x="",
       y=expression(atop("N Supply", paste("(mg ",m^{-2}~day^{-1},")")))) + 
  scale_x_discrete(labels=supply.labs) 
# scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
cr.1.supply


# Transect 2: estimated nutrient supply
cr.2.supply <- ggplot(F3.C.df,
                      aes(x=factor(grid), 
                          y=sum.N.m2)) +
  geom_boxplot(fill='#E97C07') +
  base_theme +
  labs(x="",
       y=expression(atop("N Supply", paste("(mg ",m^{-2}~day^{-1},")")))) + 
  scale_x_discrete(labels=supply.labs) 
# scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
cr.2.supply


F3.E.df <- read.csv("F3.E.df.csv")
cr.3.supply <- ggplot(F3.E.df,
                      aes(x=factor(grid), 
                          y=sum.N.m2)) +
  geom_boxplot(fill='#DD4124') +
  base_theme +
  labs(x="",
       y=expression(atop("N Supply", paste("(mg ",m^{-2}~day^{-1},")")))) + 
  scale_x_discrete(labels=supply.labs) 
# scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
cr.3.supply


# Transect 1: vertical relief at increasing distance from hotspot
cr.1.relief <- ggplot(subset(cr.hotspot, trans.cat == "1"),
                      aes(x=dist.cat, y=height_diff)) +
  geom_point(color='#EDB829') +
  geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~trans.cat) +
  base_theme +
  ylim(1, 3.5)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  ylab("Relief (m)") +
  xlab("") +
  geom_line()
cr.1.relief
fix(cr.hotspot)

# Transect 1: %N (macroalgae)
cr.1.percN <- ggplot(subset(cr.hotspot, trans.cat == "1"),
                     aes(x=factor(dist.cat), 
                         y=N)) +
  geom_boxplot(fill='#EDB829') +
  #geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
  #             position=position_dodge(.9))  +
  base_theme +
  labs(x="",
       y="Macroalgal %N") +
  scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
cr.1.percN


# Transect 2: relief
cr.2.relief <- ggplot(subset(cr.hotspot, trans.cat == "2"),
                      aes(x=dist.cat, y=height_diff)) +
  geom_point(color='#E97C07') +
  geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~trans.cat) +
  base_theme +
  ylim(1, 3.5)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  ylab("Relief (m)") +
  xlab("") +
  geom_line()
cr.2.relief

# Transect 2: %N
cr.2.percN <- ggplot(subset(cr.hotspot, trans.cat == "2"),
                     aes(x=factor(dist.cat), 
                         y=N)) +
  geom_boxplot(fill='#E97C07') +
  #geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
  #             position=position_dodge(.9))  +
  base_theme +
  labs(x="",
       y="Macroalgal %N") +
  # ylim(0.5,2) 
  scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
cr.2.percN
# fix(cr.hotspot)

# Transect 3: relief
cr.3.relief <- ggplot(subset(cr.hotspot, trans.cat == "3"),
                      aes(x=dist.cat, y=height_diff)) +
  geom_point(color='#DD4124') +
  geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
                position=position_dodge(.9))  +
  facet_wrap(~trans.cat) +
  base_theme +
  ylim(1, 3.5)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  ylab("Relief (m)") +
  xlab("")  +
  geom_line()
cr.3.relief


# Transect 3: %N
cr.3.percN <- ggplot(subset(cr.hotspot, trans.cat == "3"),
                     aes(x=factor(dist.cat), 
                         y=N)) +
  geom_boxplot(fill='#DD4124') +
  #geom_errorbar(aes(ymin=height_diff-var, ymax=height_diff+var), width=.2,
  #             position=position_dodge(.9))  +
  base_theme +
  labs(x="",
       y="Macroalgal %N") +
  scale_y_continuous(breaks = seq(0.5, 2, by=0.25), limits=c(0.5,2))
cr.3.percN

# Create Figure S4
# Plot CR1 hotspot panel plot 
# png("output/CR_panel_plot_highres.png", width = 18, height = 14, units = 'cm', res = 300)
figS4 <- cowplot::plot_grid( cr.3.relief, cr.3.supply,  cr.3.percN,
                    cr.2.relief, cr.2.supply,  cr.2.percN,
                    cr.1.relief, cr.1.supply,  cr.1.percN,
                   align = "v",
                   labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'),
                   label_size = 8,
                   ncol=3)
figS4

# Save the combined plot
ggsave(
  filename = "output/Ecol-Apps/Fig-S4.png",  # Filename for the saved plot
  plot = figS4,                    # The plot to save
  width = 12,                             # Width in cm
  height = 10,                             # Height in cm (adjust as needed)
  units = "cm",                           # Units for width and height
  dpi = 300                                # Resolution in dots per inch
)
