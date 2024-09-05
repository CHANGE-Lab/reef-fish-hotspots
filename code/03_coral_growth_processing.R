# -----------------------------------------------------------
# Quantifying fish-derived nutrient hotspots across reefscapes
# By Noelle K. Helder
# Updated March 31 2023

# Manuscript DOI: TBD 
# -----------------------------------------------------------
# The purpose of this script is to process coral monitoring data 
# from the Coral Restoration Foundation, and produce coral growth 
# plots used in Figure 6 of the associated manuscript. 
# note that the words plot & grid are used interchangeably

# Inputs
#   @coral: "CRF_coral_data.csv" - monitoring data from CRF
#   @grid_cowtag_table: "grid_cowtag_table.xlsx" - dataframe to link coral 
#                        ID (cow tag number) to corresponding survey plots (25m^2)

# Outputs
#   ndr1.cerv: ggplot of coral growth rates per plot at increasing distances from NDR hotspot, transect 1
#   ndr3.cerv: ggplot of coral growth rates per plot at increasing distances from NDR hotspot, transect 3

# Use
#   Outputs used in script '03_hotspots_panel_plots.R' to generate manuscript Figure 6
# -----------------------------------------------------------

# Load packages
library(readxl)
library(tidyverse)

# load data that matches the plots with the coral ID (cowtag number)
grid_cowtag_table <- read_excel("data/grid_cowtag_table.xlsx")
# Load CRF data, calculate growth (year measurement - month measurement)/365
coral <- read_csv("data/CRF_coral_data.csv") 
names(coral)

# Tidy
coral <- coral %>% dplyr::rename_at('Monitoring Date', ~'Date')%>% 
                          rename_at('Cow Tag', ~'Cow_tag') %>% 
                          rename_at('Coral Letter', ~'Letter')%>% 
                          rename_at('Max Diameter (cm)', ~'Diameter')

coral$Cow_tag <- as.factor(coral$Cow_tag)
coral$Letter <- as.factor(coral$Letter)
coral$Diameter <- as.numeric(coral$Diameter)

coral$Date= as.Date(coral$Date,format= "%m/%d/%Y")

# confirm that we have 1 monitoring date in June 2019 (1 month) and March 2020 (1 year) 
# to calculate growth rate
unique(coral$Date)


# Calculate coral growth rates (1 year max diameter - 1 month max diameter)/365
cGrowth <- coral %>% dplyr::select(c('Date', 'Cow_tag', 'Letter', 'Diameter')) %>% 
  group_by(Cow_tag, Letter, Date) %>% 
  arrange(Cow_tag, Letter, Date) %>% 
  # remove survey duplicates  
  distinct(.keep_all = TRUE) %>% 
  ungroup() %>% 
  group_by(Cow_tag, Letter) %>% 
  # Difference in growth between monitoring periods
  mutate(Diff_growth = Diameter - lag(Diameter)) %>%  
  drop_na() %>% 
  dplyr::select(c('Cow_tag', 'Letter', 'Diff_growth')) # growth diff in cm for each cow-tag letter coral 

# hist(cGrowth$Diff_growth)


# coral_temp <-  coral %>% select(c('Date', 'Cow_tag', 'Letter', 'Diameter')) %>% 
#   group_by(Cow_tag, Letter, Date) %>% 
#  arrange(Cow_tag, Letter, Date)

# calculate growth rate, scale to mm per day
cGrowth$Growth_rate = (cGrowth$Diff_growth)/365*10
# hist(cGrowth$Growth_rate)

# Combing growth data with survey plot numbers
# set data types for join
grid_cowtag_table$Grid <- as.factor(grid_cowtag_table$Grid)
grid_cowtag_table$Cow_tag <- as.factor(grid_cowtag_table$Cow_tag)
coral2$Cow_tag <- as.factor(coral2$Cow_tag)

# Join
cGrowth_rate <- cGrowth %>% left_join(grid_cowtag_table) %>% na.omit()
 
# looking at initial distributions and relationships with quadrats. 
# ggplot(cGrowth_rate, aes(x=factor(Grid), y=Diff_growth)) +
#   geom_boxplot(fill='#EDB829') +
#   theme_classic() +
#   labs(x="",
#        y="Diff in growth By Grid") #+ 
  # scale_x_discrete(labels=supply.labs) 

# ggplot(cGrowth_rate, aes(x=factor(Grid), y=Growth_rate)) +
#   geom_boxplot(fill='#EDB829') +
#   theme_classic() +
#   labs(x="",
#        y="Growth Rates By Grid") #+ 
# scale_x_discrete(labels=supply.labs) 


### NDR Transects CRF growth for A. cerv: exploring growth rates as a function of distance from hotspot
# 0 m = hotspot (grid == 63); other grids at 5, 10, 15 m, respectively. 
# Only have coral data for 2/3 transects.

F1.1.df <- cGrowth_rate %>% filter(Grid == "63" | Grid == "64" | Grid == "65" |Grid == "66") 
# F1.2.df <- test %>% filter(Grid == "63" |Grid == "68" | Grid == "67" | Grid == "73") 
F1.3.df <- cGrowth_rate %>% filter(Grid == "63" |Grid == "57" | Grid == "51" | Grid == "58")

# specify factor level order
F1.3.df$Grid <- factor(F1.3.df$Grid, 
                       levels=c("63", "57", "51", "58")) 
# create figure labels 
supply.labs <- c("0", "5", "10", "15")

# plot

# Define a base theme for consistency
base_theme <- theme_classic() +
  theme(
    text = element_text(size = 8),   # Adjust text size to 8 points
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 8)
  )

ndr1.cerv <- ggplot(F1.1.df,
                       aes(x=factor(Grid), 
                           y=Growth_rate)) +
              geom_boxplot(fill='#00496F') +
              base_theme +
  # theme(panel.background = element_blank(),
  #       axis.line = element_line(colour = 'black', size = 1.2),
  #       axis.text=element_text(size=18),
  #       axis.title=element_text(size=20)) +
              labs(x="",
                   y= expression(atop("Coral Growth Rate", paste("(mm "~day^{-1},")")))) + 
              scale_x_discrete(labels=supply.labs) 
ndr1.cerv



ndr3.cerv <- ggplot(F1.3.df,
                       aes(x=factor(Grid), 
                           y=Growth_rate)) +  
             geom_boxplot(fill='#94B669') +
             base_theme +
  # theme(panel.background = element_blank(),
  #       axis.line = element_line(colour = 'black', size = 1.2),
  #       axis.text=element_text(size=18),
  #       axis.title=element_text(size=20))  +
  # theme_classic() +
             labs(x="",
                  y= expression(atop("Coral Growth Rate", paste("(cm "~day^{-1},")")))) + 
             scale_x_discrete(labels=supply.labs)
ndr3.cerv

