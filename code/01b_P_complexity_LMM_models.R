# This script identifies the top model predicting P supply from fishes
#   as a function of reef complexity (relief, vrm (1 cm), vrm deviation)

# Output: Appendix S1:Table S3

# Load necessary libraries
library(lmer)
library(lmerTest) # For lmer models
library(broom.mixed) # For tidy() function on mixed models
library(performance) # For AIC
library(kableExtra) # For publication-ready tables
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
# --------------------------------------


# Define all models in a list with updated numbering
models <- list(
  # Full model with interactions (excluding vrm_1cm_sc*vrm_dev_sc)
  P.model1 = lmer(log(mean.P) ~ poly(relief, 2) * poly(vrm_1cm_sc, 2) + poly(relief, 2) * poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  
  # Two-way interactions with quadratic terms (excluding vrm_1cm_sc*vrm_dev_sc)
  P.model2 = lmer(log(mean.P) ~ poly(relief, 2) * poly(vrm_1cm_sc, 2) + poly(relief, 2) * vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model3 = lmer(log(mean.P) ~ poly(relief, 2) * vrm_1cm_sc + poly(relief, 2) * poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model4 = lmer(log(mean.P) ~ poly(relief, 2) * vrm_1cm_sc + poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model5 = lmer(log(mean.P) ~ poly(relief, 2) * vrm_dev_sc + poly(vrm_1cm_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model6 = lmer(log(mean.P) ~ poly(relief, 2) * vrm_dev_sc + vrm_1cm_sc + (1 | site), REML = FALSE, data = df),
  P.model7 = lmer(log(mean.P) ~ poly(relief, 2) + poly(vrm_1cm_sc, 2) + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model8 = lmer(log(mean.P) ~ poly(relief, 2) + poly(vrm_1cm_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model9 = lmer(log(mean.P) ~ poly(relief, 2) + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model10 = lmer(log(mean.P) ~ poly(relief, 2) + vrm_1cm_sc + (1 | site), REML = FALSE, data = df),
  P.model11 = lmer(log(mean.P) ~ poly(vrm_1cm_sc, 2) + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model12 = lmer(log(mean.P) ~ poly(vrm_1cm_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model13 = lmer(log(mean.P) ~ poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  
  # Linear terms with interactions
  P.model14 = lmer(log(mean.P) ~ relief * vrm_1cm_sc + relief * vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model15 = lmer(log(mean.P) ~ relief * vrm_1cm_sc + (1 | site), REML = FALSE, data = df),
  P.model16 = lmer(log(mean.P) ~ relief * vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model17 = lmer(log(mean.P) ~ relief + vrm_1cm_sc + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model18 = lmer(log(mean.P) ~ relief + vrm_1cm_sc + (1 | site), REML = FALSE, data = df),
  P.model19 = lmer(log(mean.P) ~ relief + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  
  # Quadratic terms only, no interactions
  P.model20 = lmer(log(mean.P) ~ poly(relief, 2) + poly(vrm_1cm_sc, 2) + poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model21 = lmer(log(mean.P) ~ poly(relief, 2) + poly(vrm_1cm_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model22 = lmer(log(mean.P) ~ poly(relief, 2) + poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model23 = lmer(log(mean.P) ~ poly(vrm_1cm_sc, 2) + poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model24 = lmer(log(mean.P) ~ poly(vrm_1cm_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model25 = lmer(log(mean.P) ~ poly(vrm_dev_sc, 2) + (1 | site), REML = FALSE, data = df),
  P.model26 = lmer(log(mean.P) ~ poly(relief, 2) + (1 | site), REML = FALSE, data = df),
  
  # Linear terms only, no interactions
  P.model27 = lmer(log(mean.P) ~ relief + vrm_1cm_sc + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model28 = lmer(log(mean.P) ~ relief + vrm_1cm_sc + (1 | site), REML = FALSE, data = df),
  P.model29 = lmer(log(mean.P) ~ relief + vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model30 = lmer(log(mean.P) ~ vrm_1cm_sc + (1 | site), REML = FALSE, data = df),
  P.model31 = lmer(log(mean.P) ~ vrm_dev_sc + (1 | site), REML = FALSE, data = df),
  P.model32 = lmer(log(mean.P) ~ relief + (1 | site), REML = FALSE, data = df)
)

# Compute AIC for each model
aic_values <- sapply(models, AIC)

# Combine model names with their AIC values
model_aic <- data.frame(
  Model = names(aic_values),
  AIC = aic_values
)

# Sort models by AIC
model_aic <- model_aic %>% arrange(AIC)

# Print the model comparison table
print(model_aic)


# Assuming models are stored in 'models'
# Extract coefficients and CIs
results <- lapply(models, function(model) {
  tidy_model <- tidy(model, conf.int = TRUE)  # Extract coefficients and CI
  aic_score <- AIC(model)  # Calculate AIC
  tidy_model$aic <- aic_score
  return(tidy_model)
})

# Combine results into one data frame
results_df <- do.call(rbind, results)

# Add model names as a column
results_df$model <- rep(names(models), sapply(results, nrow))

# Select only necessary columns
results_df <- results_df[, c("model", "term", "estimate", "conf.low", "conf.high", "aic")]

# Format and create publication-ready table
results_table <- results_df %>%
  kable(format = "html", escape = FALSE, 
        col.names = c("Model", "Term", "Estimate", "CI Lower", "CI Upper", "AIC"),
        caption = "Model Coefficients, Confidence Intervals, and AIC Scores") %>%
  kable_styling(full_width = F, position = "left", bootstrap_options = c("striped", "hover"))

# Save the table to an HTML file
writeLines(as.character(results_table), "output/Ecol-Apps/P_structure_model_comparison_table.doc")

# Alternatively, display the table in RStudio viewer
results_table
