# ========== POPULATION-LEVEL ISOTOPE NICHE ANALYSIS ==========
# Analysis of among-population variation in isotope niche with temporal and site type effects
# Uses DHGLM to model both mean isotope values and within-population variance

# Load required packages
library(brms) # Bayesian regression models using Stan
library(purrr) # Functional programming tools
library(bayestestR) # Bayesian model testing and estimation
library(ellipse) # Ellipse calculations for plotting
library(tidyverse) # Data manipulation and visualization
library(patchwork) # Combine multiple ggplot objects
library(ggridges) # Ridge density plots
source("functions.R") # Load custom helper functions


# Load and prepare population-level data
# Data contains isotope measurements across multiple populations/sites over time
DF <- read.csv("data/population_level_data.csv") %>%
  mutate(
    year = as.factor(year), # Convert year to factor
    d13C = scale(d13C), # Standardize Carbon isotope ratios
    d15N = scale(d15N), # Standardize Nitrogen isotope ratios
    date = scale(as.numeric(date)) # Standardize date (centered at mean)
  )


# Define multivariate DHGLM formulas with temporal and site type covariates
# Model both mean and variance with fixed effects of date and site type
# Random effects allow each site to have unique intercepts
form1 <- bf(d13C ~ date + site_type2 + (1 | a | site)) +
  lf(sigma ~ date + site_type2 + (1 | a | site)) # Carbon model
form2 <- bf(d15N ~ date + site_type2 + (1 | a | site)) +
  lf(sigma ~ date + site_type2 + (1 | a | site)) # Nitrogen model

# Fit the Bayesian DHGLM or load if already fitted
if (!file.exists("population_level_model.rds")) {
  model <- brm(
    formula = form1 + form2, # Joint model for both isotopes
    data = DF,
    iter = 8000, # 8000 MCMC iterations per chain
    chains = 4, # 4 independent chains
    thin = 10, # Keep every 10th sample
    cores = 4 # Parallel computation
  )
  #### Save the fitted model for future use
  saveRDS(model, file = "models/population_level_model.rds")
} else {
  # Load pre-fitted model
  model <- readRDS("models/population_level_model.rds")
}

# ========== MODEL SUMMARY AND PARAMETER EXTRACTION ==========

# Display full model summary (convergence diagnostics, parameter estimates)
summary(model)

# Extract and summarize fixed effects (date and site type effects)
fixefs <- as.data.frame(fixef(model, summary = FALSE))
summ_bayes(fixefs) # Get mean, median, and 95% HDI

# Extract and summarize random effects (among-site variance)
ranefs <- VarCorr(model, summary = FALSE)
coef_random <- as.data.frame(ranefs[["site"]][[" sd"]])^2 # Convert SD to variance
summ_bayes(coef_random)

# Function to compute intercept at mean date (date is centered at zero)
fn_int <- function(obj) {
  out <- obj[1] + obj[2] * 0 + obj[3] / 2 # Intercept + date*0 + site_type/2
  out
}

# Calculate repeatability (proportion of variance among sites)
R_pop <- repeatt(model, "site") # Custom function from functions.R
summ_bayes(R_pop)


# ========== VISUALIZATION: BLUPs AND ELLIPSES ==========

## Site-level BLUPs (Best Linear Unbiased Predictors)
# Shows site-specific niche positions and widths
p_pop <- plot_blups(model, "site") + # Custom function from functions.R
  theme(plot.title = element_text(size = 25))
p_pop
# Save high-resolution figure
ggsave(
  "population_level_blups.png",
  p_pop,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)


## Ellipse plots showing covariance structure
# Extract site-level BLUPs
blup <- ranef(model)$site[, 1, ] %>%
  data.frame() %>%
  rename_with(~ gsub("_Intercept", "", .x, fixed = TRUE)) # Clean column names

# Extract variance-covariance matrices
vars <- VarCorr(model, summary = FALSE)$site$cov %>%
  as.data.frame() %>%
  rename_with(~ gsub("_Intercept", "", .x, fixed = TRUE))


# Create 6-panel ellipse plot
pop_ellipse_plot <- plots_ell(blup, vars) # Custom function from functions.R
pop_ellipse_plot

# Save high-resolution figure
ggsave(
  "population_level_ellipses.png",
  pop_ellipse_plot,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)


# ========== SITE-SPECIFIC WITHIN-POPULATION VARIANCE RIDGE PLOT ==========

# Extract posterior draws for all parameters
comunidad1 <- as_draws_df(model)

# Select only sigma random effects (within-population variance components)
df_filtrado <- comunidad1 %>%
  select(starts_with("r_site__sigma_d"))

# Back-transform rIIV from log scale to original variance scale
# Reference: Hertel & Niemelä (2020) - guide for studying among-individual behavioral variation
community.sp <- exp(df_filtrado)
Sp.C <- community.sp[, c(1:17)] # Sites 1-17 for Carbon
Sp.N <- community.sp[, c(18:34)] # Sites 1-17 for Nitrogen

# Reshape Carbon data from wide to long format
df_Carbon <- Sp.C %>%
  pivot_longer(
    cols = everything(),
    names_to = "Population", # Site/population identifier
    values_to = "Intraindividual Variance" # Within-population variance
  ) %>%
  mutate(
    Population = gsub("r_site__sigma_d13C\\[|,Intercept\\]", "", Population)
  ) %>%
  arrange(Population) # Sort by population
df_Carbon <- as.data.frame(df_Carbon)
df_Carbon$Population <- as.factor(df_Carbon$Population)
str(df_Carbon)

# Reshape Nitrogen data from wide to long format
df_Nitrogen <- Sp.N %>%
  pivot_longer(
    cols = everything(),
    names_to = "Population",
    values_to = "Intraindividual Variance"
  ) %>%
  mutate(
    Population = gsub("r_site__sigma_d15N\\[|,Intercept\\]", "", Population)
  ) %>%
  arrange(Population) # Sort by population

df_N <- as.data.frame(df_Nitrogen)
df_N$Population <- as.factor(df_N$Population)
str(df_N)

# Combine Carbon and Nitrogen data
alls <- rbind(df_N, df_Carbon)
Isotope <- factor(c(rep("Nitrogen", 27200), rep("Carbon", 27200))) # 27200 posterior samples per isotope

# Create dataframe with isotope labels
df_largo4 <- data.frame(Isotope = Isotope)
DF_plot <- cbind(alls, df_largo4)

# Order populations by mean within-population variance (low to high)
DF_plot$Population <- factor(
  DF_plot$Population,
  levels = DF_plot %>%
    group_by(Population) %>%
    summarise(Mean_IV = mean(`Intraindividual Variance`)) %>%
    arrange(Mean_IV) %>%
    pull(Population)
)

# Create ridge density plot showing within-population variance for each site
pop_plot <- ggplot(
  DF_plot,
  aes(x = `Intraindividual Variance`, y = Population, fill = Isotope)
) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) + # Overlapping density ridges
  theme_classic() +
  scale_fill_manual(values = c("#FF6347", "#4682B4")) + # Red for Carbon, blue for Nitrogen
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "Intraindividual Variance", y = "Populations", title = "") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) +
  theme(legend.position = "none")

pop_plot
# Save high-resolution ridge plot
ggsave(
  "population_level_residual_intraunit_variance.png",
  pop_plot,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)
