# ========== COMMUNITY-LEVEL ISOTOPE NICHE ANALYSIS ==========
# Analysis of among-species variation in isotope niche position and width
# Uses DHGLM to model both mean isotope values and within-species variance

# Load required packages
library(brms) # Bayesian regression models using Stan
library(purrr) # Functional programming tools
library(bayestestR) # Bayesian model testing and estimation
library(ellipse) # Ellipse calculations for plotting
library(tidyverse) # Data manipulation and visualization
library(patchwork) # Combine multiple ggplot objects
library(ggridges) # Ridge density plots
source("functions.R") # Load custom helper functions


# Load and standardize community-level data
# Data contains isotope measurements (d13C, d15N) for multiple species
DF <- read.csv("data/community_level_data.csv") %>%
  mutate(
    mass_sc = scale(drymass), # Standardize body mass
    d13C = scale(d13c), # Standardize Carbon isotope ratios
    d15N = scale(d15n) # Standardize Nitrogen isotope ratios
  )

# Define multivariate DHGLM formulas
# Jointly model d13C and d15N, with correlation between isotopes
# Model both mean (location) and variance (dispersion) with species-level random effects
form1 <- bf(d13C ~ (1 | a | species)) + lf(sigma ~ (1 | a | species)) # Carbon model
form2 <- bf(d15N ~ (1 | a | species)) + lf(sigma ~ (1 | a | species)) # Nitrogen model

# Fit the Bayesian DHGLM or load if already fitted
if (!file.exists("models/community_level_model.rds")) {
  model <- brm(
    formula = form1 + form2 + set_rescor(TRUE), # Enable residual correlation between isotopes
    data = DF,
    iter = 4000, # 4000 MCMC iterations per chain
    chains = 4, # 4 independent chains
    thin = 10, # Keep every 10th sample
    cores = 4 # Parallel computation
  )
  #### Save the fitted model for future use
  saveRDS(model, file = "models/community_level_model.rds")
} else {
  # Load pre-fitted model
  model <- readRDS("models/community_level_model.rds")
}

# ========== MODEL SUMMARY AND PARAMETER EXTRACTION ==========

# Display full model summary (convergence diagnostics, parameter estimates)
summary(model)

# Extract and summarize fixed effects (population-level intercepts)
fixefs <- as.data.frame(fixef(model, summary = FALSE))
summ_bayes(fixefs) # Get mean, median, and 95% HDI

# Extract and summarize random effects (among-species variance)
ranefs <- VarCorr(model, summary = FALSE)
coef_random <- as.data.frame(ranefs[["species"]][["sd"]])^2 # Convert SD to variance
summ_bayes(coef_random)

# Function to compute intercept (no covariates at community level)
fn_int <- function(obj) {
  out <- obj[1] # Just the intercept
  out
}

# Calculate repeatability (proportion of variance among species)
R_com <- repeatt(model, "species") # Custom function from functions.R
summ_bayes(R_com)


# ========== VISUALIZATION: BLUPs AND ELLIPSES ==========

## Species-level BLUPs (Best Linear Unbiased Predictors)
# Shows species-specific niche positions and widths
p_com <- plot_blups(model, "species") + # Custom function from functions.R
  theme(plot.title = element_text(size = 25))
p_com

if (!file.exists("figs/community_level_blups.png")) {
  # Save high-resolution figure
  ggsave(
    "figs/community_level_blups.png",
    p_com,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300
  )
}

## Ellipse plots showing covariance structure
# Extract species-level BLUPs
blup <- ranef(model)$species[, 1, ] %>%
  data.frame() %>%
  rename_with(~ gsub("_Intercept", "", .x, fixed = TRUE)) # Clean column names

# Extract variance-covariance matrices
vars <- VarCorr(model, summary = FALSE)$species$cov %>%
  as.data.frame() %>%
  rename_with(~ gsub("_Intercept", "", .x, fixed = TRUE))

# Create 6-panel ellipse plot
com_ellipse_plot <- plots_ell(blup, vars) # Custom function from functions.R

com_ellipse_plot

if (!file.exists("figs/community_level_ellipses.png")) {
  # Save high-resolution figure
  ggsave(
    "figs/community_level_ellipses.png",
    com_ellipse_plot,
    width = 20,
    height = 15,
    units = "in",
    dpi = 300
  )
}

# ========== SPECIES-SPECIFIC WITHIN-SPECIES VARIANCE (rIIV) RIDGE PLOT ==========

# Extract posterior draws for all parameters
comunidad1 <- as_draws_df(model)

# Select only sigma random effects (within-species variance components)
df_filtrado <- comunidad1 %>%
  select(starts_with("r_species__sigma_d"))

# Back-transform rIIV from log scale to original variance scale
# Reference: Hertel & Niemelä (2020) - guide for studying among-individual behavioral variation
community.sp <- exp(df_filtrado)
Sp.C <- community.sp[, c(1:10)] # Species 1-10 for Carbon
Sp.N <- community.sp[, c(11:20)] # Species 1-10 for Nitrogen

# Reshape Carbon data from wide to long format
df_Carbon <- Sp.C %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sp", # Species identifier
    values_to = "Intraindividual Variance" # Within-species variance
  ) %>%
  mutate(Sp = gsub("r_species__sigma_d13c\\[|,Intercept\\]", "", Sp)) %>%
  arrange(Sp) # Sort by species

df_Carbon <- as.data.frame(df_Carbon)
df_Carbon$Sp <- as.factor(df_Carbon$Sp)
str(df_Carbon)

# Reshape Nitrogen data from wide to long format
df_Nitrogen <- Sp.N %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sp",
    values_to = "Intraindividual Variance"
  ) %>%
  mutate(Sp = gsub("r_species__sigma_d15n\\[|,Intercept\\]", "", Sp)) %>%
  arrange(Sp) # Sort by species

df_N <- as.data.frame(df_Nitrogen)
df_N$Sp <- as.factor(df_N$Sp)

# Combine Carbon and Nitrogen data
alls <- rbind(df_N, df_Carbon)
Isotope <- factor(c(rep("Nitrogen", 8000), rep("Carbon", 8000))) # 8000 posterior samples per isotope

# Create dataframe with isotope labels
df_largo4 <- data.frame(Isotope = Isotope)
DF_plot <- cbind(alls, df_largo4)
str(DF_plot)

# Order species by mean within-species variance (low to high)
DF_plot$Sp <- factor(
  DF_plot$Sp,
  levels = DF_plot %>%
    group_by(Sp) %>%
    summarise(Mean_IV = mean(`Intraindividual Variance`)) %>%
    arrange(Mean_IV) %>%
    pull(Sp)
)

# Create ridge density plot showing within-species variance for each species
spp_plot <- ggplot(
  DF_plot,
  aes(x = `Intraindividual Variance`, y = Sp, fill = Isotope)
) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) + # Overlapping density ridges
  theme_classic() +
  scale_fill_manual(values = c("#FF6347", "#4682B4")) + # Red for Carbon, blue for Nitrogen
  scale_x_continuous(limits = c(0, 3)) +
  labs(x = "Intraindividual Variance", y = "Species", title = "") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) +
  theme(legend.position = "none")

spp_plot

# Save high-resolution ridge plot

if (!file.exists("figs/community_level_intraspp_variance.png")) {
  ggsave(
    "figs/community_level_intraspp_variance.png",
    spp_plot,
    width = 12,
    height = 15,
    units = "in",
    dpi = 300
  )

}