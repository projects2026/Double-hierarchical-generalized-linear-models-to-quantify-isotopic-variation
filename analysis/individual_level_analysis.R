# Load required packages for Bayesian modeling, data manipulation, and visualization
library(brms) # Bayesian regression models using Stan
library(purrr) # Functional programming tools
library(bayestestR) # Bayesian model testing and estimation
library(ellipse) # Ellipse calculations for plotting
library(tidyverse) # Data manipulation and visualization (dplyr, ggplot2, tidyr, etc.)
library(patchwork) # Combine multiple ggplot objects
library(ggridges) # Ridge density plots
source("functions.R") # Load custom helper functions (plot_blups, repeatt, summ_bayes, plots_ell)


# Load individual-level isotope data (d13C and d15N measurements with individual ID, sex, and foraging strategy)
DF <- read.csv("data/individual_level_data.csv")

# Define DHGLM (Double Hierarchical Generalized Linear Model) formulas
# Model both the mean (location) and variance (dispersion) of d13C
# form1: Mean structure - d13C varies by sex and foraging strategy, with random intercepts by individual
# form2: Variance structure - sigma (within-individual variance) also varies by sex, foraging, and individual
form1 <- bf(d13C ~ Sexo + foraging + (1 | a | ID)) +
  lf(sigma ~ Sexo + foraging + (1 | a | ID))
form2 <- bf(d13C ~ Sexo + foraging + (1 | a | ID)) +
  lf(sigma ~ Sexo + foraging + (1 | a | ID))

# Fit the Bayesian DHGLM or load it if already fitted
# This avoids re-running the computationally expensive model fitting
if (!file.exists("individual_level_model.rds")) {
  model <- brm(
    formula = form1 + form2, # Combined mean + variance structure
    data = DF,
    iter = 8000, # 8000 MCMC iterations per chain
    chains = 4, # 4 independent MCMC chains for convergence diagnostics
    thin = 10, # Keep every 10th sample (reduces autocorrelation)
    cores = 4 # Parallel computation across 4 cores
  )
  #### Save the fitted model to disk for future use
  saveRDS(model, file = "models/individual_level_model.rds")
} else {
  # Load pre-fitted model to save time
  model <- readRDS("models/individual_level_model.rds")
}

# ========== MODEL SUMMARY AND PARAMETER EXTRACTION ==========

# Display full model summary (fixed effects, random effects, convergence diagnostics)
summary(model)

# Extract and summarize fixed effects (population-level effects of sex and foraging)
# Returns posterior samples for intercepts and coefficients
fixefs <- as.data.frame(fixef(model, summary = FALSE))
summ_bayes(fixefs) # Custom function to compute Bayesian credible intervals and means

# Extract and summarize random effects (among-individual variance components)
# VarCorr extracts variance-covariance matrices for random effects
ranefs <- VarCorr(model, summary = FALSE)
coef_random <- as.data.frame(ranefs[["ID"]][["sd"]])^2 # Convert SD to variance
summ_bayes(coef_random)

# Custom function to compute average effect (including intercept and average fixed effects)
fn_int <- function(obj) {
  out <- obj[1] + obj[2] / 2 + obj[3] / 2 + obj[4] / 4
  out
}

# Calculate repeatability (proportion of total variance explained by among-individual differences)
# High repeatability = consistent individual differences in trait values
RAA <- repeatt(model, "ID") # Custom function from functions.R
summ_bayes(RAA)


# ========== VISUALIZATION: BLUPs (Best Linear Unbiased Predictors) ==========

# Plot individual-level random effects (BLUPs) showing how each individual deviates from population mean
# Creates scatter plot of individual estimates for mean and variance components
p_ind <- plot_blups(model, "ID") + # Custom function from functions.R
  ggtitle("A. australis") +
  theme(plot.title = element_text(size = 25))
p_ind
# Save high-resolution figure
ggsave(
  "individual_level_blups.png",
  p_ind,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)


# ========== VISUALIZATION: ELLIPSES (Bivariate Individual Variation) ==========

# Extract BLUPs (random effects) for each individual
# Shows individual deviations in both mean isotope values and within-individual variance
blup <- ranef(model)$ID[, 1, ] %>%
  data.frame() %>%
  rename_with(~ gsub("_Intercept", "", .x, fixed = TRUE)) # Clean column names

# Extract variance-covariance matrices for random effects
# Used to draw confidence ellipses around individual estimates
vars <- VarCorr(model, summary = FALSE)$ID$cov %>%
  as.data.frame() %>%
  rename_with(~ gsub("_Intercept", "", .x, fixed = TRUE))

# Create ellipse plot showing bivariate distribution of individual effects
ind_ellipse_plot <- plots_ell(blup, vars) # Custom function from functions.R
ind_ellipse_plot
# Save high-resolution figure
ggsave(
  "individual_level_ellipses.png",
  ind_ellipse_plot,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300
)


# ========== EXTRACT INDIVIDUAL-LEVEL INTRA-INDIVIDUAL VARIANCE (rIIV) ==========

# Convert model output to draws dataframe containing all posterior samples
comunidad1 <- as_draws_df(model)

# Select only columns for sigma (within-individual variance) random effects
# These represent individual-specific variance components for d13C and d15N
df_filtrado <- comunidad1 %>%
  select(starts_with("r_ID__sigma_d"))

# Back-transform rIIV from log scale to original scale
# In DHGLMs, sigma is modeled on log scale, so exp() returns to sd scale
community.sp <- exp(df_filtrado)

# Split by isotope: columns 1-37 are Carbon (d13C), columns 38-74 are Nitrogen (d15N)
Sp.C <- community.sp[, c(1:37)] # Individual-specific variance for Carbon
Sp.N <- community.sp[, c(38:74)] # Individual-specific variance for Nitrogen
str(Sp.C)
str(Sp.N)

# ========== RESHAPE DATA FOR VISUALIZATION ==========

# Reshape Carbon data from wide to long format (one row per posterior draw per individual)
df_Carbon <- Sp.C %>%
  pivot_longer(
    cols = everything(),
    names_to = "ID", # Individual ID
    values_to = "Intraindividual Variance" # Posterior samples of within-individual variance
  ) %>%
  mutate(ID = gsub("r_ID__sigma_d13C\\[|,Intercept\\]", "", ID)) %>% # Clean ID names
  arrange(ID) # Sort by individual ID

df_Carbon <- as.data.frame(df_Carbon)
df_Carbon$ID <- as.factor(df_Carbon$ID)

# Reshape Nitrogen data from wide to long format
df_Nitrogen <- Sp.N %>%
  pivot_longer(
    cols = everything(),
    names_to = "ID",
    values_to = "Intraindividual Variance"
  ) %>%
  mutate(ID = gsub("r_ID__sigma_d15N\\[|,Intercept\\]", "", ID)) %>% # Clean ID names
  arrange(ID)

df_N <- as.data.frame(df_Nitrogen)
df_N$ID <- as.factor(df_N$ID)

# Combine Carbon and Nitrogen data into single dataframe
alls <- rbind(df_N, df_Carbon)
nrow(df_N) # Should be 296,000 rows (8000 posterior samples × 37 individuals)
nrow(df_Carbon) # Should be 296,000 rows

# Create isotope type identifier (Nitrogen first, then Carbon)
Isotope <- factor(c(rep("Nitrogen", 296000), rep("Carbon", 296000)))

# Create dataframe with isotope labels
df_largo4 <- data.frame(Isotope = Isotope)

# Combine variance estimates with isotope labels
DF_plot <- cbind(alls, df_largo4)
str(DF_plot)

# ========== SEX-SPECIFIC RIDGE PLOTS: FEMALES ==========

# Define list of female individual IDs
hembras <- c(
  "167",
  "168",
  "169",
  "170",
  "172",
  "Aa156",
  "Aa159",
  "Aa181",
  "Aa182",
  "Aa184",
  "Aa185",
  "Aa186",
  "Aa187",
  "Aa188",
  "Aa189",
  "Aa190",
  "Aah10",
  "Aah11",
  "Aah13",
  "Aah14",
  "Aah9"
)

# Subset data to include only female individuals
hembras1 <- subset(DF_plot, ID %in% hembras)
head(hembras1)

# Calculate mean intra-individual variance for each female
# Used for ordering individuals by their average variance
DF_plot_fem <- hembras1 %>%
  group_by(ID) %>%
  mutate(Mean_IV = mean(`Intraindividual Variance`)) %>%
  ungroup()

# Order individuals by their mean intra-individual variance (low to high)
# This makes the ridge plot easier to interpret
DF_plot_fem$ID <- factor(
  DF_plot_fem$ID,
  levels = DF_plot_fem %>%
    group_by(ID) %>%
    summarise(Mean_IV = mean(`Intraindividual Variance`)) %>%
    arrange(Mean_IV) %>%
    pull(ID)
)


# Create ridge density plot for females
# Shows full posterior distribution of intra-individual variance for each female
# Separate distributions for Carbon (red) and Nitrogen (blue) isotopes
fem_plot <- ggplot(
  DF_plot_fem,
  aes(x = `Intraindividual Variance`, y = ID, fill = Isotope)
) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) + # Overlapping density ridges
  theme_classic() +
  scale_fill_manual(values = c("#FF6347", "#4682B4")) + # Tomato red for Carbon, steel blue for Nitrogen
  scale_x_continuous(limits = c(0, 4)) +
  labs(x = "Intraindividual Variance", y = "Individuals (ID)", title = "") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) +
  theme(legend.position = "none") # Remove legend (will be in combined plot)

fem_plot

# ========== SEX-SPECIFIC RIDGE PLOTS: MALES ==========

# Define list of male individual IDs
machos <- c(
  "Aam1",
  "Aam10",
  "Aam11",
  "Aam12",
  "Aam17",
  "Aam2",
  "Aam3",
  "Aam4",
  "Aam5",
  "Aam6",
  "Aam7",
  "Aam8",
  "Aam9",
  "190509",
  "270709",
  "40209"
)

# Subset data to include only male individuals
machos1 <- subset(DF_plot, ID %in% machos)
head(machos1)

# Calculate mean intra-individual variance for each male
DF_plot_male <- machos1 %>%
  group_by(ID) %>%
  mutate(Mean_IV = mean(`Intraindividual Variance`)) %>%
  ungroup()

# Order individuals by their mean intra-individual variance (low to high)
DF_plot_male$ID <- factor(
  DF_plot_male$ID,
  levels = DF_plot_male %>%
    group_by(ID) %>%
    summarise(Mean_IV = mean(`Intraindividual Variance`)) %>%
    arrange(Mean_IV) %>%
    pull(ID)
)

# Create ridge density plot for males (matching style of female plot)
males_plot <- ggplot(
  DF_plot_male,
  aes(x = `Intraindividual Variance`, y = ID, fill = Isotope)
) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) +
  theme_classic() +
  scale_fill_manual(values = c("#FF6347", "#4682B4")) + # Same colors as female plot
  scale_x_continuous(limits = c(0, 4)) +
  labs(x = "Intraindividual Variance", y = "Individuals (males)", title = "") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) +
  theme(axis.title.y = element_blank()) # Remove y-axis label for cleaner combined plot
males_plot

# ========== COMBINE AND SAVE FINAL FIGURE ==========

# Combine female and male ridge plots side-by-side using patchwork
df_plot <- fem_plot + males_plot
df_plot

# Save combined figure as high-resolution TIFF for publication
ggsave(
  "individual level.tiff",
  df_plot,
  units = "px",
  width = 8000, # 8000 pixels width
  height = 5000, # 5000 pixels height
  device = tiff,
  dpi = 350 # 350 dots per inch for publication quality
)
