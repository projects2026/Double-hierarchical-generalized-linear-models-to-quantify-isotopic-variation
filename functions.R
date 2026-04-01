# ========== BAYESIAN SUMMARY FUNCTION ==========

# Compute summary statistics for Bayesian posterior distributions
# Returns mean, median, and 95% highest density interval (HDI) for each parameter
# HDI = shortest interval containing 95% of posterior probability
summ_bayes <- function(object) {
  Reduce(
    merge,
    list(
      point_estimate(object, "mean"), # Posterior mean
      point_estimate(object, "median"), # Posterior median
      hdi(object, ci = 0.95) # 95% credible interval (HDI)
    )
  )
}

###############################################################################

# ========== REPEATABILITY CALCULATION ==========

# Calculate repeatability (R) for each isotope in a DHGLM
# Repeatability = proportion of total variance due to among-unit differences
# R = Var_among / (Var_among + Var_within)
# High R indicates consistent individual/species/site differences
#
# Arguments:
#   model: brms DHGLM model object
#   r_e: random effect grouping variable name (e.g., "ID", "species", "site")
#
# Returns:
#   Data frame with repeatability estimates (r2_c for Carbon, r2_n for Nitrogen)
repeatt <- function(model, r_e) {
  # Extract among-unit variance (random effect SD squared) for each parameter
  # Converts SD to variance and renames "Intercept" to "var" for clarity
  vars <- VarCorr(model, summary = FALSE)[[r_e]]$sd^2 %>%
    as.data.frame() %>%
    rename_with(~ gsub("Intercept", "var", .x, fixed = TRUE))

  # Extract fixed effects for sigma (within-unit variance parameters)
  vr_exp <- fixef(model, summary = FALSE) %>%
    data.frame() %>%
    select(starts_with("sigma"))

  # Get sigma estimates for Carbon (d13C)
  sig_d13C <- select(vr_exp, starts_with("sigma_d13C")) %>%
    fn_int() %>% # Apply intercept function (defined in analysis scripts)
    rename(sig_d13C = "sigma_d13C_Intercept")

  # Get sigma estimates for Nitrogen (d15N)
  sig_d15N <- select(vr_exp, starts_with("sigma_d15N")) %>%
    fn_int() %>%
    rename(sig_d15N = "sigma_d15N_Intercept")

  # Calculate repeatability for both isotopes
  var_pd <- data.frame(vars) %>%
    cbind(sig_d13C) %>%
    cbind(sig_d15N) %>%
    mutate(
      # Back-transform within-unit variance from log scale
      # exp(log_sigma + var_log_sigma/2)^2 accounts for lognormal distribution
      v_r_c = exp(sig_d13C + sigma_d13C_var / 2)^2,
      v_r_n = exp(sig_d15N + sigma_d15N_var / 2)^2,
      # Repeatability = among-unit variance / total variance
      r2_c = d13C_var / (d13C_var + v_r_c), # Carbon repeatability
      r2_n = d15N_var / (d15N_var + v_r_n) # Nitrogen repeatability
    )
}


# ========== BLUP PLOTTING FUNCTION ==========

# Create bivariate plot of BLUPs (Best Linear Unbiased Predictors) for mean isotope values
# Shows each unit (individual/species/site) in isotope bi-plot space
# with error bars representing within-unit variance (sigma)
#
# Arguments:
#   model: brms DHGLM model object
#   r_e: random effect grouping variable name
#
# Returns:
#   ggplot object showing d13C vs d15N with variance error bars
plot_blups <- function(model, r_e) {
  # Extract population-level (fixed effect) estimates
  fix_coef <- fixef(model, summary = TRUE)[, 1]

  # Get population mean for Carbon (d13C)
  int_d13C <- select(data.frame(t(fix_coef)), starts_with("d13C")) %>%
    fn_int() %>% # Apply intercept function for covariate effects
    unlist()

  # Get population mean for Nitrogen (d15N)
  int_d15N <- select(data.frame(t(fix_coef)), starts_with("d15N")) %>%
    fn_int() %>%
    unlist()

  # Get population-level within-unit SD for Carbon
  sig_d13C <- select(data.frame(t(fix_coef)), starts_with("sigma_d13C")) %>%
    fn_int() %>%
    unlist()

  # Get population-level within-unit SD for Nitrogen
  sig_d15N <- select(data.frame(t(fix_coef)), starts_with("sigma_d15N")) %>%
    fn_int() %>%
    unlist()

  # Extract BLUPs (unit-specific deviations from population mean)
  blup <- ranef(model, summary = TRUE)[[r_e]]
  blups <- blup[, 1, 1:4] # Get estimates for 4 parameters per unit

  # Combine BLUPs with fixed effects to get unit-specific predictions
  bl_fix <- data.frame(blups) %>%
    mutate(
      # Unit-specific mean isotope values (BLUP + fixed effect)
      int_d13C = d13C_Intercept + int_d13C,
      int_d15N = d15N_Intercept + int_d15N,
      # Unit-specific within-unit variance (back-transformed from log scale)
      sig_d13C = exp(sigma_d13C_Intercept + sig_d13C),
      sig_d15N = exp(sigma_d15N_Intercept + sig_d15N)
    )

  # Create bi-plot with error bars representing within-unit variance
  ggplot(bl_fix, aes(x = int_d13C, y = int_d15N)) +
    geom_point() + # Points = unit-specific mean isotope values
    # Horizontal error bars = within-unit variance in Carbon
    geom_errorbar(
      aes(
        xmin = int_d13C - sig_d13C,
        xmax = int_d13C + sig_d13C
      ),
      color = "darkblue"
    ) +
    # Vertical error bars = within-unit variance in Nitrogen
    geom_errorbar(
      aes(
        ymin = int_d15N - sig_d15N,
        ymax = int_d15N + sig_d15N
      ),
      color = "coral"
    ) +
    xlab("d13C") +
    ylab("d15N") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 25, face = "bold"),
      axis.title.y = element_text(size = 25, face = "bold")
    ) +
    theme(axis.text = element_text(size = 18))
}


# ========== ELLIPSE PLOTTING FUNCTION ==========

# Create bivariate plot with confidence ellipses showing uncertainty in random effects
# Ellipses represent covariance structure between two parameters
#
# Arguments:
#   v: variance-covariance matrix (not directly used; function uses 'vars' from parent environment)
#   df_blups: data frame of BLUPs (random effects)
#   coefs: character vector of length 2 with parameter names to plot
#   n_ell: number of posterior draws to use for bootstrap ellipses
#   alpha_ell: transparency level for bootstrap ellipses (default 0.05)
#   x_axs, y_axs: axis labels
#   title: plot title (not currently used)
#   axs_lim: axis limits (will be -axs_lim to +axs_lim)
#
# Returns:
#   ggplot object with mean ellipse, bootstrap ellipses, and BLUP points
plot_blup_ell <- function(
  v,
  df_blups,
  coefs,
  n_ell,
  alpha_ell = 0.05,
  x_axs,
  y_axs,
  title,
  axs_lim
) {
  # Extract 2×2 variance-covariance matrix for the two parameters of interest
  # Selects var(param1), cov(param1,param2), cov(param2,param1), var(param2)
  v_est <- vars %>%
    select(
      paste(coefs[1], coefs[1], sep = "."), # Variance of parameter 1
      paste(coefs[1], coefs[2], sep = "."), # Covariance
      paste(coefs[2], coefs[1], sep = "."), # Covariance (symmetric)
      paste(coefs[2], coefs[2], sep = ".") # Variance of parameter 2
    )

  # Determine how many posterior draws to use for bootstrap ellipses
  n_boot <- ifelse(!missing(n_ell), n_ell, nrow(v_est))
  samples <- sample(seq_len(nrow(v_est)), n_boot, replace = FALSE)

  # Number of coordinate points to define each ellipse curve
  n_ellpoints <- 100

  # Create mean ellipse using average variance-covariance matrix
  # This represents the central estimate of the covariance structure
  df_ell <- as_tibble(ellipse::ellipse(
    matrix(colMeans(v_est), 2, 2), # 2×2 covariance matrix
    npoints = n_ellpoints
  )) # Using ellipse package (not car package)

  # Initialize empty data frame for bootstrap ellipses
  # Will contain ellipses from multiple posterior draws to show uncertainty
  df_bootell <- tibble(
    Replicate = numeric(0),
    x = numeric(0),
    y = numeric(0)
  )

  # Loop through selected posterior draws to create bootstrap ellipses
  # Each ellipse represents uncertainty in the covariance structure
  for (i in samples) {
    # Convert this posterior draw to 2×2 covariance matrix
    mat_tmp <- matrix(unlist(v_est[i, ]), 2, 2)
    # Generate ellipse coordinates for this draw
    df_tmp <- cbind(
      Replicate = rep(i, n_ellpoints),
      as_tibble(ellipse::ellipse(mat_tmp, npoints = n_ellpoints))
    )

    # Append to bootstrap ellipse data frame
    df_bootell <- bind_rows(df_bootell, df_tmp)
  }

  # Create ggplot with ellipses and BLUPs
  gg_blup_ellipse <- ggplot() +
    # Add reference lines at zero (centered parameters)
    geom_hline(
      yintercept = 0,
      colour = "grey30",
      alpha = 0.4,
      linetype = "dotted"
    ) +
    geom_vline(
      xintercept = 0,
      colour = "grey30",
      alpha = 0.4,
      linetype = "dotted"
    ) +
    # Add bootstrap ellipses (light grey, semi-transparent) showing uncertainty
    geom_path(
      data = df_bootell,
      aes(
        x = x,
        y = y,
        group = factor(Replicate) # Each replicate is a separate ellipse
      ),
      colour = "grey50",
      alpha = alpha_ell # Low alpha for many overlapping ellipses
    ) +
    # Add mean ellipse (thick black line) showing central estimate
    geom_path(
      data = df_ell,
      aes(
        x = x,
        y = y
      ),
      linewidth = 1.2
    ) +
    # Add BLUP points showing unit-specific estimates
    geom_point(
      data = df_blups,
      aes(x = .data[[coefs[1]]], y = .data[[coefs[2]]]),
      alpha = 0.4
    ) +
    # Set symmetric axis limits
    xlim(-axs_lim, axs_lim) +
    ylim(-axs_lim, axs_lim) +
    labs(
      x = x_axs,
      y = y_axs,
    ) +
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 25, face = "bold")
    )

  # Return the plot
  gg_blup_ellipse
}

# ========== MULTI-PANEL ELLIPSE PLOT ==========

# Create 6-panel figure showing all pairwise relationships between DHGLM parameters
# Shows covariance structure between mean and variance components for both isotopes
#
# Arguments:
#   blup: data frame of BLUPs (random effects) with columns d13C, d15N, sigma_d13C, sigma_d15N
#   vars: variance-covariance matrix for random effects
#
# Returns:
#   Combined patchwork plot with 6 panels showing all parameter pairs
plots_ell <- function(blup, vars) {
  # Panel 1: Mean d13C vs mean d15N (niche position)
  p1 <- plot_blup_ell(
    vars,
    blup,
    c("d13C", "d15N"),
    200, # Use 200 bootstrap ellipses
    x_axs = "d13C",
    y_axs = "d15N",
    axs_lim = 3.5
  )

  # Panel 2: Mean d13C vs within-unit variance in d15N
  p2 <- plot_blup_ell(
    vars,
    blup,
    c("d13C", "sigma_d15N"),
    200,
    x_axs = "d13C",
    y_axs = "sigma d15N",
    axs_lim = 3.5
  )

  # Panel 3: Mean d13C vs within-unit variance in d13C
  p3 <- plot_blup_ell(
    vars,
    blup,
    c("d13C", "sigma_d13C"),
    200,
    x_axs = "d13C",
    y_axs = "sigma d13C",
    axs_lim = 3.5
  )

  # Panel 4: Within-unit variance in d13C vs mean d15N
  p4 <- plot_blup_ell(
    vars,
    blup,
    c("sigma_d13C", "d15N"),
    200,
    x_axs = "sigma_d13C",
    y_axs = "d15N",
    axs_lim = 3.5
  )

  # Panel 5: Within-unit variance in d13C vs within-unit variance in d15N (niche width)
  p5 <- plot_blup_ell(
    vars,
    blup,
    c("sigma_d13C", "sigma_d15N"),
    200,
    x_axs = "sigma_d13C",
    y_axs = "sigma_d15N",
    axs_lim = 3.5
  )

  # Panel 6: Mean d15N vs within-unit variance in d15N
  p6 <- plot_blup_ell(
    vars,
    blup,
    c("d15N", "sigma_d15N"),
    200,
    x_axs = "d15N",
    y_axs = "sigma_d15N",
    axs_lim = 3.5
  )

  # Combine all 6 panels using patchwork
  (p1 + p2 + p3 + p4 + p5 + p6)
}
