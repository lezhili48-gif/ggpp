# ============================================================================
# Three-Level Meta-Analysis: Biochar Effects on DNRA and Denitrification
# ============================================================================
# Author: [Your Name]
# Date: 2025-02-25
# Description: Complete R code for three-level meta-analysis
# ============================================================================

# ============================================================================
# 1. PACKAGE INSTALLATION AND LOADING
# ============================================================================

# Install required packages (run once)
# install.packages(c("metafor", "ggplot2", "orchaRd", "clubSandwich", 
#                    "dplyr", "tidyr", "readxl", "writexl", "gridExtra"))

# Load packages
library(metafor)      # Main meta-analysis package
library(ggplot2)      # Visualization
library(orchaRd)      # Orchard plots
library(clubSandwich) # Robust variance estimation
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(readxl)       # Read Excel files
library(writexl)      # Write Excel files
library(gridExtra)    # Combine plots

# Set working directory (adjust as needed)
# setwd("/path/to/your/data")

# ============================================================================
# 2. DATA IMPORT AND PREPARATION
# ============================================================================

# Read data extraction sheet
data <- read_excel("Data_Extraction_Completed.xlsx", sheet = "Data Extraction", skip = 1)

# View data structure
str(data)
head(data)

# ============================================================================
# 3. EFFECT SIZE CALCULATION
# ============================================================================

# Function to calculate lnRR and variance
calculate_lnRR <- function(mean_t, mean_c, sd_t, sd_c, n_t, n_c) {
  # Natural log response ratio
  lnRR <- log(mean_t / mean_c)
  
  # Variance of lnRR
  var_lnRR <- (sd_t^2 / (n_t * mean_t^2)) + (sd_c^2 / (n_c * mean_c^2))
  
  return(data.frame(lnRR = lnRR, var_lnRR = var_lnRR))
}

# Calculate effect sizes for DNRA
data_dnra <- data %>%
  filter(!is.na(DNRA_Ctrl_Mean) & !is.na(DNRA_Trt_Mean)) %>%
  mutate(
    lnRR_DNRA = log(DNRA_Trt_Mean / DNRA_Ctrl_Mean),
    var_lnRR_DNRA = (DNRA_Trt_SD^2 / (DNRA_Trt_n * DNRA_Trt_Mean^2)) + 
                    (DNRA_Ctrl_SD^2 / (DNRA_Ctrl_n * DNRA_Ctrl_Mean^2)),
    outcome = "DNRA"
  )

# Calculate effect sizes for Denitrification
data_denit <- data %>%
  filter(!is.na(Denit_Ctrl_Mean) & !is.na(Denit_Trt_Mean)) %>%
  mutate(
    lnRR_Denit = log(Denit_Trt_Mean / Denit_Ctrl_Mean),
    var_lnRR_Denit = (Denit_Trt_SD^2 / (Denit_Trt_n * Denit_Trt_Mean^2)) + 
                     (Denit_Ctrl_SD^2 / (Denit_Ctrl_n * Denit_Ctrl_Mean^2)),
    outcome = "Denitrification"
  )

# Calculate effect sizes for N2O emissions
data_n2o <- data %>%
  filter(!is.na(N2O_Ctrl_Mean) & !is.na(N2O_Trt_Mean)) %>%
  mutate(
    lnRR_N2O = log(N2O_Trt_Mean / N2O_Ctrl_Mean),
    var_lnRR_N2O = (N2O_Trt_SD^2 / (N2O_Trt_n * N2O_Trt_Mean^2)) + 
                   (N2O_Ctrl_SD^2 / (N2O_Ctrl_n * N2O_Ctrl_Mean^2)),
    outcome = "N2O"
  )

# Calculate effect sizes for functional genes (using log-transformed abundances)
data_nrfA <- data %>%
  filter(!is.na(nrfA_Ctrl_Mean) & !is.na(nrfA_Trt_Mean)) %>%
  mutate(
    lnRR_nrfA = log(nrfA_Trt_Mean / nrfA_Ctrl_Mean),
    var_lnRR_nrfA = (nrfA_Trt_SD^2 / (nrfA_Trt_n * nrfA_Trt_Mean^2)) + 
                    (nrfA_Ctrl_SD^2 / (nrfA_Ctrl_n * nrfA_Ctrl_Mean^2)),
    outcome = "nrfA"
  )

data_nirS <- data %>%
  filter(!is.na(nirS_Ctrl_Mean) & !is.na(nirS_Trt_Mean)) %>%
  mutate(
    lnRR_nirS = log(nirS_Trt_Mean / nirS_Ctrl_Mean),
    var_lnRR_nirS = (nirS_Trt_SD^2 / (nirS_Trt_n * nirS_Trt_Mean^2)) + 
                    (nirS_Ctrl_SD^2 / (nirS_Ctrl_n * nirS_Ctrl_Mean^2)),
    outcome = "nirS"
  )

data_nirK <- data %>%
  filter(!is.na(nirK_Ctrl_Mean) & !is.na(nirK_Trt_Mean)) %>%
  mutate(
    lnRR_nirK = log(nirK_Trt_Mean / nirK_Ctrl_Mean),
    var_lnRR_nirK = (nirK_Trt_SD^2 / (nirK_Trt_n * nirK_Trt_Mean^2)) + 
                    (nirK_Ctrl_SD^2 / (nirK_Ctrl_n * nirK_Ctrl_Mean^2)),
    outcome = "nirK"
  )

data_nosZ <- data %>%
  filter(!is.na(nosZ_Ctrl_Mean) & !is.na(nosZ_Trt_Mean)) %>%
  mutate(
    lnRR_nosZ = log(nosZ_Trt_Mean / nosZ_Ctrl_Mean),
    var_lnRR_nosZ = (nosZ_Trt_SD^2 / (nosZ_Trt_n * nosZ_Trt_Mean^2)) + 
                    (nosZ_Ctrl_SD^2 / (nosZ_Ctrl_n * nosZ_Ctrl_Mean^2)),
    outcome = "nosZ"
  )

# ============================================================================
# 4. THREE-LEVEL META-ANALYSIS MODEL
# ============================================================================

# Function to fit three-level meta-analysis model
fit_three_level_model <- function(data, yi_col, vi_col, study_col = "Study_ID") {
  
  # Fit three-level random-effects model
  model <- rma.mv(
    yi = data[[yi_col]],
    V = data[[vi_col]],
    random = ~ 1 | data[[study_col]] / rownames(data),  # Level 3: Study, Level 2: Observation
    method = "REML",
    data = data,
    control = list(optimizer = "optim", optmethod = "Nelder-Mead")
  )
  
  return(model)
}

# Fit models for each outcome
model_dnra <- fit_three_level_model(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA")
model_denit <- fit_three_level_model(data_denit, "lnRR_Denit", "var_lnRR_Denit")
model_n2o <- fit_three_level_model(data_n2o, "lnRR_N2O", "var_lnRR_N2O")
model_nrfA <- fit_three_level_model(data_nrfA, "lnRR_nrfA", "var_lnRR_nrfA")
model_nirS <- fit_three_level_model(data_nirS, "lnRR_nirS", "var_lnRR_nirS")
model_nirK <- fit_three_level_model(data_nirK, "lnRR_nirK", "var_lnRR_nirK")
model_nosZ <- fit_three_level_model(data_nosZ, "lnRR_nosZ", "var_lnRR_nosZ")

# Print model summaries
cat("\n=== DNRA Model Summary ===\n")
summary(model_dnra)

cat("\n=== Denitrification Model Summary ===\n")
summary(model_denit)

cat("\n=== N2O Model Summary ===\n")
summary(model_n2o)

# ============================================================================
# 5. HETEROGENEITY ANALYSIS
# ============================================================================

# Function to calculate I-squared for three-level model
calculate_I2 <- function(model) {
  # Total variance
  W <- diag(1 / model$vi)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  
  # Variance components
  sigma2_study <- model$sigma2[1]  # Between-study variance (Level 3)
  sigma2_obs <- model$sigma2[2]     # Within-study variance (Level 2)
  
  # Total sampling variance
  total_sampling_var <- (model$k - model$p) / sum(diag(P))
  
  # I-squared calculations
  I2_total <- (sigma2_study + sigma2_obs) / (sigma2_study + sigma2_obs + total_sampling_var) * 100
  I2_study <- sigma2_study / (sigma2_study + sigma2_obs + total_sampling_var) * 100
  I2_obs <- sigma2_obs / (sigma2_study + sigma2_obs + total_sampling_var) * 100
  
  return(data.frame(
    I2_total = I2_total,
    I2_study = I2_study,
    I2_obs = I2_obs
  ))
}

# Calculate I-squared for each model
I2_dnra <- calculate_I2(model_dnra)
I2_denit <- calculate_I2(model_denit)
I2_n2o <- calculate_I2(model_n2o)

cat("\n=== Heterogeneity (I²) ===\n")
cat("DNRA - Total:", round(I2_dnra$I2_total, 1), "%, Study:", round(I2_dnra$I2_study, 1), "%, Obs:", round(I2_dnra$I2_obs, 1), "%\n")
cat("Denitrification - Total:", round(I2_denit$I2_total, 1), "%, Study:", round(I2_denit$I2_study, 1), "%, Obs:", round(I2_denit$I2_obs, 1), "%\n")
cat("N2O - Total:", round(I2_n2o$I2_total, 1), "%, Study:", round(I2_n2o$I2_study, 1), "%, Obs:", round(I2_n2o$I2_obs, 1), "%\n")

# ============================================================================
# 6. PUBLICATION BIAS ASSESSMENT
# ============================================================================

# Function to assess publication bias
assess_pub_bias <- function(data, yi_col, vi_col) {
  yi <- data[[yi_col]]
  vi <- data[[vi_col]]
  
  # Funnel plot
  funnel(yi, vi, main = "Funnel Plot")
  
  # Egger's regression test
  egger <- regtest(yi, vi, model = "lm")
  
  # Rosenberg's fail-safe number
  fsnum <- fsn(yi, vi)
  
  return(list(
    egger_test = egger,
    fail_safe_N = fsnum
  ))
}

# Assess publication bias for each outcome
cat("\n=== Publication Bias Assessment ===\n")

pub_bias_dnra <- assess_pub_bias(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA")
cat("\nDNRA Egger's test p-value:", pub_bias_dnra$egger_test$pval, "\n")
cat("DNRA Fail-safe N:", pub_bias_dnra$fail_safe_N$fsnum, "\n")

pub_bias_denit <- assess_pub_bias(data_denit, "lnRR_Denit", "var_lnRR_Denit")
cat("\nDenitrification Egger's test p-value:", pub_bias_denit$egger_test$pval, "\n")
cat("Denitrification Fail-safe N:", pub_bias_denit$fail_safe_N$fsnum, "\n")

# ============================================================================
# 7. MODERATOR ANALYSIS (META-REGRESSION)
# ============================================================================

# Function for moderator analysis
moderator_analysis <- function(data, yi_col, vi_col, moderator, study_col = "Study_ID") {
  
  formula <- as.formula(paste("~", moderator))
  
  model <- rma.mv(
    yi = data[[yi_col]],
    V = data[[vi_col]],
    mods = formula,
    random = ~ 1 | data[[study_col]] / rownames(data),
    method = "REML",
    data = data
  )
  
  return(model)
}

# Categorical moderator: Soil pH class
data_dnra$soil_pH_class <- cut(data_dnra$Soil_pH, 
                                breaks = c(-Inf, 6.0, 7.5, Inf),
                                labels = c("Acidic", "Neutral", "Alkaline"))

mod_ph_dnra <- moderator_analysis(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "soil_pH_class")
cat("\n=== Moderator: Soil pH (DNRA) ===\n")
summary(mod_ph_dnra)

# Continuous moderator: Soil pH
mod_ph_cont <- moderator_analysis(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "Soil_pH")
cat("\n=== Moderator: Soil pH (continuous) ===\n")
summary(mod_ph_cont)

# Continuous moderator: Biochar application rate
mod_rate <- moderator_analysis(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "App_Rate_t_ha")
cat("\n=== Moderator: Application Rate ===\n")
summary(mod_rate)

# Continuous moderator: MAT
mod_mat <- moderator_analysis(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "MAT_C")
cat("\n=== Moderator: Mean Annual Temperature ===\n")
summary(mod_mat)

# ============================================================================
# 8. SUBGROUP ANALYSIS
# ============================================================================

# Function for subgroup analysis
subgroup_analysis <- function(data, yi_col, vi_col, subgroup_col, study_col = "Study_ID") {
  
  subgroups <- unique(data[[subgroup_col]])
  results <- data.frame()
  
  for (sg in subgroups) {
    if (is.na(sg)) next
    
    sg_data <- data[data[[subgroup_col]] == sg, ]
    if (nrow(sg_data) < 3) next  # Skip if too few studies
    
    model <- rma.mv(
      yi = sg_data[[yi_col]],
      V = sg_data[[vi_col]],
      random = ~ 1 | sg_data[[study_col]] / rownames(sg_data),
      method = "REML",
      data = sg_data
    )
    
    results <- rbind(results, data.frame(
      Subgroup = sg,
      k = model$k,
      lnRR = model$b[1],
      ci_lb = model$ci.lb,
      ci_ub = model$ci.ub,
      p_value = model$pval[1],
      I2 = calculate_I2(model)$I2_total
    ))
  }
  
  return(results)
}

# Subgroup analysis by biochar feedstock
subgroup_feedstock <- subgroup_analysis(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "Feedstock")
cat("\n=== Subgroup: Biochar Feedstock ===\n")
print(subgroup_feedstock)

# Subgroup analysis by pyrolysis temperature
data_dnra$pyro_class <- cut(data_dnra$Pyro_T_C, 
                             breaks = c(-Inf, 400, 600, Inf),
                             labels = c("Low (<400)", "Medium (400-600)", "High (>600)"))

subgroup_pyro <- subgroup_analysis(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "pyro_class")
cat("\n=== Subgroup: Pyrolysis Temperature ===\n")
print(subgroup_pyro)

# ============================================================================
# 9. VISUALIZATION
# ============================================================================

# 9.1 Forest Plot for Overall Effects
forest_overall <- function(models, labels) {
  
  results <- data.frame(
    Outcome = labels,
    lnRR = sapply(models, function(m) m$b[1]),
    ci_lb = sapply(models, function(m) m$ci.lb),
    ci_ub = sapply(models, function(m) m$ci.ub),
    k = sapply(models, function(m) m$k)
  )
  
  results$percent_change <- (exp(results$lnRR) - 1) * 100
  
  p <- ggplot(results, aes(x = lnRR, y = reorder(Outcome, lnRR))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.2, linewidth = 0.8) +
    geom_point(size = 4, color = "#2E7D32") +
    geom_text(aes(label = paste0("+", round(percent_change, 1), "%")), 
              vjust = -0.8, size = 3.5, fontface = "bold", color = "#2E7D32") +
    geom_text(aes(label = paste0("n=", k)), 
              vjust = 1.5, size = 3, color = "gray40") +
    labs(x = "Effect size (lnRR)", y = "") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12)
    )
  
  return(p)
}

# Create forest plot
models_list <- list(model_dnra, model_denit, model_n2o, model_nrfA, model_nosZ)
labels_list <- c("DNRA rate", "Denitrification rate", "N2O emission", "nrfA abundance", "nosZ abundance")

p_forest <- forest_overall(models_list, labels_list)
ggsave("Forest_Plot_Overall.png", p_forest, width = 8, height = 5, dpi = 300)

# 9.2 Orchard Plot (using orchaRd package)
orchard_plot <- orchaRd::orchard_plot(
  model_dnra,
  mod = "1",
  group = "Study_ID",
  xlab = "Effect size (lnRR)",
  transfm = "none"
)
ggsave("Orchard_Plot_DNRA.png", orchard_plot, width = 8, height = 6, dpi = 300)

# 9.3 Bubble Plot for Continuous Moderators
bubble_plot <- function(data, yi_col, vi_col, moderator_col, xlabel) {
  
  p <- ggplot(data, aes_string(x = moderator_col, y = yi_col)) +
    geom_point(aes(size = 1/sqrt(data[[vi_col]])), alpha = 0.5, color = "#1565C0") +
    geom_smooth(method = "lm", color = "#C62828", se = TRUE, fill = "#FFCDD2") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = xlabel, y = "Effect size (lnRR)", size = "Precision (1/SE)") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

p_bubble_ph <- bubble_plot(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "Soil_pH", "Soil pH")
ggsave("Bubble_Plot_pH.png", p_bubble_ph, width = 7, height = 5, dpi = 300)

p_bubble_rate <- bubble_plot(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "App_Rate_t_ha", "Biochar application rate (t/ha)")
ggsave("Bubble_Plot_Rate.png", p_bubble_rate, width = 7, height = 5, dpi = 300)

# 9.4 Funnel Plot
funnel_plot <- function(data, yi_col, vi_col, title) {
  yi <- data[[yi_col]]
  vi <- data[[vi_col]]
  sei <- sqrt(vi)
  
  p <- ggplot(data, aes_string(x = yi, y = "1/sqrt(data[[vi_col]])")) +
    geom_point(alpha = 0.5, color = "#1565C0") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "Effect size (lnRR)", y = "Precision (1/SE)", title = title) +
    theme_minimal()
  
  return(p)
}

p_funnel_dnra <- funnel_plot(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA", "Funnel Plot: DNRA")
ggsave("Funnel_Plot_DNRA.png", p_funnel_dnra, width = 6, height = 6, dpi = 300)

# ============================================================================
# 10. SENSITIVITY ANALYSIS
# ============================================================================

# 10.1 Leave-one-out analysis
leave_one_out <- function(data, yi_col, vi_col, study_col = "Study_ID") {
  
  studies <- unique(data[[study_col]])
  results <- data.frame()
  
  for (study in studies) {
    subset_data <- data[data[[study_col]] != study, ]
    
    model <- rma.mv(
      yi = subset_data[[yi_col]],
      V = subset_data[[vi_col]],
      random = ~ 1 | subset_data[[study_col]] / rownames(subset_data),
      method = "REML",
      data = subset_data
    )
    
    results <- rbind(results, data.frame(
      Excluded_Study = study,
      k = model$k,
      lnRR = model$b[1],
      ci_lb = model$ci.lb,
      ci_ub = model$ci.ub
    ))
  }
  
  return(results)
}

# Run leave-one-out analysis
loo_dnra <- leave_one_out(data_dnra, "lnRR_DNRA", "var_lnRR_DNRA")
cat("\n=== Leave-One-Out Analysis (DNRA) ===\n")
print(loo_dnra)

# 10.2 Exclude outliers (|lnRR| > 3)
data_dnra_no_outliers <- data_dnra %>% filter(abs(lnRR_DNRA) <= 3)
model_dnra_no_outliers <- fit_three_level_model(data_dnra_no_outliers, "lnRR_DNRA", "var_lnRR_DNRA")

cat("\n=== Sensitivity: Excluding Outliers ===\n")
cat("Original model - lnRR:", round(model_dnra$b[1], 3), "CI:", round(model_dnra$ci.lb, 3), "-", round(model_dnra$ci.ub, 3), "\n")
cat("No outliers - lnRR:", round(model_dnra_no_outliers$b[1], 3), "CI:", round(model_dnra_no_outliers$ci.lb, 3), "-", round(model_dnra_no_outliers$ci.ub, 3), "\n")

# ============================================================================
# 11. EXPORT RESULTS
# ============================================================================

# Create summary table
summary_table <- data.frame(
  Outcome = c("DNRA rate", "Denitrification rate", "N2O emission", 
              "nrfA abundance", "nirS abundance", "nirK abundance", "nosZ abundance"),
  k = c(model_dnra$k, model_denit$k, model_n2o$k,
        model_nrfA$k, model_nirS$k, model_nirK$k, model_nosZ$k),
  lnRR = c(model_dnra$b[1], model_denit$b[1], model_n2o$b[1],
           model_nrfA$b[1], model_nirS$b[1], model_nirK$b[1], model_nosZ$b[1]),
  ci_lb = c(model_dnra$ci.lb, model_denit$ci.lb, model_n2o$ci.lb,
            model_nrfA$ci.lb, model_nirS$ci.lb, model_nirK$ci.lb, model_nosZ$ci.lb),
  ci_ub = c(model_dnra$ci.ub, model_denit$ci.ub, model_n2o$ci.ub,
            model_nrfA$ci.ub, model_nirS$ci.ub, model_nirK$ci.ub, model_nosZ$ci.ub),
  p_value = c(model_dnra$pval[1], model_denit$pval[1], model_n2o$pval[1],
              model_nrfA$pval[1], model_nirS$pval[1], model_nirK$pval[1], model_nosZ$pval[1]),
  I2_total = c(I2_dnra$I2_total, I2_denit$I2_total, I2_n2o$I2_total, NA, NA, NA, NA)
)

summary_table$percent_change <- (exp(summary_table$lnRR) - 1) * 100

# Export to Excel
write_xlsx(
  list(
    Summary_Results = summary_table,
    Subgroup_Feedstock = subgroup_feedstock,
    Subgroup_Pyrolysis = subgroup_pyro,
    Leave_One_Out = loo_dnra
  ),
  "Meta_Analysis_Results.xlsx"
)

cat("\n=== Results exported to Meta_Analysis_Results.xlsx ===\n")
print(summary_table)

# ============================================================================
# 12. SESSION INFO
# ============================================================================

sessionInfo()

# ============================================================================
# END OF SCRIPT
# ============================================================================
