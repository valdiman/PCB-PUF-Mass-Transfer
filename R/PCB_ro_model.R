# PCB AIR-PUF TRANSFER COEFFICIENT MODEL

# Install packages
{
  install.packages("dplyr")
  install.packages("ggplot2")
}

# Upload libraries
{
  library(dplyr)
  library(ggplot2)
}

# 1. Input data -----------------------------------------------------------
# Independently estimated uptake coefficients with 95% CI
ro_data <- data.frame(
  congener = c("PCB4", "PCB18", "PCB52"),
  MW = c(188.65, 257.54, 291.99),
  ro = c(1619.59, 1507.26, 1206.89),
  SE = c(75.63, 133.03, 97.38),
  CI_lower = c(1445.18, 1200.49, 982.32),
  CI_upper = c(1793.99, 1814.04, 1431.46)
)

# Preserve congener order
ro_data$congener <- factor(
  ro_data$congener,
  levels = c("PCB4", "PCB18", "PCB52")
)

# 2. Plot 1: ro vs PCB ----------------------------------------------------
p_ro_pcb <- ggplot(ro_data, aes(x = congener, y = ro)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = round(ro, 1)), vjust = -0.5,
            size = 4) +
  scale_y_continuous(limits = c(0, 1900),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Air-PUF mass-transfer coefficients",
       x = "PCB congener", y = expression(r[o] ~ "(cm/day)")) +
  theme_bw(base_size = 12)

# see plot
p_ro_pcb

# 3. Plot 2: ro vs MW -----------------------------------------------------
p_ro_MW <- ggplot(ro_data, aes(x = MW, y = ro)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(aes(label = round(ro, 1)), vjust = -0.5,
            size = 4) +
  scale_y_continuous(limits = c(0, 1900),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Air-PUF mass-transfer coefficients",
       x = "MW", y = expression(r[o] ~ "(cm/day)")) +
  theme_bw(base_size = 12)

p_ro_MW

# 4. MW based air diffusivity ---------------------------------------------
# Air diffusion coefficients
# Reference parameters
D_H2O_air <- 0.2743615  # cm2/s
MW_H2O <- 18.015         # g/mol

# Calculate PCB diffusion coefficients in air
ro_data <- ro_data %>%
  mutate(D_air = D_H2O_air * (MW / MW_H2O)^(-0.5))

print(ro_data)

# 5. Plot 3: ro vs MW-based air diffusivity -------------------------------
p_ro_Dair <- ggplot(ro_data, aes(x = D_air, y = ro)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(aes(label = round(ro, 1)), vjust = -0.5,
            size = 4) +
  scale_y_continuous(limits = c(0, 1900),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Air-PUF mass-transfer coefficient vs. air diffusivity",
       x = expression(D[air] ~ "(cm"^2*"/s)"), y = expression(r[o] ~ "(cm/day)")) +
  theme_bw(base_size = 12)

p_ro_Dair

# 6. Fuller diffusion volumes ---------------------------------------------
# Molecular composition
ro_data <- ro_data %>%
  mutate(
    n_C = 12,
    n_H = c(8, 7, 6),
    n_Cl = c(2, 3, 4),
    n_rings = 2)

# Fuller atomic diffusion-volume contributions
V_C <- 16.5
V_H <- 2.0
V_Cl <- 19.5
V_ring <- -20.2

# Calculate diffusion volumes
ro_data <- ro_data %>%
  mutate(
    V_D = n_C * V_C +
      n_H * V_H +
      n_Cl * V_Cl +
      n_rings * V_ring)

# Experimental conditions
T_K <- 298.15
P_atm <- 1

# Air properties
MW_air <- 28.97
V_D_air <- 20.1

# Fuller molecular diffusivity
# Equation 18-44
ro_data <- ro_data %>%
  mutate(
    D_air_Fuller =
      (1e-3 * T_K^1.75 *
         sqrt(1 / MW_air + 1 / MW)) /
      (P_atm *
         (V_D_air^(1/3) + V_D^(1/3))^2))

print(ro_data)

# 7. Plot 4: ro vs Fuller air diffusivity ---------------------------------
ro_data$congener <- factor(ro_data$congener,
                           levels = c("PCB4", "PCB18", "PCB52"))

p_ro_Fuller <- ggplot(ro_data, aes(x = congener, y = ro)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text(aes(label = paste0("D = ", round(D_air_Fuller, 4))),
            vjust = -1, size = 3.8) +
  scale_y_continuous(limits = c(1000, 1800)) +
  labs(
    title = "Air-PUF mass-transfer coefficients",
    subtitle = "Air diffusivity estimated using Fuller et al. (1966)",
    x = "PCB congener", y = expression(r[o] ~ "(cm/day)")) +
  theme_bw(base_size = 12)

p_ro_Fuller

# 8. Theoretical Air-PUF transfer Model -----------------------------------
# ro_i = ro_PCB4 * (D_i / D_PCB4)^0.67
# This is a theoretical comparison, not a new fit.
ro_ref <- 1619.59

D_ref <- ro_data$D_air_Fuller[ro_data$congener == "PCB4"]

# Theoretical diffusivity exponent
n_theoretical <- 0.67

# Calculate predicted coefficients
ro_data <- ro_data %>%
  mutate(ro_diffusion = ro_ref * (D_air_Fuller / D_ref)^n_theoretical,
         difference = ro - ro_diffusion,
         difference_percent = 100 * difference / ro_diffusion)

print(ro_data %>%
        select(congener, D_air_Fuller, ro, ro_diffusion,
               difference_percent))

# 9. Common effective air-film model --------------------------------------
# Common effective air-film thickness (cm)
# Estimated previously using all three congeners
delta_common_cm <- 3.212345

# PUF exposed area (cm2)
Apuf_cm2 <- 7.07

# Calculate common diffusion-based ro for all congeners
ro_data <- ro_data %>%
  mutate(
    # Predicted air-PUF coefficient (cm/day)
    ro_common = 86400 * D_air_Fuller / delta_common_cm,
    # Equivalent sampling rate (m3/day)
    Rs_common_m3_day = ro_common * Apuf_cm2 / 1e6,
    # Difference between individually fitted and predicted ro
    difference_common = ro - ro_common,
    # Percentage difference
    difference_common_percent =
      100 * (ro - ro_common) / ro_common)

# Display results
cat("\nCOMMON DIFFUSION MODEL\n")
cat("Common effective thickness:",
    delta_common_cm * 10, "mm\n")

print(ro_data %>%
        select(congener,
               D_air_Fuller,
               ro,
               ro_common,
               Rs_common_m3_day,
               difference_common_percent))

# 10. Plot 5: Fitted vs Common diffusion model ----------------------------
# Prepare comparison data
p_ro_comparison <- ggplot(ro_data, aes(x = congener)) +
  # Approximate 95% CI for individually fitted ro
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.12, linewidth = 0.8, color = "steelblue") +
  # Individually fitted ro
  geom_point(aes(y = ro, color = "Individually fitted"),
             size = 4, shape = 16) +
  # Common diffusion model
  geom_point(aes(y = ro_common, color = "Common diffusion model"),
             size = 4, shape = 17) +
  scale_color_manual(values = c("Individually fitted" = "steelblue",
                                "Common diffusion model" = "red")) +
  scale_y_continuous(limits = c(0, 2200), breaks = seq(0, 2200, 200),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(title = "PCB air-PUF mass-transfer coefficients",
       subtitle = "Individual estimates (95% CI) vs. common diffusion model",
       caption = "Common effective diffusion length = 32.12 mm",
       x = "PCB congener", y = expression(r[o] ~ "(cm/day)"), color = "Model") +
  theme_bw(base_size = 12)

print(p_ro_comparison)

# 12. Weighted common diffusion model -------------------------------------
# Weighted objective function
objective_delta_weighted <- function(delta_cm) {
  ro_predicted <- 86400 * ro_data$D_air_Fuller / delta_cm
  sum(((ro_data$ro - ro_predicted) / ro_data$SE)^2)
}

# Estimate weighted common thickness
fit_delta_weighted <- optimize(
  objective_delta_weighted,
  interval = c(0.1, 10))

delta_weighted_cm <- fit_delta_weighted$minimum
delta_weighted_mm <- delta_weighted_cm * 10

cat("\nWEIGHTED COMMON DIFFUSION MODEL\n")
cat("Unweighted thickness:", delta_common_cm * 10, "mm\n")
cat("Weighted thickness:", delta_weighted_mm, "mm\n")
cat("Weighted objective:", fit_delta_weighted$objective, "\n")

# 13. Uncertainty in weighted thickness -----------------------------------
# Predictor and inverse-variance weights
x <- 86400 * ro_data$D_air_Fuller
w <- 1 / ro_data$SE^2

# Weighted estimate of inverse thickness
beta_hat <- sum(w * x * ro_data$ro) / sum(w * x^2)

# Standard error
SE_beta <- sqrt(1 / sum(w * x^2))

# Approximate 95% confidence interval
beta_lower <- beta_hat - qnorm(0.975) * SE_beta
beta_upper <- beta_hat + qnorm(0.975) * SE_beta

# Convert to effective thickness
delta_weighted_cm <- 1 / beta_hat

if (beta_lower <= 0) {
  stop("Confidence interval includes nonpositive inverse thickness.")
}

delta_lower_mm <- 10 / beta_upper
delta_upper_mm <- 10 / beta_lower

cat("\nWEIGHTED EFFECTIVE THICKNESS\n")
cat("Estimate:", round(delta_weighted_cm * 10, 2), "mm\n")
cat("Approximate 95% CI:", round(delta_lower_mm, 2), "to",
    round(delta_upper_mm, 2), "mm\n")

# 14. Weighted model predictions ------------------------------------------
ro_data <- ro_data %>%
  mutate(ro_weighted = 86400 * D_air_Fuller / delta_weighted_cm,
         difference_weighted_percent = 100 * (ro - ro_weighted) / ro_weighted)

print(ro_data %>%
        select(congener,
               ro,
               CI_lower,
               CI_upper,
               ro_common,
               ro_weighted,
               difference_weighted_percent))

# 15. Plot 6: Weighted vs unweighted model --------------------------------
p_ro_weighted <- ggplot(ro_data, aes(x = congener)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.12,
                linewidth = 0.8, color = "steelblue") +
  geom_point(aes(y = ro, color = "Individually fitted"), size = 4,
             shape = 16) +
  geom_point(aes(y = ro_common, color = "Unweighted diffusion"),
             size = 4, shape = 17) +
  geom_point(aes(y = ro_weighted, color = "Weighted diffusion"),
             size = 4, shape = 15) +
  scale_color_manual(values = c("Individually fitted" = "steelblue",
                                "Unweighted diffusion" = "red",
                                "Weighted diffusion" = "darkgreen")) +
  scale_y_continuous(limits = c(0, 2200), breaks = seq(0, 2200, 200),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(title = "PCB air-PUF mass-transfer coefficients",
       subtitle = "Individual estimates (95% CI) vs. common diffusion models",
       caption = paste0("Unweighted thickness = ", round(delta_common_cm * 10, 2),
                        " mm | Weighted thickness = ",
                        round(delta_weighted_cm * 10, 2), " mm"),
       x = "PCB congener", y = expression(r[o] ~ "(cm/day)"), color = "Model") +
  theme_bw(base_size = 12)

print(p_ro_weighted)


