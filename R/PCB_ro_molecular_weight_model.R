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

# Independently estimated uptake coefficients
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

# Plot 1: ro vs PCB
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

# plot 2: ro vs MW
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

# Air diffusion coefficients
# Reference parameters
D_H2O_air <- 0.2743615  # cm2/s
MW_H2O <- 18.015         # g/mol

# Calculate PCB diffusion coefficients in air
ro_data <- ro_data %>%
  mutate(D_air = D_H2O_air * (MW / MW_H2O)^(-0.5))

print(ro_data)

# Plot 3: ro vs air diffusivity
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
ro_data <- ro_data %>%
  mutate(
    D_air_Fuller =
      (1e-3 * T_K^1.75 *
         sqrt(1 / MW_air + 1 / MW)) /
      (P_atm *
         (V_D_air^(1/3) + V_D^(1/3))^2))

print(ro_data)

# Plot 4: ro vs Fuller molecular diffusivity
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


# Theoretical Air-PUF transfer model
# Reference PCB4 values
ro_ref <- 1619.601

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
      100 * (ro - ro_common) / ro_common
  )

# Display results
cat("\nCOMMON DIFFUSION MODEL\n")

cat("Common effective thickness:",
    delta_common_cm * 10, "mm\n")

print(
  ro_data %>%
    select(
      congener,
      D_air_Fuller,
      ro,
      ro_common,
      Rs_common_m3_day,
      difference_common_percent
    )
)

# Prepare comparison data
ro_comparison <- ro_data %>%
  select(congener, ro, ro_common) %>%
  tidyr::pivot_longer(
    cols = c(ro, ro_common),
    names_to = "model",
    values_to = "ro_cm_day"
  ) %>%
  mutate(
    model = recode(
      model,
      ro = "Individually fitted",
      ro_common = "Common diffusion model"
    ),
    congener = factor(
      congener,
      levels = c("PCB4", "PCB18", "PCB52")
    )
  )

# Plot
p_ro_comparison <- ggplot(
  ro_comparison,
  aes(
    x = congener,
    y = ro_cm_day,
    color = model
  )
) +
  
  geom_point(
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 1800),
    breaks = seq(0, 1800, 200),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Air-PUF mass-transfer coefficients",
    subtitle = "Individual kinetic estimates vs. common diffusion model",
    x = "PCB congener",
    y = expression(r[o] ~ "(cm/day)"),
    color = "Model"
  ) +
  
  theme_bw(base_size = 12)

print(p_ro_comparison)
