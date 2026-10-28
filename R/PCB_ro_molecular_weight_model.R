# ============================================================
# PCB AIR-PUF TRANSFER COEFFICIENT MODEL
# Molecular-weight dependence
# ============================================================

library(dplyr)
library(ggplot2)

# Independently estimated uptake coefficients
ro_data <- data.frame(
  congener = c("PCB4", "PCB18+30", "PCB52"),
  MW = c(188.65, 257.54, 291.99),
  ro = c(1619.601, 1507.263, 1206.890)
)

# Reference molecular weight
MW_ref <- 188.65

# ------------------------------------------------------------
# 1. Empirical power-law regression
# ------------------------------------------------------------

fit_power <- lm(
  log(ro) ~ log(MW / MW_ref),
  data = ro_data
)

print(summary(fit_power))

# Extract fitted parameters
ro_ref_fitted <- exp(coef(fit_power)[1])
b_fitted <- -coef(fit_power)[2]

cat("\nFitted reference ro:", ro_ref_fitted, "cm/day\n")
cat("Fitted exponent b:", b_fitted, "\n")

# ------------------------------------------------------------
# 2. Compare with theoretical relationship
# ------------------------------------------------------------

ro_data <- ro_data %>%
  mutate(
    ro_empirical =
      ro_ref_fitted * (MW / MW_ref)^(-b_fitted),
    
    ro_theoretical =
      1619.601 * (MW / MW_ref)^(-0.335),
    
    residual_empirical = ro - ro_empirical,
    
    residual_theoretical = ro - ro_theoretical
  )

print(ro_data)

# ------------------------------------------------------------
# 3. Generate smooth prediction curves
# ------------------------------------------------------------

MW_sequence <- seq(
  min(ro_data$MW),
  max(ro_data$MW),
  length.out = 200
)

prediction_curves <- data.frame(
  MW = MW_sequence
) %>%
  mutate(
    empirical =
      ro_ref_fitted * (MW / MW_ref)^(-b_fitted),
    
    theoretical =
      1619.601 * (MW / MW_ref)^(-0.335)
  )

# ------------------------------------------------------------
# 4. Plot results
# ------------------------------------------------------------

p_ro <- ggplot() +
  
  geom_line(
    data = prediction_curves,
    aes(MW, empirical),
    color = "black",
    linewidth = 1
  ) +
  
  geom_line(
    data = prediction_curves,
    aes(MW, theoretical),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_point(
    data = ro_data,
    aes(MW, ro),
    color = "blue",
    size = 3
  ) +
  
  geom_text(
    data = ro_data,
    aes(MW, ro, label = congener),
    vjust = -1,
    size = 3.5
  ) +
  
  labs(
    x = "Molecular weight (g/mol)",
    y = expression(r[o] ~ "(cm/day)"),
    title = "PCB air-PUF mass-transfer coefficients",
    subtitle = paste0(
      "Empirical exponent b = ",
      round(b_fitted, 3),
      " | Theoretical exponent = 0.335"
    )
  ) +
  
  theme_bw(base_size = 12)

print(p_ro)

# ------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------

dir.create(
  "Results/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Tables",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  ro_data,
  "Results/Tables/PCB_ro_molecular_weight.csv",
  row.names = FALSE
)

ggsave(
  "Results/Figures/PCB_ro_molecular_weight.png",
  p_ro,
  width = 8,
  height = 5,
  dpi = 300
)