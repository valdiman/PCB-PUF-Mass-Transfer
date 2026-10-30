# PCB52 PUF uptake
# Run from project root containing Data/Processed/*.csv
# Model B is primary; Model A is a reversible diagnostic.

# Install packages
{
  install.packages("dplyr")
  install.packages("ggplot2")
  install.packages("deSolve")
}

# Upload libraries
{
  library(dplyr)
  library(ggplot2)
  library(deSolve)
}

# 1. Inputs ---------------------------------------------------------------
water_file <- 'Data/Processed/water_data.csv'
puf_file   <- 'Data/Processed/puf_data.csv'
congener_id <- 'PCB52'
Cw0_ng_L <- 4538.16                 # measured initial concentration
Va_cm3 <- 135.98                     # modeled lower air compartment
Vpuf_cm3 <- 29
Aaw_cm2 <- 56.74502
Apuf_cm2 <- 7.07
Kaw <- 10^(-1.88455)              # air/water concentration ratio
Kpuf <- 10^(0.6366 * 7.724555 - 3.1774) * 21300
# Weak-stirring water-air transfer coefficient, calculated per congener.
# Same reference assumptions as the PCB4 analysis (Schmidt scaling).
MW_pcb <- 291.99
MW_H2O <- 18.015
MW_CO2 <- 44.01
D_H2O_air_cm2_s <- 0.2743615
D_CO2_water_cm2_s <- 1.67606e-5
nu_water_cm2_s <- 0.010072884
KG_H2O_m_s <- 0.003
film_CO2_um <- 127
D_pcb_air <- D_H2O_air_cm2_s * (MW_pcb / MW_H2O)^(-0.5)
D_pcb_water <- D_CO2_water_cm2_s * (MW_pcb / MW_CO2)^(-0.5)
Sc_pcb <- nu_water_cm2_s / D_pcb_water
KG_pcb_m_s <- KG_H2O_m_s * (D_pcb_air / D_H2O_air_cm2_s)^0.67
KL_CO2_m_s <- (D_CO2_water_cm2_s / (film_CO2_um * 1e-4)) / 100
KL_pcb_m_s <- KL_CO2_m_s * (Sc_pcb / 600)^(-0.5)
kaw_o_cm_day <- (1 / (1 / (KG_pcb_m_s * Kaw) + 1 / KL_pcb_m_s)) * 100 * 86400
cat('Calculated weak-stirring kaw_o:', kaw_o_cm_day, 'cm/day\n')
fit_end_hr <- 74                  # includes observations around 73.1 h
ro_bounds <- c(100, 3000)
# Compare the independently fitted ro with the PCB4 diffusivity-scaled prediction.
# PCB4-based diffusivity prediction is a comparison only, NOT the fitted value.
ro_reference <- 1619.601 * (MW_pcb / 188.65)^(-0.335)

water_data <- read.csv(water_file)
puf_data <- read.csv(puf_file)
stopifnot(all(c('congener','time_hr','Cw_ng_L','sample_id') %in% names(water_data)))
stopifnot(all(c('congener','time_hr','PCB_value') %in% names(puf_data)))
# Verify the exact congener label in both input files.
if (!(congener_id %in% water_data$congener) ||
    !(congener_id %in% puf_data$congener)) {
  stop(paste0('Congener label not found: ', congener_id,
              '\nWater labels: ', paste(unique(water_data$congener), collapse=', '),
              '\nPUF labels: ', paste(unique(puf_data$congener), collapse=', ')))
}

# 2. Water data: full-period, grouped by sampling period ------------------
water_obs <- water_data %>%
  filter(congener == congener_id, is.finite(time_hr), is.finite(Cw_ng_L)) %>%
  mutate(time_group = sub('^(T[^_]+)_.*$', '\\1', sample_id))
water_mean <- water_obs %>% group_by(time_group) %>%
  summarise(time_hr = mean(time_hr), Cw_ng_L = mean(Cw_ng_L),
            .groups = 'drop') %>%
  bind_rows(data.frame(time_group = 'T0', time_hr = 0,
                       Cw_ng_L = Cw0_ng_L)) %>% arrange(time_hr)
stopifnot(all(water_mean$Cw_ng_L > 0), all(diff(water_mean$time_hr) > 0))
# Piecewise exponential interpolation: linear interpolation in log(Cw).
log_Cw <- approxfun(water_mean$time_hr, log(water_mean$Cw_ng_L), rule = 2)
Cw_exp <- function(time_hr) exp(log_Cw(time_hr))

puf_obs <- puf_data %>% filter(congener == congener_id,
  is.finite(time_hr), is.finite(PCB_value)) %>% arrange(time_hr)
puf_fit <- puf_obs %>% filter(time_hr <= fit_end_hr)
stopifnot(nrow(puf_fit) > 0)

# 3. Full transient air + PUF model ---------------------------------------
# t in days; concentration in ng/L; transfer coefficients in cm/day.
# Water is a measured forcing function, not a predicted compartment.
kinetic_rhs <- function(t, state, parms) {
  Ca <- unname(state['Ca']); Cpuf <- unname(state['Cpuf'])
  ro <- parms$ro
  Cw <- parms$Cw_fun(t * 24)
  driving <- if (parms$reversible) Ca - Cpuf / Kpuf else Ca
  flux_wa <- kaw_o_cm_day * Aaw_cm2 * (Cw - Ca / Kaw) / 1000 # ng/day
  flux_ap <- ro * Apuf_cm2 * driving / 1000                  # ng/day
  dCa <- 1000 * (flux_wa - flux_ap) / Va_cm3
  dCpuf <- 1000 * flux_ap / Vpuf_cm3
  list(c(dCa, dCpuf))
}

simulate_puf <- function(times_hr, ro, model = 'B', Cw_fun = Cw_exp) {
  stopifnot(model %in% c('A','B'), is.finite(ro), ro > 0)
  # Include all water interpolation knots to handle slope changes.
  target <- sort(unique(c(0, times_hr)))
  knots <- water_mean$time_hr
  times <- sort(unique(c(target, knots[knots > 0 & knots < max(target)]))) / 24
  raw <- ode(y = c(Ca = 0, Cpuf = 0), times = times,
             func = kinetic_rhs,
             parms = list(ro = ro, reversible = model == 'A', Cw_fun = Cw_fun),
             method = 'lsoda', rtol = 1e-7, atol = 1e-9)
  out <- as.data.frame(raw)
  out$time_hr <- out$time * 24
  out$Mpuf_ng <- out$Cpuf * Vpuf_cm3 / 1000
  out$Ma_ng <- out$Ca * Va_cm3 / 1000
  if (any(!is.finite(out$Mpuf_ng))) stop('Non-finite simulation')
  out
}

predict_at <- function(sim, times_hr) {
  approx(sim$time_hr, sim$Mpuf_ng, xout = times_hr, rule = 2)$y
}
objective <- function(ro, model = 'B') {
  if (!is.finite(ro) || ro <= 0) return(1e20)
  sim <- tryCatch(simulate_puf(puf_fit$time_hr, ro, model),
                  error = function(e) NULL)
  if (is.null(sim)) return(1e20)
  pred <- predict_at(sim, puf_fit$time_hr)
  if (any(!is.finite(pred))) return(1e20)
  sum((puf_fit$PCB_value - pred)^2)
}

# 4. Fit 0-73 h, predict full period (no late-time refit ------------------
fit_A <- optimize(function(ro) objective(ro, 'A'), ro_bounds)
fit_B <- optimize(function(ro) objective(ro, 'B'), ro_bounds)
results <- data.frame(model = c('B uptake only','A reversible'),
                      ro_cm_day = c(fit_B$minimum, fit_A$minimum),
                      SS_fit = c(fit_B$objective, fit_A$objective))
print(results)
cat('PCB4-scaled predicted ro (not fitted):', ro_reference, 'cm/day\n')
cat('Difference (fitted minus scaled prediction):', fit_B$minimum - ro_reference, 'cm/day\n')

# Use the independently fitted Model B estimate for this congener.
ro_A <- fit_A$minimum
ro_B <- fit_B$minimum

# 4.1 Uncertainty in fitted ro --------------------------------------------
# Number of PUF observations used for fitting
n_obs <- nrow(puf_fit)

# Residual degrees of freedom
df_res <- n_obs - 1

if (df_res <= 0) {
  stop("At least two PUF observations are required.")
}

# Residual variance from Model B
sigma2 <- fit_B$objective / df_res

# Numerical curvature of the objective function
h <- max(1, ro_B * 0.001)

SS_left  <- objective(ro_B - h, "B")
SS_mid   <- objective(ro_B, "B")
SS_right <- objective(ro_B + h, "B")

curvature <- (SS_right - 2 * SS_mid + SS_left) / h^2

if (!is.finite(curvature) || curvature <= 0) {
  stop("Cannot estimate uncertainty: invalid objective curvature.")
}

# Approximate standard error of ro
SE_ro <- sqrt(2 * sigma2 / curvature)

# Approximate 95% confidence interval
t_critical <- qt(0.975, df = df_res)

ro_lower <- ro_B - t_critical * SE_ro
ro_upper <- ro_B + t_critical * SE_ro

# Display results
cat("\nPCB52 - MODEL B UNCERTAINTY\n")
cat("Fitted ro:", round(ro_B, 2), "cm/day\n")
cat("Standard error:", round(SE_ro, 2), "cm/day\n")
cat("95% CI:", round(ro_lower, 2), "to",
    round(ro_upper, 2), "cm/day\n")

Rs_m3_day <- ro_B * Apuf_cm2 / 1e6
cat(sprintf('Model B equivalent air-PUF sampling rate: %.6f m3/day (%.3f L/day)\n',
            Rs_m3_day, Rs_m3_day * 1000))

full_times <- seq(0, max(water_mean$time_hr), length.out = 1000)
sim_A <- simulate_puf(full_times, ro_A, 'A')
sim_B <- simulate_puf(full_times, ro_B, 'B')

# 4.2 Save ro uncertainty results -----------------------------------------
ro_uncertainty <- data.frame(
  congener = congener_id,
  ro_cm_day = ro_B,
  SE_cm_day = SE_ro,
  CI_lower = ro_lower,
  CI_upper = ro_upper)

# 5. Diagnostic summaries -------------------------------------------------
selected_hr <- c(6, 73, 121, 241, 480)
selected_hr <- selected_hr[selected_hr <= max(full_times)]
backtransfer <- data.frame(time_hr = selected_hr,
  Ca = approx(sim_A$time_hr, sim_A$Ca, xout = selected_hr)$y,
  Cpuf = approx(sim_A$time_hr, sim_A$Cpuf, xout = selected_hr)$y) %>%
  mutate(B = ifelse(Ca > 0, Cpuf / (Kpuf * Ca), NA_real_),
         backtransfer_percent = 100 * B)
print(backtransfer)

residuals <- bind_rows(lapply(c('A','B'), function(m) {
  sim <- if (m == 'A') sim_A else sim_B
  data.frame(model = m, time_hr = puf_obs$time_hr,
             observed_ng = puf_obs$PCB_value,
             predicted_ng = predict_at(sim, puf_obs$time_hr))
})) %>% mutate(residual_ng = observed_ng - predicted_ng,
               period = ifelse(time_hr <= fit_end_hr, 'fit', 'prediction'))
print(residuals)

# 6. Plots ----------------------------------------------------------------
# Format fitted parameters
ro_A_label <- round(ro_A, 1)
ro_B_label <- round(ro_B, 1)

# Equivalent air-PUF sampling rate
Rs_B <- ro_B * Apuf_cm2 / 1e6

# 6.1 Water depletion
p_water <- ggplot() +
  geom_point(data = water_obs, aes(time_hr, Cw_ng_L),
             color = "blue", size = 2.5) +
  geom_point(data = water_mean, aes(time_hr, Cw_ng_L),
             color = "red", size = 2.5) +
  geom_line(data = data.frame(time_hr = full_times,
                              Cw = Cw_exp(full_times)),
            aes(time_hr, Cw), color = "black", linewidth = 1) +
  labs(title = "PCB52: Water depletion",
       subtitle = "Piecewise exponential interpolation",
       x = "Time (h)",
       y = "PCB52 concentration in water (ng/L)") +
  theme_bw(base_size = 12)

print(p_water)

# 6.2 Combine Model A and Model B predictions
model_lines <- bind_rows(
  transform(sim_A, model = "Model A: Reversible"),
  transform(sim_B, model = "Model B: Uptake only"))

# 6.3 Model equations for figure annotations
# Model A:
# dCpuf/dt = (ro * Apuf / Vpuf) *
#            (Ca - Cpuf / Kpuf)

# Model B:
# dCpuf/dt = (ro * Apuf / Vpuf) * Ca

equation_A <- paste0(
  "Model A: dC[PUF]/dt == ",
  "(r[o]*A[PUF]/V[PUF])",
  "(C[a]-C[PUF]/K[PUF])")

equation_B <- paste0(
  "Model B: dC[PUF]/dt == ",
  "(r[o]*A[PUF]/V[PUF])*C[a]")

# 6.4 Initial PUF uptake plot
p_initial <- ggplot() +
  geom_point(data = puf_fit, aes(time_hr, PCB_value), color = "blue",
             size = 3) +
  geom_line(data = subset(model_lines, time_hr <= fit_end_hr),
            aes(time_hr, Mpuf_ng, linetype = model), linewidth = 1) +
  annotate("label", x = 50, y = max(puf_fit$PCB_value) * 0.1,
           label = paste0("Model A: ro = ", ro_A_label, " cm/day\n",
                          "Model B: ro = ", ro_B_label, " cm/day\n",
                          "Rs (Model B) = ", round(Rs_B, 5), " m3/day"),
           hjust = 0, size = 3.5, fill = "white") +
  labs(title = "PCB52: Initial PUF uptake kinetics",
       x = "Time (h)", y = "PCB52 mass in PUF (ng/PUF)",
       linetype = "Kinetic model") +
  theme_bw(base_size = 12)

# see plot
print(p_initial)

# save plot
ggsave("Output/Plots/Model/PCB52_InitialUptake.png", plot = p_initial,
       width = 10, height = 5, dpi = 500)

# 6.5 Full-period PUF uptake plot
p_full <- ggplot() +
  geom_point(data = puf_obs, aes(time_hr, PCB_value), color = "blue",
             size = 3) +
  geom_line(data = model_lines, aes(time_hr, Mpuf_ng,
                                    linetype = model), linewidth = 1) +
  geom_vline(xintercept = fit_end_hr, linetype = "dotted",
             color = "gray40") +
  annotate("label", x = 400, y = max(puf_obs$PCB_value) * 0.45,
           label = paste0("Model A: ro = ", ro_A_label, " cm/day\n",
                          "Model B: ro = ", ro_B_label, " cm/day\n",
                          "Rs (Model B) = ", round(Rs_B, 5), " m3/day"),
           size = 3.5, fill = "white") +
  labs(title = "PCB52: Full-period PUF uptake",
       subtitle = paste0("Fitted using 0-73 h | ",
                         "Predicted through 480 h"),
       x = "Time (h)", y = "PCB52 mass in PUF (ng/PUF)",
       linetype = "Kinetic model") +
  theme_bw(base_size = 12)

# see plot
print(p_full)

# 6.6 Display the equations
p_full <- p_full +
  annotate("text", x = 400, y = 80,
           label = paste0("frac(d*C[PUF],dt)==",
                          "frac(r[o]*A[PUF],V[PUF])*",
                          "(C[a]-frac(C[PUF],K[PUF]))"),
           parse = TRUE, size = 4)

# see plot
print(p_full)

p_full <- p_full +
  annotate("text", x = 400, y = 40,
           label = paste0("frac(d*C[PUF],dt)==",
                          "frac(r[o]*A[PUF],V[PUF])*C[a]"),
           parse = TRUE, size = 4)

# see plot
print(p_full)

# save plot
ggsave("Output/Plots/Model/PCB52.png", plot = p_full, width = 10,
       height = 5, dpi = 500)

