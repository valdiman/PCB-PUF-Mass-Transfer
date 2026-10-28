# PCB52 PUF uptake | clean reproducible analysis
# Run from project root containing Data/Processed/*.csv
# Model B is primary; Model A is a reversible diagnostic.

library(dplyr)
library(ggplot2)
library(deSolve)

# ---- 1. Inputs (edit here) ----
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
# For a combined PCB18+30 peak, the CSV label may differ from 'PCB18.30'.
if (!(congener_id %in% water_data$congener) ||
    !(congener_id %in% puf_data$congener)) {
  stop(paste0('Congener label not found: ', congener_id,
              '\nWater labels: ', paste(unique(water_data$congener), collapse=', '),
              '\nPUF labels: ', paste(unique(puf_data$congener), collapse=', ')))
}

# ---- 2. Water data: full-period, grouped by sampling period ----
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

# ---- 3. Full transient air + PUF model ----
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

# ---- 4. Fit 0-73 h, predict full period (no late-time refit) ----
fit_B <- optimize(function(ro) objective(ro, 'B'), ro_bounds)
fit_A <- optimize(function(ro) objective(ro, 'A'), ro_bounds)
results <- data.frame(model = c('B uptake only','A reversible'),
                      ro_cm_day = c(fit_B$minimum, fit_A$minimum),
                      SS_fit = c(fit_B$objective, fit_A$objective))
print(results)
cat('PCB4-scaled predicted ro (not fitted):', ro_reference, 'cm/day\n')
cat('Difference (fitted minus scaled prediction):', fit_B$minimum - ro_reference, 'cm/day\n')

# Use the independently fitted Model B estimate for this congener.
ro_B <- fit_B$minimum
ro_A <- fit_A$minimum
Rs_m3_day <- ro_B * Apuf_cm2 / 1e6
cat(sprintf('Model B equivalent air-PUF sampling rate: %.6f m3/day (%.3f L/day)\n',
            Rs_m3_day, Rs_m3_day * 1000))

full_times <- seq(0, max(water_mean$time_hr), length.out = 1000)
sim_B <- simulate_puf(full_times, ro_B, 'B')
sim_A <- simulate_puf(full_times, ro_A, 'A')

# ---- 5. Diagnostic summaries ----
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

# ---- 6. Figures ----
p_water <- ggplot() +
  geom_point(data = water_obs, aes(time_hr, Cw_ng_L), color = 'blue') +
  geom_point(data = water_mean, aes(time_hr, Cw_ng_L), color = 'red') +
  geom_line(data = data.frame(time_hr = full_times, Cw = Cw_exp(full_times)),
            aes(time_hr, Cw), linewidth = 0.9) +
  labs(title = 'PCB52: measured water depletion',
       subtitle = 'Blue: observations | Red: means | Black: piecewise exponential',
       x = 'Time (h)', y = 'PCB52 in water (ng/L)') + theme_bw()

model_lines <- bind_rows(transform(sim_A, model = 'A: reversible'),
                         transform(sim_B, model = 'B: uptake only'))
p_initial <- ggplot() +
  geom_point(data = puf_fit, aes(time_hr, PCB_value), color = 'blue', size = 2.5) +
  geom_line(data = subset(model_lines, time_hr <= fit_end_hr),
            aes(time_hr, Mpuf_ng, linetype = model), linewidth = 0.9) +
  labs(title = 'PCB52: initial PUF uptake (fit)', x = 'Time (h)',
       y = 'PCB52 in PUF (ng/PUF)') + theme_bw()
p_full <- ggplot() +
  geom_point(data = puf_obs, aes(time_hr, PCB_value), color = 'blue', size = 2.5) +
  geom_line(data = model_lines, aes(time_hr, Mpuf_ng, linetype = model), linewidth = 0.9) +
  geom_vline(xintercept = fit_end_hr, linetype = 'dotted') +
  labs(title = 'PCB52: full-period PUF uptake',
       subtitle = 'Coefficients fitted only through ~73 h; water is measured input',
       x = 'Time (h)', y = 'PCB52 in PUF (ng/PUF)') + theme_bw()
print(p_water); print(p_initial); print(p_full)

print(results)
