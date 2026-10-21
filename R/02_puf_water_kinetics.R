# ============================================================
# PCB-PUF MASS TRANSFER
# WATER-AIR-PUF KINETIC MODEL
#
# SCRIPT 02: MODELING
# ============================================================


# Load libraries ----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(deSolve)


# Read data ---------------------------------------------------------------

puf_data <- read.csv("Data/Processed/puf_data.csv")
water_data <- read.csv("Data/Processed/water_data.csv")


# Experimental parameters -------------------------------------------------

# Compartment volumes (cm3)

Vw <- 200
Vpuf <- 29

# Headspace geometry (cm)

h <- 4.8
D_top <- 3.1
D_bottom <- 8.5

# Headspace volume (cm3)

Va <- (pi * h / 12) *
  (D_top^2 + D_top * D_bottom + D_bottom^2)

# Interfacial areas (cm2)

Aaw <- pi * (D_bottom / 2)^2
Apuf <- 7.07

# Experimental conditions

T <- 298.15
shaking_rpm <- 60

# Display geometry

cat("Water volume:", Vw, "cm3\n")
cat("Air volume:", Va, "cm3\n")
cat("PUF volume:", Vpuf, "cm3\n")
cat("Water-air area:", Aaw, "cm2\n")
cat("PUF area:", Apuf, "cm2\n")


# Air-water partition coefficients ---------------------------------------

partition_data <- data.frame(
  congener = c("PCB4", "PCB18.30", "PCB52"),
  logKaw = c(-1.87155, -1.88355, -1.88455)
)

partition_data$Kaw <- 10^partition_data$logKaw

print(partition_data)


# PUF-air partition coefficient ------------------------------------------

# PCB4 octanol-air partition coefficient

logKoa_PCB4 <- 6.789298102

# PUF density (g/m3)

rho_puf <- 21300

# PCB4 PUF-air partition coefficient
# Dimensionless after multiplication by PUF density

Kpuf_PCB4 <- 10^(
  0.6366 * logKoa_PCB4 - 3.1774
) * rho_puf

cat("PCB4 Kpuf:", Kpuf_PCB4, "\n")


# Measured initial water concentrations ----------------------------------

initial_concentrations <- data.frame(
  congener = c("PCB4", "PCB18.30", "PCB52"),
  Cw0_ng_L = c(
    3204.195,
    3190.523,
    4538.160
  )
)

# Initial PCB mass in water (ng)

initial_concentrations <- initial_concentrations %>%
  mutate(
    Mw0_ng = Cw0_ng_L * Vw / 1000
  )

print(initial_concentrations)


# Prepare experimental observations --------------------------------------

# Convert time from hours to days

water_data <- water_data %>%
  mutate(
    time_day = time_hr / 24
  )

puf_data <- puf_data %>%
  mutate(
    time_day = time_hr / 24
  )

# Water observations

water_obs <- water_data %>%
  transmute(
    congener,
    time_day,
    sample_type = "Water",
    observed = Cw_ng_L
  )

# PUF observations

puf_obs <- puf_data %>%
  transmute(
    congener,
    time_day,
    sample_type = "PUF",
    observed = PCB_value
  )

# Combine observations

observations <- bind_rows(
  water_obs,
  puf_obs
)

print(head(observations))


# ============================================================
# WATER-AIR-PUF KINETIC MODEL
# ============================================================

# State variables:
#
# Cw   = PCB concentration in water (ng/L)
# Ca   = PCB concentration in air (ng/L)
# Cpuf = PCB concentration in PUF (ng/L PUF)
#
# Transfer coefficients:
#
# kaw_o = water-air mass-transfer coefficient (cm/day)
# ro    = air-PUF mass-transfer coefficient (cm/day)
#
# Partition coefficients:
#
# Kaw  = air-water partition coefficient
# Kpuf = PUF-air partition coefficient


water_air_puf_model <- function(t, state, parameters) {
  
  # State variables
  
  Cw <- unname(state["Cw"])
  Ca <- unname(state["Ca"])
  Cpuf <- unname(state["Cpuf"])
  
  # Parameters
  
  Vw <- unname(parameters["Vw"])
  Va <- unname(parameters["Va"])
  Vpuf <- unname(parameters["Vpuf"])
  
  Aaw <- unname(parameters["Aaw"])
  Apuf <- unname(parameters["Apuf"])
  
  Kaw <- unname(parameters["Kaw"])
  Kpuf <- unname(parameters["Kpuf"])
  
  kaw_o <- unname(parameters["kaw_o"])
  ro <- unname(parameters["ro"])
  
  # Differential equations (ng/L/day)
  
  dCw <- -(kaw_o * Aaw / Vw) *
    (Cw - Ca / Kaw)
  
  dCa <- (kaw_o * Aaw / Va) *
    (Cw - Ca / Kaw) -
    (ro * Apuf / Va) *
    (Ca - Cpuf / Kpuf)
  
  dCpuf <- (ro * Apuf / Vpuf) *
    (Ca - Cpuf / Kpuf)
  
  return(
    list(c(dCw, dCa, dCpuf))
  )
}


# ============================================================
# MODEL SIMULATION
# ============================================================

simulate_pcb <- function(times,
                         Cw0,
                         kaw_o,
                         ro,
                         Kaw,
                         Kpuf) {
  
  # Initial conditions (ng/L)
  
  initial_state <- c(
    Cw = unname(Cw0),
    Ca = 0,
    Cpuf = 0
  )
  
  # Parameters
  
  parameters <- c(
    Vw = unname(Vw),
    Va = unname(Va),
    Vpuf = unname(Vpuf),
    Aaw = unname(Aaw),
    Apuf = unname(Apuf),
    Kaw = unname(Kaw),
    Kpuf = unname(Kpuf),
    kaw_o = unname(kaw_o),
    ro = unname(ro)
  )
  
  # Simulation times (days)
  
  times <- sort(unique(c(0, times)))
  
  # Solve differential equations
  
  simulation <- tryCatch(
    suppressWarnings(
      deSolve::ode(
        y = initial_state,
        times = times,
        func = water_air_puf_model,
        parms = parameters,
        method = "lsoda"
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(simulation)) {
    return(NULL)
  }
  
  simulation <- as.data.frame(simulation)
  
  # Verify simulation
  
  required_columns <- c(
    "time", "Cw", "Ca", "Cpuf"
  )
  
  if (!all(required_columns %in% names(simulation))) {
    return(NULL)
  }
  
  if (nrow(simulation) != length(times)) {
    return(NULL)
  }
  
  if (!all(is.finite(as.matrix(simulation)))) {
    return(NULL)
  }
  
  # Calculate compartment masses (ng)
  
  simulation <- simulation %>%
    mutate(
      Mw_ng = Cw * Vw / 1000,
      Ma_ng = Ca * Va / 1000,
      Mpuf_ng = Cpuf * Vpuf / 1000,
      
      Mtotal_ng = Mw_ng + Ma_ng + Mpuf_ng
    )
  
  return(simulation)
}


# ============================================================
# TEST MODEL SIMULATION: PCB4
# ============================================================

# Partition coefficients

Kaw_PCB4 <- partition_data$Kaw[
  partition_data$congener == "PCB4"
]

# Measured initial water concentration

Cw0_PCB4 <- initial_concentrations$Cw0_ng_L[
  initial_concentrations$congener == "PCB4"
]

# Experimental sampling times (days)

times_test <- c(
  0, 0.5, 6, 73, 121, 241, 480
) / 24

# Hypothetical transfer coefficients (cm/day)

kaw_o_test <- 10
ro_test <- 10

# Test simulation

test_simulation <- simulate_pcb(
  times = times_test,
  Cw0 = Cw0_PCB4,
  kaw_o = kaw_o_test,
  ro = ro_test,
  Kaw = Kaw_PCB4,
  Kpuf = Kpuf_PCB4
)

print(test_simulation)


# Mass conservation -------------------------------------------------------

if (!is.null(test_simulation)) {
  
  mass_error <- max(
    abs(
      test_simulation$Mtotal_ng -
        test_simulation$Mtotal_ng[1]
    )
  )
  
  cat("Maximum mass-balance error:", mass_error, "ng\n")
}


# ============================================================
# OBJECTIVE FUNCTION
# ============================================================

# Estimate:
#
# kaw_o = water-air transfer coefficient (cm/day)
# ro    = air-PUF transfer coefficient (cm/day)
#
# Initial concentration is fixed at the measured value.


objective_function <- function(log_parameters,
                               congener_name,
                               Kaw_value,
                               Kpuf_value,
                               Cw0_value) {
  
  # Convert parameters to original scale
  
  pars <- exp(log_parameters)
  
  kaw_o <- unname(pars[1])
  ro <- unname(pars[2])
  
  # Select observations
  
  obs <- observations %>%
    filter(congener == congener_name)
  
  # Run simulation
  
  sim <- simulate_pcb(
    times = obs$time_day,
    Cw0 = Cw0_value,
    kaw_o = kaw_o,
    ro = ro,
    Kaw = Kaw_value,
    Kpuf = Kpuf_value
  )
  
  # Penalize failed simulations
  
  if (is.null(sim)) {
    return(1e20)
  }
  
  # Match predictions with observations
  
  idx <- match(obs$time_day, sim$time)
  
  if (any(is.na(idx))) {
    return(1e20)
  }
  
  predicted <- ifelse(
    obs$sample_type == "Water",
    sim$Cw[idx],
    sim$Mpuf_ng[idx]
  )
  
  # Separate scaling factors
  
  water_scale <- median(
    obs$observed[obs$sample_type == "Water"],
    na.rm = TRUE
  )
  
  puf_scale <- median(
    obs$observed[obs$sample_type == "PUF"],
    na.rm = TRUE
  )
  
  scales <- ifelse(
    obs$sample_type == "Water",
    water_scale,
    puf_scale
  )
  
  if (any(!is.finite(predicted)) ||
      any(!is.finite(scales)) ||
      any(scales <= 0)) {
    return(1e20)
  }
  
  # Weighted sum of squared residuals
  
  sum(
    ((obs$observed - predicted) / scales)^2
  )
}


# ============================================================
# FIT PCB4
# ============================================================

# Initial parameter guesses (cm/day)

initial_guess <- c(
  kaw_o = 10,
  ro = 10
)

# Provisional parameter bounds (cm/day)
# These are numerical search limits, not physical constraints.

lower <- c(
  kaw_o = 1e-4,
  ro = 1e-4
)

upper <- c(
  kaw_o = 1000,
  ro = 1000
)

# Optimization

fit_PCB4 <- optim(
  par = log(initial_guess),
  fn = objective_function,
  
  congener_name = "PCB4",
  Kaw_value = Kaw_PCB4,
  Kpuf_value = Kpuf_PCB4,
  Cw0_value = Cw0_PCB4,
  
  method = "L-BFGS-B",
  
  lower = log(lower),
  upper = log(upper),
  
  control = list(maxit = 10000)
)

# Estimated parameters

estimated_parameters <- setNames(
  exp(fit_PCB4$par),
  names(initial_guess)
)

cat("\nESTIMATED PARAMETERS: PCB4\n")

print(estimated_parameters)

cat("\nOptimization convergence:\n")
print(fit_PCB4$convergence)

cat("\nObjective function value:\n")
print(fit_PCB4$value)


# ============================================================
# MODEL PREDICTIONS
# ============================================================

# Simulation times (days)

times_plot <- seq(
  0,
  480 / 24,
  length.out = 500
)

# Simulate fitted model

sim_PCB4 <- simulate_pcb(
  times = times_plot,
  Cw0 = Cw0_PCB4,
  kaw_o = unname(estimated_parameters["kaw_o"]),
  ro = unname(estimated_parameters["ro"]),
  Kaw = Kaw_PCB4,
  Kpuf = Kpuf_PCB4
)

# Convert predicted time to hours for plotting

if (!is.null(sim_PCB4)) {
  
  sim_PCB4 <- sim_PCB4 %>%
    mutate(
      time_hr = time * 24
    )
}


# ============================================================
# OBSERVED VS PREDICTED: WATER
# ============================================================

p_water_fit <- ggplot() +
  
  geom_point(
    data = water_data %>%
      filter(congener == "PCB4"),
    
    aes(
      x = time_hr,
      y = Cw_ng_L
    ),
    
    size = 3,
    color = "blue"
  ) +
  
  geom_line(
    data = sim_PCB4,
    
    aes(
      x = time_hr,
      y = Cw
    ),
    
    linewidth = 1,
    color = "black"
  ) +
  
  labs(
    x = "Time (h)",
    y = "PCB4 concentration in water (ng/L)",
    title = "PCB4: Water depletion"
  ) +
  
  theme_bw(base_size = 12)

print(p_water_fit)


# ============================================================
# OBSERVED VS PREDICTED: PUF
# ============================================================

p_puf_fit <- ggplot() +
  
  geom_point(
    data = puf_data %>%
      filter(congener == "PCB4"),
    
    aes(
      x = time_hr,
      y = PCB_value
    ),
    
    size = 3,
    color = "blue"
  ) +
  
  geom_line(
    data = sim_PCB4,
    
    aes(
      x = time_hr,
      y = Mpuf_ng
    ),
    
    linewidth = 1,
    color = "black"
  ) +
  
  labs(
    x = "Time (h)",
    y = "PCB4 mass in PUF (ng/PUF)",
    title = "PCB4: PUF accumulation"
  ) +
  
  theme_bw(base_size = 12)

print(p_puf_fit)


# ============================================================
# END OF SCRIPT 02
# ============================================================