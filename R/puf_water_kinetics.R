# ============================================================
# PCB-PUF MASS TRANSFER
# WATER-AIR-PUF KINETIC EXPERIMENT
#
# STEP 1: DATA IMPORT AND PREPARATION
# ============================================================


# ------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# ------------------------------------------------------------
# 2. IMPORT RAW DATA
# ------------------------------------------------------------

data_raw <- read.csv(
  "Data/puf_water.csv",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Inspect data
str(data_raw)
head(data_raw)
names(data_raw)
nrow(data_raw)


# ------------------------------------------------------------
# 3. EXPERIMENTAL PARAMETERS
# ------------------------------------------------------------

# Water volume (L)
Vw <- 0.200

# Water volume extracted for analysis (L)
V_extracted <- 0.100

# Headspace geometry (cm)
h <- 4.8
D_top <- 3.1
D_bottom <- 8.5

# Calculate headspace volume (cm3)
Va_cm3 <- (pi * h / 12) *
  (D_top^2 + D_top * D_bottom + D_bottom^2)

# Convert to liters
Va <- Va_cm3 / 1000

# PUF volume (L)
Vpuf <- 0.029

# Shaking speed (rpm)
shaking_rpm <- 60

# Display experimental parameters
cat("\nEXPERIMENTAL PARAMETERS\n")
cat("-----------------------\n")
cat("Water volume:", Vw, "L\n")
cat("Water extracted:", V_extracted, "L\n")
cat("Headspace volume:", round(Va, 4), "L\n")
cat("PUF volume:", Vpuf, "L\n")
cat("Shaking speed:", shaking_rpm, "rpm\n")


# ------------------------------------------------------------
# 4. CONVERT PCB DATA TO LONG FORMAT
# ------------------------------------------------------------

data_long <- data_raw %>%
  pivot_longer(
    cols = c(PCB4, PCB18.30, PCB52),
    names_to = "congener",
    values_to = "PCB_value"
  )


# ------------------------------------------------------------
# 5. IDENTIFY SAMPLE TYPES AND REPLICATES
# ------------------------------------------------------------

data_long <- data_long %>%
  mutate(
    
    sample_type = case_when(
      str_detect(sample_id, "^BLK_PUF") ~ "Blank",
      str_detect(sample_id, "_PUF_") ~ "PUF",
      str_detect(sample_id, "_W_") ~ "Water",
      TRUE ~ "Unknown"
    ),
    
    replicate = str_extract(
      sample_id,
      "R[0-9]+$"
    )
    
  )


# ------------------------------------------------------------
# 6. ASSIGN NOMINAL SAMPLING TIMES
# ------------------------------------------------------------

# Actual time_hr values are preserved.
# Nominal times are used only to organize replicates.

data_long <- data_long %>%
  mutate(
    
    time_nominal = case_when(
      is.na(time_hr) ~ NA_real_,
      time_hr < 1   ~ 0.5,
      time_hr < 12  ~ 6,
      time_hr < 100 ~ 72,
      time_hr < 180 ~ 120,
      time_hr < 350 ~ 240,
      TRUE ~ 480
    )
    
  )


# ------------------------------------------------------------
# 7. SEPARATE WATER, PUF AND BLANK DATA
# ------------------------------------------------------------

water_data <- data_long %>%
  filter(sample_type == "Water")

puf_data <- data_long %>%
  filter(sample_type == "PUF")

blank_data <- data_long %>%
  filter(sample_type == "Blank")


# ------------------------------------------------------------
# 8. CALCULATE WATER CONCENTRATION
# ------------------------------------------------------------

# Original water data: ng/100 mL
#
# Cw (ng/L) = PCB_value / 0.100 L
#
# Total water PCB mass:
# Mw (ng) = Cw * 0.200 L

water_data <- water_data %>%
  mutate(
    
    Cw_ng_L = PCB_value / V_extracted,
    
    Mw_ng = Cw_ng_L * Vw
    
  )


# ------------------------------------------------------------
# 9. ORGANIZE PUF MEASUREMENTS
# ------------------------------------------------------------

# Original PUF data: ng/PUF
#
# No blank correction at this stage.

puf_data <- puf_data %>%
  mutate(
    M_PUF_ng = PCB_value
  )


# ------------------------------------------------------------
# 10. ORGANIZE PUF BLANKS
# ------------------------------------------------------------

blank_data <- blank_data %>%
  mutate(
    M_blank_ng = PCB_value
  )


# ------------------------------------------------------------
# 11. CHECK DATA STRUCTURE
# ------------------------------------------------------------

cat("\nSAMPLE COUNTS\n")
cat("-------------\n")

print(
  data_long %>%
    count(sample_type)
)

cat("\nSAMPLING TIMES\n")
cat("--------------\n")

print(
  data_long %>%
    filter(sample_type != "Blank") %>%
    distinct(time_nominal, time_hr) %>%
    arrange(time_nominal, time_hr)
)


# ------------------------------------------------------------
# 12. COMBINE WATER AND PUF MEASUREMENTS
# ------------------------------------------------------------

# IMPORTANT:
# This pairing assumes that water and PUF samples
# with the same replicate ID came from the same reactor.
#
# If they came from different reactors, this pairing
# must not be used for replicate-level mass balances.

mass_balance_data <- water_data %>%
  select(
    congener,
    time_nominal,
    replicate,
    time_water_hr = time_hr,
    Cw_ng_L,
    Mw_ng
  ) %>%
  left_join(
    
    puf_data %>%
      select(
        congener,
        time_nominal,
        replicate,
        time_puf_hr = time_hr,
        M_PUF_ng
      ),
    
    by = c(
      "congener",
      "time_nominal",
      "replicate"
    )
    
  ) %>%
  mutate(
    M_observed_ng = Mw_ng + M_PUF_ng
  )


# ------------------------------------------------------------
# 13. CHECK THE COMBINED DATA
# ------------------------------------------------------------

cat("\nWATER-PUF DATA\n")
cat("--------------\n")

print(mass_balance_data)


# ------------------------------------------------------------
# 14. PREPARE KINETIC DATA
# ------------------------------------------------------------

# Proposed kinetic fitting window:
# approximately 0.5 to 200 hours.
#
# Actual measured points within this interval:
# 0.5, 6, 73 and 121 hours.

water_kinetic <- water_data %>%
  filter(
    time_hr >= 0.4,
    time_hr <= 200
  )

puf_kinetic <- puf_data %>%
  filter(
    time_hr >= 0.4,
    time_hr <= 200
  )

library(deSolve)

water_air_puf_model <- function(t, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    # Concentrations (ng/L)
    Cw <- Mw / Vw
    Ca <- Ma / Va
    
    # Water-air transfer (ng/h)
    Jwa <- Kwa * (Cw - Ca / Kaw)
    
    # Air-PUF transfer (ng/h)
    Jap <- kua * Ca - ke * Mpuf
    
    # Mass balances (ng/h)
    dMw <- -Jwa
    dMa <- Jwa - Jap
    dMpuf <- Jap
    
    list(
      c(dMw, dMa, dMpuf),
      Cw = Cw,
      Ca = Ca,
      Jwa = Jwa,
      Jap = Jap
    )
  })
}

# ------------------------------------------------------------
# STEP 4A: MODEL SIMULATION
# ------------------------------------------------------------

simulate_pcb <- function(
    times,
    Cw0,
    Kwa,
    kua,
    Kaw,
    ke = 0
) {
  
  # Initial conditions (ng)
  initial_state <- c(
    Mw = Cw0 * Vw,
    Ma = 0,
    Mpuf = 0
  )
  
  # Model parameters
  parameters <- c(
    Vw = Vw,
    Va = Va,
    Kwa = Kwa,
    Kaw = Kaw,
    kua = kua,
    ke = ke
  )
  
  # Simulation times
  times <- sort(unique(c(0, times)))
  
  # Solve differential equations
  simulation <- deSolve::ode(
    y = initial_state,
    times = times,
    func = water_air_puf_model,
    parms = parameters,
    method = "lsoda"
  )
  
  # Convert to data frame and calculate concentrations
  simulation <- as.data.frame(simulation) %>%
    mutate(
      Cw_pred = Mw / Vw,
      Ca_pred = Ma / Va,
      Mpuf_pred = Mpuf
    )
  
  return(simulation)
}

