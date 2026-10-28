# PCB4 air–PUF uptake model: theory and analysis record

## Purpose
Estimate the air–PUF mass-transfer coefficient, **r₀** (cm/day), from PCB4 uptake in a sealed water–air–PUF experiment. The measured water concentration declines substantially, so the model treats **Cw(t)** as a measured, time-varying input. The main fit uses the initial ~0–73 h observations; 73–480 h is a conditional prediction period.

## Experimental inputs and units

| Quantity | Value | Units / basis |
|---|---:|---|
| Water volume | 200 | cm³ |
| Modeled air volume | 135.98 | cm³, lower region between water and PUF |
| PUF volume | 29 | cm³ |
| Water–air interface area | 56.74502 | cm² |
| PUF transfer area | 7.07 | cm² |
| Initial measured PCB4 water concentration | 3204.195 | ng/L |
| Air–water partition coefficient, Kaw | 10^(-1.87155) | Ca/Cw, dimensionless |
| PUF–air partition coefficient, Kpuf | 10^(0.6366 × 6.789298102 − 3.1774) × 21300 | Cpuf/Ca, dimensionless (empirical estimate) |
| Fixed overall water–air transfer coefficient | 52.82 | cm/day; provisional weak-stirring scenario |
| Temperature and orbital shaking | 298.15; 60 | K; rpm |

**Geometry qualification:** The modeled 135.98 cm³ is the lower gas region, not necessarily the entire sealed flask gas volume. The nominal flask capacity and liquid fill should be checked independently; this model uses the previously agreed compartment approximation.

## Governing equations

Let **Ca** and **Cpuf** be concentrations in ng/L of air and ng/L of PUF volume, respectively; time **t** is in days. Water concentration **Cw(t)** is supplied from measurements. With transfer velocities in cm/day, areas in cm² and volumes in cm³:

**Model A: reversible exchange**

    dCa/dt   = (kaw,o Aaw/Va) [Cw(t) − Ca/Kaw] − (r₀ Apuf/Va) [Ca − Cpuf/Kpuf]
    dCpuf/dt = (r₀ Apuf/Vpuf) [Ca − Cpuf/Kpuf]

**Model B: initial uptake (primary model)**

    dCa/dt   = (kaw,o Aaw/Va) [Cw(t) − Ca/Kaw] − (r₀ Apuf/Va) Ca
    dCpuf/dt = (r₀ Apuf/Vpuf) Ca

Both models retain transient headspace dynamics: **neither assumes Ca = Kaw Cw nor dCa/dt = 0**. Initial conditions are Ca(0) = 0 and Cpuf(0) = 0. Predicted mass is **Mpuf = Cpuf × Vpuf / 1000** (ng).

## Water concentration interpolation

Water measurements are averaged within experimental sampling groups; the independently measured time-zero concentration is inserted. Between consecutive mean concentrations, the preferred interpolation is **linear in log concentration**, equivalent to a piecewise exponential:

    Cw(t) = Cw(ti) exp{ [ln(Cw(ti+1)/Cw(ti)) / (ti+1−ti)] (t−ti) }.

This passes through the measured means and remains positive. It is **not** a fitted global first-order depletion law. Piecewise linear interpolation was also tested as a sensitivity case. The time-series mean is used as an external forcing function; this model does not reconstruct the full closed-reactor mass balance.

## Parameter fitting and results

The script minimizes the **unweighted sum of squared residuals** between measured PUF mass (ng/PUF) and predicted Mpuf for observations through ~73 h, using one fitted parameter, r₀. The water–air coefficient is fixed at 52.82 cm/day. Fitting is performed with `optimize()`; late-time PUF observations are excluded from parameter estimation.

Previous results from the analysis:

| Water interpolation | Model | Estimated r₀ (cm/day) | Fit SS (ng²) |
|---|---|---:|---:|
| Piecewise linear | B, uptake only | 1369.140 | 5157.285 |
| Piecewise exponential | B, uptake only | **1619.601** | **4671.426** |

The exponential estimate is ~18.3% larger; SS is ~9.4% lower. The updated script refits the exponential case and prints the difference from 1619.601 cm/day as a reproducibility check. Small numerical differences may occur with solver settings and interpolation knots.

## Equivalent sampling rate

The **air–PUF exchange conductance**, expressed as a volume/time, is:

    Rs = r₀ × Apuf / 10^6   [m³/day].

At r₀ = 1619.601 cm/day and Apuf = 7.07 cm², **Rs = 0.01145 m³/day = 11.45 L/day**. This is an air–PUF equivalent sampling rate, **not** the overall water-to-PUF sampling rate.

## Interpretation and limitations

The early reversible-model diagnostic gave a back-transfer/forward-term ratio of about 1.1% at 73 h, supporting uptake-only Model B in the fitting window. In a previous full-period reversible simulation, this ratio became much larger later (about 12% at 121 h, 53% at 241 h and 34% at 480 h), so **Model B should not be interpreted as evidence of negligible late-time back-transfer**. Both formulations can generate an apparent plateau when measured water concentrations decline; a plateau alone does not prove thermodynamic equilibrium.

The fitted r₀ is **conditional** on the assumed 52.82 cm/day water–air transfer coefficient, PUF geometry, initial conditions, and interpolated water trajectory. The water–air coefficient comes from a weak-stirring literature-based scenario, not a direct calibration of this 60-rpm orbital flask. The original water/PUF mass-balance discrepancy is unresolved, and measured headspace concentrations are unavailable. The model is therefore a practical uptake characterization rather than a complete mechanistic closure of all PCB mass.

## Script organization

`PCB4_PUF_uptake_updated.R` loads data, constructs the full-period piecewise exponential water input, fits r₀ on early PUF measurements for Models A and B, predicts the complete period, calculates back-transfer and residual diagnostics, and plots water depletion, early PUF uptake and full-period uptake. All user-editable physical parameters appear near the top.

## Source context

Kanwisher (1963), *On the Exchange of Gases Between the Atmosphere and the Sea*, provided diffusion-film measurements used in earlier discussions of the weak-stirring water-side transfer scenario. The fixed overall 52.82 cm/day value is a modeling assumption derived in our preceding comparison, not a directly reported coefficient from Kanwisher.
