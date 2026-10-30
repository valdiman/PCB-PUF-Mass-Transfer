# Review and visualize data

# Install packages
{
  install.packages("tidyr")
  install.packages("dplyr")
  install.packages("ggplot2")
}

# Load libraries
{
  library(tidyr)
  library(dplyr)
  library(ggplot2)
}

# Read data ---------------------------------------------------------------
data_raw <- read.csv("Data/puf_water.csv")
# Check data
str(data_raw)

# Convert data to long format ---------------------------------------------
data_long <- data_raw %>%
  pivot_longer(cols = c(PCB4, PCB18, PCB52),
               names_to = "congener",
               values_to = "PCB_value")
# Check data
str(data_long)

# PUF blank evaluation ----------------------------------------------------
# Separate blanks
blank_data <- data_long %>%
  filter(grepl("^BLK_PUF", sample_id))

# Normal Q-Q plots
p_qq <- ggplot(blank_data, aes(sample = PCB_value)) +
  stat_qq(size = 2.5) +
  stat_qq_line(color = "red") +
  facet_wrap(~ congener, scales = "free") +
  labs(title = "Normal Q-Q plots of PUF blanks", x = "Theoretical quantiles",
       y = "Observed quantiles (ng/PUF)") +
  theme_bw(base_size = 12)

print(p_qq)
# Normal scale looks better

blank_normality <- blank_data %>%
  group_by(congener) %>%
  summarise(n = n(),
            W_original = unname(shapiro.test(PCB_value)$statistic),
            p_original = shapiro.test(PCB_value)$p.value,
            W_log10 = unname(shapiro.test(log10(PCB_value))$statistic),
            p_log10 = shapiro.test(log10(PCB_value))$p.value, .groups = "drop")

print(blank_normality)
# Go with normal scale

# Calculate loq (95% CI)
blank_LOQ <- blank_data %>%
  group_by(congener) %>%
  summarise(n = n(),
            mean_blank = mean(PCB_value),
            sd_blank = sd(PCB_value),
            SE = sd_blank / sqrt(n),
            t_critical = qt(0.975, df = n - 1),
            LOQ = mean_blank + t_critical * SE, .groups = "drop")

# Compare PUF measurements with loq
puf_data <- data_long %>%
  filter(grepl("_PUF_", sample_id),
         !grepl("^BLK_", sample_id))

puf_LOQ <- puf_data %>%
  left_join(blank_LOQ %>% select(congener, LOQ),
    by = "congener"
  ) %>%
  mutate(LOQ_status = if_else(PCB_value >= LOQ,
                              "Above LOQ",
                              "Below LOQ"))

# All values above the loq, no need to correct.

# Plots -------------------------------------------------------------------
# PUF
p.puf <- ggplot(puf_data, aes(x = time_hr, y = PCB_value,
                              color = congener)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "Time (h)", y = "PCB mass (ng/PUF)",
       color = "PCB congener") +
  scale_color_discrete(breaks = c("PCB4", "PCB18", "PCB52")) +
  theme_bw(base_size = 12)

# see plot
p.puf

# save plot
ggsave("Output/Plots/General/Puf_PCB.png", plot = p.puf, width = 10,
       height = 5, dpi = 500)

# Water
water_data <- data_long %>%
  filter(grepl("_W_", sample_id)) %>%
  mutate(Cw_ng_L = PCB_value * 10)

p.water <- ggplot(water_data, aes(x = time_hr, y = Cw_ng_L,
                                  color = congener)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "Time (h)", y = "PCB concentration in water (ng/L)",
       color = "PCB congener") +
  scale_color_discrete(breaks = c("PCB4", "PCB18", "PCB52")) +
  theme_bw(base_size = 12)

# see plot
p.water

# save plot
ggsave("Output/Plots/General/Water_PCB.png", plot = p.water, width = 10,
       height = 5, dpi = 500)

# Export data
write.csv(puf_data,"Data/Processed/puf_data.csv", row.names = FALSE)
write.csv(water_data,"Data/Processed/water_data.csv", row.names = FALSE)

