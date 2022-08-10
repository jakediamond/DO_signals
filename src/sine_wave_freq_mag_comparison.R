#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with estimates from streammetabolizer
# Date: 9 August 2022
# 

# Load libraries
library(furrr)
library(tidyverse)

# Set parallel computing 
plan(multisession, workers = 8)
# Set seed for reproducibility
set.seed(42)

# Load the model
source(file.path("src", "1D_DO_model.R"))

# Load the functions
source(file.path("src", "model_streammetabolizer_comparison_functions.R"))

# First get the integration length (3*v/K) in meters
L <- 3 * (parms[["Q"]] / (parms[["w"]] * parms[["d"]]) * 86400) / #to get from m/s to m/d
  parms[["K"]] #1/d

# For the 1000 iterations version; amp must be less than GPP_max (=2 in base version)
trts <- tibble(amp = runif(1000, min = 0, max = 1),
               phase = runif(1000, min = -1, max = 1),
               f = sample(1:16, 1000, TRUE))

# Look at the sampling space
ggplot(trts,aes(x=amp/parms[["GPP_max"]], y=(parms[["L"]] / f) / L)) + 
  geom_point() +
  labs(y = expression(lambda/L),
       x = expression(sigma))

# For the iterations
mods <- trts %>%
  mutate(out = future_pmap(list(f = f, phase = phase, amp = amp), 
                           func_mod2, .progress = T,
                           .options = furrr_options(seed = TRUE)))

# Convert to usable dataframe
res <- mods %>%
  mutate(data = future_map(out, df_fun, .progress = T))

# Save this data so can skip these steps in the future
# saveRDS(res, file.path("data", "data_sine_freq_mag_test.RDS"))
# res <- readRDS(file.path("data", "data_sine_freq_mag_test.RDS"))

# Get data in nice format for easy plotting
res_p <- select(res, -out) %>%
  mutate(data = map(data, ~ .x %>% 
                             filter(reach == 98) %>%
                             mutate_all(as.character))) %>% 
  tidyfast::dt_unnest(data) %>%
  type_convert()

# get the amplitude and timing of peak DO
res_sine <- res_p %>%
  ungroup() %>%
  filter(type == "DO",
         between(time_hr, 24, 72))  %>%
  mutate(value = as.numeric(value),
         day = ifelse(time_hr < 48, 1, 2),
         hr_day = time_hr %% 24) %>%
  group_by(phase, f, amp, day) %>%
  summarize(range = max(value) - min(value),
            max_time = hr_day[which(value == max(value))]) %>%
  ungroup() %>%
  group_by(phase, f, amp) %>%
  summarize(range = mean(range),
            max_time = mean(max_time))

ggplot(data = res_sine,
       aes(x = 1/f,
           y = phase,
           color = max_time)) +
  geom_point()+
  scale_color_viridis_c() +
  theme_bw()


# Base case model
df_base <- readRDS(file.path("data", "DO_model_base_case.RDS"))

# Calculate the same thing, diel range and time for peak DO
df_base_sum <- df_base %>%
  ungroup() %>%
  filter(type == "DO",
         between(as.numeric(time_hr), 24, 72))  %>%
  mutate(value = as.numeric(value),
         day = ifelse(time_hr < 48, 1, 2),
         hr_day = time_hr %% 24) %>%
  group_by(day) %>%
  summarize(range = max(value) - min(value),
            max_time = hr_day[which(value == max(value))]) %>%
  ungroup() %>%
  summarize(range = mean(range),
            max_time = mean(max_time))


# Get comparison to modeled output
df_p <- res_sine %>%
  mutate(rel_range = range / df_base_sum$range,
         rel_max_time = max_time - df_base_sum$max_time)

# Plot relative change
ggplot(data = df_p,
       aes(x = phase,
           y = amp,
           color = rel_max_time)) +
           # fill = rel_range)) +
  # geom_raster(interpolate = T) +
  geom_point()+
  # stat_density_2d() +
  scale_color_steps2(midpoint = 1, n.breaks = 10) +
  theme_bw()
