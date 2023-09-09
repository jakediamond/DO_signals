#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with different parameter sets
# Date: November 19, 2021
# 

# Load libraries
library(furrr) # for running in parallel
library(tidyverse)
library(patchwork)

# Set up parallelization
plan(multisession, workers = 4)

# Load the model and necessary functions
# The model code
source(file.path("src", "1D_DO_model.R"))

# Internal functions for the model
source(file.path("src", "functions_for_model.R"))

# model time step, parameters, and initial conditions
source(file.path("src", "model_initialization.R"))

# Functions that extract change model parameters format results
source(file.path("src", "model_treatment_functions.R"))

# Now get the different treatments ----------------------------------------
# Create dataframe of treatments
# 7 treatments of ramping of GPP
trts <- tibble(
  gpp_choice = 2, # ramped
  ramp_rate = c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1)
)

# Apply 1D DO model to the treatments, "out" is the results of the model
mods <- trts |>
  mutate(
    # "future" runs these in parallel to improve speed
    out = future_pmap(
      list(
        # This is the list of parameters to adjust from the "trts" dataframe
        ramp_rate = ramp_rate,
        gpp_choice = gpp_choice),
      # This is a function that changes the parameters and runs the model
      func_mod))

# Create dataframes for with input data for different reaches
# "mod_df" is the nicely formatted dataframe results of the model
res <- mods |>
  mutate(mod_df = future_map(out, df_fun))

# Save these data
saveRDS(res, file.path("results", "ramped.RDS"))

# Unnest the data so that all the results are in one dataframe
df_res <- select(res, -out) |>
  tidytable::unnest(mod_df)

# Plot over time
# Get legend info for plotting lines
line_info <- tibble(rate = unique(df_res$ramp_rate),
                    lab = c("2-0", "1.67-0.33", "1.33-0.67",
                            "none", "0.67-1.33", "0.33-1.67", "0-2"),
                    color = c("pink", "red", "darkred", "black", "darkblue",
                              "blue", "lightblue"))

# Look at variation in GPP at noon
p_gpp <- ggplot(data = filter(df_res, type == "GPP", time == 36),
                aes(x = dist,
                    y = value,
                    color = as.factor(ramp_rate))) +
  geom_line() +
  scale_color_manual(name = "",
                     values = line_info$color,
                     labels = line_info$lab) +
  theme_classic() +
  scale_x_continuous(breaks = c(0, 10000, 20000)) +
  labs(x = "distance (m)",
       y = expression(GPP~"("*g~O[2]~m^{-2}~d^{-1}*")")) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill='transparent', color=NA))
p_gpp

# Plot over time
p_t <- ggplot() +
  geom_line(data = filter(df_res, type == "DO", reach == 90,
                          between(time, 48, 96)), #middle two days
            aes(x = time - 48, y = value, 
                color = as.factor(ramp_rate)),
                # color = as.factor(round(ramp_rate,2)), group = ramp_rate),
            size = 1.5) +
  scale_color_manual(name = "ramp rate",
                     values = line_info$color,
                     labels = line_info$lab) +
  # scale_color_viridis_d(expression(K~"("*d^{-1}*")")) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  scale_y_continuous(limits = c(6, 14)) +
  theme_classic(base_size = 10) +
  labs(x = "time (hr)",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.5, 0.1),
        legend.direction = "horizontal")
p_t

# Plot the data over space, can pick multiple times
p_d <- ggplot(data = filter(df_res, time %in% c(24, 28, 32, 36, 40, 44),
                     type == "DO",
                     ramp_rate %in% c(-1, 0, 1)),
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(paste0(as.numeric(time %% 24), ":00")),
           group = time)) +
  geom_line() +
  facet_grid(rows = vars(ramp_rate)) +
  scale_color_viridis_d(name = "time") +
  theme_classic() +
  xlab("distance (m)") +
  ylab(expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.45, 0.65),
        legend.background = element_rect(fill = "transparent"),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.2, "cm"))

p_d

layout <- c(
  area(t = 1, l = 1, b = 8, r = 12),
  area(t = 1, l = 6, b = 2, r = 9),
  area(t = 1, l = 13, b = 8, r = 18)
)

p <- p_t + p_gpp + p_d + plot_layout(design = layout)
p

ggsave(plot = p,
       filename = file.path("results", "replicate_fig3.png"),
       dpi = 300,
       units = "cm",
       height = 16,
       width = 18.4)