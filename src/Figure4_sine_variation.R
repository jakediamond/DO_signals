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
# 6 treatments of varying GPP frequency and phase, last is normal
trts <- tibble(
  gpp_choice = c(rep(3, 6), 1), #sine then constant
  f = c(1, 1, 1, 2, 4, 8, 0),
  phase = c(0, 0.25, -0.25, 0, 0, 0, 0)
)

# Apply 1D DO model to the treatments
mods <- trts |>
  mutate(
    # "future" runs these in parallel to improve speed
    out = future_pmap(
      list(
        # This is the list of parameters to adjust from the "trts" dataframe
        gpp_choice = gpp_choice,
        f = f,
        phase = phase),
      # This is a function that changes the parameters and runs the model
      func_mod))

# Create dataframe for with input data for different reaches
# "mod_df" is the nicely formatted dataframe results of the model
res <- mods |>
  mutate(mod_df = future_map(out, df_fun))

# Save these data
saveRDS(res, file.path("results", "sine.RDS"))

# Clean the results up for plotting
df_res <- select(res, -out) |>
  tidytable::unnest(mod_df) |>
  mutate(line = if_else(phase != 0,
                        paste0("freq ", f, ", ", "phase ", phase),
                        paste("freq", f)),
         line = if_else(f == 0, "uniform", line))

# Plot over time
# Get legend info for plotting lines
line_info <- tibble(line = unique(df_res$line),
                    color = c(rep("darkred", 3), "darkblue", "blue", "lightblue", "black"),
                    linetype = c(1, 2, 3, rep(1,4)))

# Look at variation in GPP at noon
p_gpp <- ggplot(data = filter(df_res, type == "GPP", time == 36,
                     line != "uniform") |>
         mutate(fac = if_else(f == 1, "b", "a")),
       aes(x = dist,
           y = value,
           color = line,
           linetype = line)) +
  geom_line() +
  scale_color_manual(name = "",
                     values = line_info$color,
                     labels = line_info$line) +
  scale_linetype_manual(name = "",
                        values = line_info$linetype,
                        labels = line_info$line) +
  facet_grid(rows = vars(fac)) +
  theme_classic() +
  labs(x = "distance (m)",
       y = expression(GPP~"("*g~O[2]~m^{-2}~d^{-1}*")")) +
  theme(legend.position = "none",
        strip.text = element_blank())
p_gpp

# Plot over time
p_t <- ggplot() +
  geom_line(data = filter(df_res, type == "DO", reach == 95),
            aes(x = time, y = value, 
                color = line, 
                linetype = line,
                group = line),
            size = 1) +
  scale_color_manual(name = "",
                     values = line_info$color,
                     labels = line_info$line) +
  scale_linetype_manual(name = "",
                        values = line_info$linetype,
                        labels = line_info$line) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "time (hr)",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.15, 0.8),
        legend.background = element_rect(color = "transparent",
                                         fill = "transparent"))

p_t

# Plot the data over space, can pick multiple times
p_d <- ggplot(data = filter(df_res, time %in% c(24, 28, 32, 36, 40, 44),
                     type == "DO", line != "uniform"),
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(paste0(as.numeric(time %% 24), ":00")),
           group = time)) +
  geom_line() +
  facet_grid(rows = vars(line)) +
  scale_color_viridis_d(name = "") +
  # scale_y_continuous()
  theme_classic() +
  xlab("distance (m)") +
  ylab(expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.5, 0.15),
        legend.background = element_rect(fill = "transparent"),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.2, "cm"))
  # guides(colour = guide_legend(nrow = 1))
p_d

layout <- c(
  area(t = 1, l = 1, b = 2, r = 4),
  area(t = 3, l = 1, b = 7, r = 4),
  area(t = 1, l = 5, b = 7, r = 6)
)

p <- p_gpp + p_t + p_d + plot_layout(design = layout)
p

ggsave(plot = p,
       filename = file.path("results", "replicate_fig4.png"),
       dpi = 300,
       units = "cm",
       height = 16,
       width = 18.4)
