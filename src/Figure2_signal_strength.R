#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with different parameter sets
# Date: November 19, 2021
# 

# Load libraries
library(furrr) # for running in parallel
library(tidyverse)

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
# 12 treatments of varying K, D, and storage area, final is plug flow
trts <- tibble(
  K = c(1, 5, 20, rep(2, 9)),
  D = c(4, 4, 4, 2, 8, 32, 64, 4, 4, 4, 4, 0) * 3600,
  A_stor_frac = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.01, 0.1, 0.5, 2, 0.2),
  alpha = c(rep(0.252, 11), 0)
  )

# Apply 1D DO model to the treatments
mods <- trts |>
  mutate(
    # "future" runs these in parallel to improve speed
    out = future_pmap(
      list(
        # This is the list of parameters to adjust from the "trts" dataframe
        K = K,
        D = D,
        A_stor_frac = A_stor_frac),
      # This is a function that changes the parameters and runs the model
      func_mod))

# Create dataframe for with input data for different reaches
# "mod_df" is the nicely formatted dataframe results of the model
res <- mods |>
  mutate(mod_df = future_map(out, df_fun))

# Save these data
saveRDS(res, file.path("results", "fig2_treatments.RDS"))

# Unnest the data so that all the results are in one dataframe
df_res <- select(res, -out) |>
  tidytable::unnest(mod_df)

# changes in K
p_k <- ggplot() +
  geom_line(data = filter(df_res, D == 4*3600, A_stor_frac == 0.2,
                          type == "DO", reach == 50),
            aes(x = time, y = value, color = as.factor(K), group = K),
            size = 1.5) +
  scale_color_viridis_d(expression(K~"("*d^{-1}*")")) +
  # geom_line(data = filter(df_res, D == 0, type == "DO", reach == 50),
  #           aes(x = time, y = value), color = "black",
  #           size = 1) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.12, 0.7),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = "transparent"))
p_k

# changes in D
p_d <- ggplot() +
  geom_line(data = filter(df_res, K == 2, A_stor_frac == 0.2, D!=0,
                          type == "DO", reach == 50) |>
              mutate(D = D / 3600),
            aes(x = time, y = value, color = as.factor(D), group = D),
            size = 1.5) +
  scale_color_viridis_d(expression(D~"("*m^{2}~s^{-1}*")")) +
  # geom_line(data = filter(df_res, D == 0, type == "DO", reach == 50),
  #           aes(x = time, y = value), color = "black",
  #           size = 1) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.12, 0.7),
        legend.key.size = unit(0.25, "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))
p_d

# Changes in storage
p_s <- ggplot() +
  geom_line(data = filter(df_res, K == 2, D == 4*3600, type == "DO", reach == 50),
            aes(x = time, y = value, color = as.factor(A_stor_frac), group = A_stor_frac),
            size = 1.5) +
  scale_color_viridis_d(expression(A[S]~"(-)")) +
  # geom_line(data = filter(df_res, D == 0, type == "DO", reach == 50),
  #           aes(x = time, y = value), color = "black",
  #           size = 1) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.12, 0.7),
        legend.key.size = unit(0.25, "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))
p_s

p <- p_k | p_d | p_s
p

ggsave(plot = p,
       filename = file.path("results", "replicate_fig2.png"),
       dpi = 300,
       units = "cm",
       height = 9,
       width = 18.4)
