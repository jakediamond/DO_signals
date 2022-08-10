#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with estimates from streammetabolizer
# Date: 9 August 2022
# 

# Load libraries
library(furrr)
library(scales)
library(patchwork)
library(tidyverse)

# Set parallel computing 
plan(multisession, workers = 8)
# Set seed for reproducibility
set.seed(42)

# Load the model
source(file.path("src", "1D_DO_model.R"))

# Load the functions
source(file.path("src", "model_streammetabolizer_comparison_functions.R"))

# Get a table of the different models to run
trts <- tibble(gpp_choice = c(rep("sine",6), "constant"),
               phase = c(0,0.25,-0.25,0,0,0,0),
               f = c(1,1,1,2,4,8,0),
               name = c("f==1", "{f==1}~{phi==+0.25}",
                        "{f==1}~{phi==-0.25}", "f==2",
                        "f==4", "f==8", "uniform"))

# Do the models
mods <- trts %>%
  mutate(out = future_pmap(list(f = f, phase = phase, gpp_choice = gpp_choice), 
                           func_mod2, .options = furrr_options(seed = TRUE)))

# Convert to usable dataframe
res <- mods %>%
  mutate(data = future_map(out, df_fun))

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

# Colors
cols <- c("dark red", "dark red", "dark red",
        "dark blue", "blue", "light blue", "black")

# Linetypes
lts <- c("solid", "dashed", "dotted", "solid", "solid",
         "solid", "solid")

# Function to parse labels
parse.labels <- function(x) parse(text = x, srcfile = NULL)

# First panel, GPP vs distance 
p_a <- select(res, -out) %>%
  mutate(data = map(data, ~ .x %>% 
                      filter(time_hr == 60,
                             type == "GPP") %>%
                      mutate_all(as.character))) %>% 
  tidyfast::dt_unnest(data) %>%
  type_convert() %>%
  mutate(part = if_else(f==1, "part2", "part1")) %>%
  ggplot(aes(x = dist,
             y = value,
             color = name,
             linetype = name)) +
  geom_line(size = 1.2)+
  scale_x_continuous(breaks = seq(0,20000, 5000),
                     expand = c(0,0)) +
  scale_color_manual(name = "",
                     values = cols,
                     labels = parse.labels) +
  scale_linetype_manual(name = "",
                        values = lts,
                        labels = parse.labels) +
  facet_grid(rows = vars(part)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  labs(x = "distance (m)",
       y = expression("GPP (g "*O[2]~m^{-2}~d^{-1}*")"))
p_a


p_b <- res_p %>%
  filter(between(time_hr, 48, 96),
         type == "DO") %>%
  ggplot(aes(x = time_hr-48,
           y = value,
           color = name,
           linetype = name)) +
  geom_line(size = 1.2)+
  scale_x_continuous(breaks = c(0,24,48),
                     expand = c(0,0)) +
  scale_color_manual(name = "",
                     values = cols,
                     labels = parse.labels) +
  scale_linetype_manual(name = "",
                        values = lts,
                        labels = parse.labels) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.text.align = 0,
        legend.position = c(0.12, 0.84),
        legend.background = element_rect(fill = "transparent")) +
  labs(x = "time (hour)",
       y = expression("dissolved oxygen (mg "*L^{-1}*")"))
p_b

p_c <- select(res, -out) %>%
  mutate(data = map(data, ~ .x %>% 
                      filter(time_hr %in% c(48, 52, 56, 60, 64, 68),
                             type == "DO") %>%
                      mutate_all(as.character))) %>% 
  tidyfast::dt_unnest(data) %>%
  type_convert() %>%
  filter(name != "uniform") %>%
  mutate(t = factor(paste0(time_hr%%24,":","00"),
                       levels = c("0:00", "4:00", "8:00", "12:00", "16:00", "20:00")),
         name = fct_relevel(name, c("f==8", "f==4","f==2", "f==1",
                                    "{f==1}~{phi==-0.25}","{f==1}~{phi==+0.25}"))) %>%
  ggplot(aes(x = dist,
             y = value,
             color = t,
             group = t)) +
  geom_line(size = 1.2)+
  scale_x_continuous(breaks = seq(0,20000, 5000),
                     expand = c(0,0)) +
  scale_color_viridis_d(name = "",
                        option = "magma") +
  facet_grid(rows = vars(name)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal") +
  labs(x = "distance (m)",
       y = expression("dissolved oxygen (mg "*L^{-1}*")"))
p_c


p_fig4 <- (p_a/p_b) | p_c
p_fig4
