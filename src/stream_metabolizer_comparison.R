#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To define functions for comparison of DO model results with streammetabolizer
# Date: November 19, 2021
# 

# Load libraries
library(tidyverse)

# Load functions
source(file.path("src", "model_initialization.R"))
source(file.path("src", "streammetabolizer_helper_functions.R"))
source(file.path("src", "functions_for_model.R"))

# Get data ready to model ------------------------------------------------------
# Load results
df_ramp <- readRDS(file.path("results", "ramped.RDS"))
df_sine <- readRDS(file.path("results", "sine.RDS"))

# Now add a column for the streammetabolizer input data
df_sm <- df_sine |> #df_ramp |>
  mutate(
    # This takes the modeled O2 and creates an input data frame 
    # for the streamMetabolizer program. The second argument is a reach number
    sm_data = map2(mod_df, 95, sm_input_fun))

# Now run stream metabolizer on that data ---------------------------------
df_sm <- df_sm |>
  mutate(sm = map(sm_data, sm_fun, 
                  modtype = "mle", # can prescribe mle or bayesian ("bayes")
                  gpp_fun = "linlight",  #can have gpp be a saturating function of light or linear ("linlight)
                  er_fun = "constant", #er is constant, can be a function of temp ("q10temp")
                  k = "none" #prescribe initial K600 value if you want (1/d), leave as 'none" you want the model to choose
                  )
         ) 

# Dataframe of comparisons between streammetabolizer and DO model
df_com <- df_sm |>
  mutate(temp = map_dbl(sm_data, ~mean(.$temp.water)), #stream temperature from model
         Sc = 1801 - 120.1 * temp + 3.782 * temp^2 -0.0476 * temp^3, #schmidt number based on temp
         k600_mod = parms["K"] / sqrt(600 / Sc), #get k600 from DO model instead of KO2
         modinfo = map2(mod_df, 95, do_mod_res)) |>
  select(gpp_choice, f, phase, sm, temp, k600_mod, modinfo) |> #ramp_rate,
  unnest(cols = c(sm, modinfo))

# Save these results so you don't have to run again
saveRDS(df_com, file.path("results", "sine_sm_mle.RDS"))

# Get into good format for plotting
df_com_l <- df_com |>
  filter(!str_detect(msgs.fit, "W|E")) |> #no days with errors in streammetabolizer
  select(-GPP.lower, -GPP.upper, -ER.lower, -ER.upper, -date_mod, -msgs.fit,
         -warnings, -errors, -temp, -date) |>
  rename(GPP_sm = GPP, ER_sm = ER, k600_sm = K600.daily) |>
  pivot_longer(cols = GPP_sm:GPP_mod,
               names_to = c("type", "model"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = model, values_from = value, values_fn = mean) %>%
  mutate(sm = if_else(type == "ER", sm * -1, sm)) |>
  mutate(pctdiff_loc = (sm - mod) / mod,
         pctdiff_glob = (sm - 12) / 12)

# Plot the data
p_loc <- ggplot(data = filter(df_com_l, f != 0, phase == 0), #ramp_rate != -1),
       aes(x = f, #ramp_rate,
           y = pctdiff_loc * 100)) +
  facet_wrap(~type) +
  geom_hline(yintercept = 0) +
  geom_point() +
  theme_classic() +
  labs(x = "frequency (cycles per integration length)",
       y = "diff. from model (%)",
       title = "difference from model in local reach")
p_loc

p_glob <- ggplot(data = filter(df_com_l, f != 0, phase == 0, #ramp_rate != -1,
                               type == "GPP"),
                 aes(x = f,
                     y = pctdiff_glob * 100)) +
  # facet_wrap(~type) +
  geom_hline(yintercept = 0) +
  geom_point() +
  theme_classic() +
  labs(x = "frequency (cycles per integration length)",
       y = "GPP diff. from model (%)",
       title = "difference from model average")
p_glob

p <- p_glob / p_loc
p
ggsave(plot = p,
       filename = file.path('results', "sm_sine_comparison.png"),
       dpi = 300,
       units = "cm",
       height = 12,
       width = 18)
