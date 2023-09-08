#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To define functions for comparison of DO model results with streammetabolizer
# Date: November 19, 2021
# 

# Load libraries
library(deSolve)
library(tidyverse)
library(ReacTran)
library(lubridate)

# Function to run the 1D model with choices of treatments
func_mod <- function(...){
  
  # Get the changes to parameters
  arguments <- list(...)
  
  # Model with user specifications
  ode.1D(y = yini,
         times = times,
         func = model,
         parms = replace(parms,
                         names(arguments),
                         unlist(arguments)),
         nspec = 2,
         dimens = with(as.list(parms), L / dx))
}

# Turn model outputs into dataframes to compare with streammetabolizer
# function to do so
df_fun <- function(mod, reach_choice, pars = parms){
  as_tibble(mod) |>
    tidytable::pivot_longer(cols = -time, names_to = "key", values_to = "value") |>
    tidytable::separate_wider_regex(key, c(type = "[A-Za-z]+",  reach = "\\d+")) |>
    tidytable::mutate(reach = as.numeric(reach),
                      dist = as.numeric(reach) * with(as.list(pars), dx)) |>
    tidytable::mutate(type_plot = recode(type,
                                         `DO` = "DO~(mg~L^{-1})",
                                         `DO_stor` = "DO[storage]~(mg~L^{-1})",
                                         `ER.DO_stor` = "ER~(g~O^2~m^{-2}~h^{-1})",
                                         `ER` = "ER~(g~O^2~m^{-2}~h^{-1})",
                                         `GPP` = "GPP~(g~O^2~m^{-2}~h^{-1})",
                                         `PAR` = "PAR~({`mu`}*mol~m^{-2}~s^{-1})"))
}

