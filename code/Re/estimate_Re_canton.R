###########################################################
# Estimating Re for B117 and non-B117 per canton
# Author: JS Huisman

# In contrast to the other script this one only writes
# values to csv (rather than plotting)
###########################################################
library(tidyverse)
library(lubridate)
library(viridis)
library(EpiEstim)

# This code assumes you have cloned the git repository
# https://github.com/covid-19-Re/shiny-dailyRe in the location
# app_dir
app_dir = '<shiny-dailyRe>'
data_dir = '../data'
plot_dir = '../figures/cantons'

source(paste0(app_dir,'/app/otherScripts/2_utils_getInfectionIncidence.R'))
source(paste0(app_dir,'/app/otherScripts/3_utils_doReEstimates.R'))
source('Re_functions.R')

###########################################################
## Some functions for readability ####
deconvolve_b117 <- function(data, n_boot = 1000){
  orig_deconv <- deconvolveIncidence(data, incidence_var = 'original',
                                     getCountParams('incubation'), 
                                     getCountParams('confirmed'),
                                     smooth_param = TRUE, n_boot = n_boot)
  
  b117_deconv <- deconvolveIncidence(data, incidence_var = 'b117',
                                     getCountParams('incubation'), 
                                     getCountParams('confirmed'),
                                     smooth_param = TRUE, n_boot = n_boot)
  
  all_deconv <- bind_rows(orig_deconv, b117_deconv) %>%
    select(data_type, date, replicate, value) %>%
    mutate(data_type = ifelse(data_type == 'infection_b117', 'B117', 'Non-B117'))
  
  return(list(orig = orig_deconv, 
              b117 = b117_deconv,
              all = all_deconv))
}

estimate_Re_b117 <- function(deconv_result, interval_ends,
                             estimate_types = c("step")){
  orig_Re <- getReBootstrap(deconv_result$orig, 
                            interval_ends = interval_ends,
                            estimate_types = estimate_types)
  b117_Re <- getReBootstrap(deconv_result$b117, interval_ends = interval_ends,
                            estimate_types = estimate_types)
  
  all_Re <- bind_rows(orig_Re, b117_Re) %>%
    select(date, data_type, estimate_type, median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
    mutate(data_type = ifelse(data_type == 'infection_b117', 'B117', 'Non-B117'))
  
  return(all_Re)
}
###########################################################
## General settings ####
start_date = as_date("2021-01-01")
data_end_date = as_date("2021-02-11")
##########################
## Read raw Data ####

data_files <- list.files(path = data_dir, full.names = TRUE)

Jan_estimates <- data.frame()
for (data_file in data_files){
  raw_data = read_csv(data_file)
  data <- raw_data %>%
    filter(date <= data_end_date ) %>%
    #mutate(orig_data = TRUE) %>%
    complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
    mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) 
  
  ## Deconvolve
  deconv_results <- deconvolve_b117(data, n_boot = 100)
  
  # Estimate Re
  intervals = list('CHE' = c(start_date-1, "2021-01-17"))
  all_Re <- estimate_Re_b117(deconv_results, interval_ends = intervals,
                             estimate_types = c("step"))
  
  new_Jan_estimates <- all_Re %>%
    filter(date == start_date) %>%
    select(data_type, R = median_R_mean) %>%
    mutate(dataset = gsub("../data/estimated_case_numbers_|.csv", "", data_file) )
    
  Jan_estimates <- bind_rows(Jan_estimates, new_Jan_estimates)
}

write_csv(Jan_estimates, paste0(plot_dir, '/', 'Jan_estimates.csv'))

