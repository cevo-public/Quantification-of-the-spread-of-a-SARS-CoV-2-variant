###########################################################
#  (Helper) Functions to estimate Re for B117 and non-B117
# Author: JS Huisman
###########################################################

getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

# For citations on these parameters see Huisman et al
getCountParams <- function(obs_type){
  switch(obs_type,
         incubation = getGammaParams(5.3, 3.2),
         zero = list(shape = 0, scale = 0),
         death = getGammaParams(15.0, 6.9),
         hospitalisation = getGammaParams(5.1, 4.2),
         confirmed = getGammaParams(5.5, 3.8)
  )
}

# This is for compatibility with the code from shiny-dailyRe
addUselessColumns <- function(df, inc_var = 'original'){
  
  observation_df <- df %>%
    dplyr::select(date, value = all_of(inc_var)) %>%
    mutate(data_type = inc_var,
           source = 'ETH',
           region = 'CHE',
           variable = 'incidence',
           country = 'Switzerland',
           date_type = 'report',
           local_infection = TRUE)
  
  return(observation_df)
}

# Wrapper around code from shiny-dailyRe
deconvolveIncidence <- function(df, incidence_var = 'original',
                                IncubationParams, OnsetToCountParams,
                                smooth_param = TRUE, n_boot = 100){
  infection_df <- addUselessColumns(df, inc_var = incidence_var)
  
  constant_delay_distributions <- list("Cases" = get_vector_constant_waiting_time_distr(
    IncubationParams$shape, IncubationParams$scale,
    OnsetToCountParams$shape, OnsetToCountParams$scale),
    "Symptoms" = get_vector_constant_waiting_time_distr(
      IncubationParams$shape, IncubationParams$scale,
      0, 0))
  
  estimatedInfections <- get_infection_incidence_by_deconvolution(
    infection_df,
    is_local_cases = T,
    constant_delay_distribution = constant_delay_distributions[['Cases']],
    constant_delay_distribution_incubation = constant_delay_distributions[["Symptoms"]],
    max_iterations = 100,
    smooth_incidence = smooth_param,
    empirical_delays = tibble(),
    n_bootstrap = n_boot,
    verbose = FALSE)
  
  return(estimatedInfections)
}

# Wrapper around code from shiny-dailyRe
getReBootstrap <- function(deconvoluted_data, 
                           interval_ends = list('CHE' = c("2020-12-21", "2021-01-17", "2021-02-17")),
                           estimate_types = c("slidingWindow", "step")){
  
  all_delays <- lapply(unique(deconvoluted_data$data_type), function(x){ c(Cori = 0)})
  names(all_delays) <- unique(deconvoluted_data$data_type)
  
  truncations <- list(left = c(Cori = 5),
                      right = c(Cori = 0))
  
  rawReEstimates <- doAllReEstimations(
    deconvoluted_data,
    slidingWindow = 3,
    methods = c("Cori"),
    variationTypes = estimate_types,
    all_delays,
    truncations,
    interval_ends = interval_ends
  )
  
  cleanEstimates <- cleanCountryReEstimate(rawReEstimates, method = 'bootstrap',
                                           rename_types = F, report_sd = T, alpha=0.95)
  
  return(cleanEstimates)
}
