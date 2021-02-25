###########################################################
# Estimating Re for B117 and non-B117
# Author: JS Huisman
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
plot_dir = '../figures'

source(paste0(app_dir,'/app/otherScripts/2_utils_getInfectionIncidence.R'))
source(paste0(app_dir,'/app/otherScripts/3_utils_doReEstimates.R'))
###########################################################
source('Re_functions.R')
###########################################################
## Read raw Data ####

##########################
## Just rerun everything below with the different setting; 
# will save separate plots

#data_source = 'viollier'
#data_source = 'risch'
#start_date = as_date("2021-01-01")
#start_date = as_date("2020-12-15")

start_dates = c(as_date("2020-12-22"), as_date("2021-01-01"))

for (data_source in c('viollier', 'risch')){
  for (start_date_ind in 1:2 ){
    start_date = start_dates[start_date_ind]
    
    ##########################
    # there is roughly a 10 day delay between infection and case confirmation
    data_end_date = as_date("2021-02-11")
    plot_end_date = as_date("2021-02-01")
    
    raw_data = read_csv(paste0(data_dir, '/estimated_case_numbers_', data_source, '.csv'))
    data <- raw_data %>%
      filter(date <= data_end_date )
    
    if (start_date == as_date("2021-01-01")){
      plot_dir = '../figures/Jan1'
    } else {
      plot_dir = '../figures/Dec22'
    }
    
    ##########################
    plot_data = data %>% 
      pivot_longer(cols = c('b117', 'original'))
    
    ggplot(plot_data, aes(x=date, y = value, fill = name)) +
      geom_bar(position = 'stack', stat = 'identity', alpha = 1) +
      #geom_bar(position = 'identity', stat = 'identity', alpha = 0.5) +
      facet_wrap(vars(name), ncol = 1, scale = 'free_y') +
      labs(x = 'Date' , y='Cases per day', fill = 'Variant') +
      scale_x_date(date_breaks = "1 weeks", 
                   date_labels = '%b\n%d',
                   limits = c(start_date, plot_end_date)) +
      #scale_y_continuous(trans = 'log') +
      scale_fill_manual(values = c("#0D4A70", "#67916E")) + 
      theme_minimal() +
      theme(
        strip.text.y= element_text(size=20),
        strip.text.x= element_text(size=20),
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)
      )
    
    ggsave(paste0(plot_dir, '/', 'B117_raw_', data_source, '.pdf'), height = 11, width = 16)
    # ggsave(paste0(plot_dir, '/', 'B117_raw_', data_source, '.png'), height = 11, width = 16)
    
    
    ## Deconvolve Data ####
    orig_deconv <- deconvolveIncidence(data, incidence_var = 'original',
                                       getCountParams('incubation'),
                                       getCountParams('confirmed'),
                                       smooth_param = TRUE, n_boot = 1000)
    
    b117_deconv <- deconvolveIncidence(data, incidence_var = 'b117',
                                       getCountParams('incubation'),
                                       getCountParams('confirmed'),
                                       smooth_param = TRUE, n_boot = 1000)
    
    all_deconv <- bind_rows(orig_deconv, b117_deconv) %>%
      select(data_type, date, replicate, value) %>%
      mutate(data_type = ifelse(data_type == 'infection_b117', 'B117', 'Non-B117'))
    
    ## Plot Deconvolved cases ####
    mean_deconv <- all_deconv %>%
      group_by(date, data_type) %>%
      summarise(sd = sd(value),
                value = mean(value),
                .groups = 'drop')
    
    write_csv(mean_deconv, paste0(plot_dir, '/', 'mean_deconv_', data_source, '.csv'))
    write_csv(all_deconv, paste0(plot_dir, '/', 'all_deconv_', data_source, '.csv'))
    
    ggplot() +
      geom_errorbar(data = mean_deconv,
                    aes(x=date, ymin = value -sd,  ymax = value +sd, colour = data_type)) +
      labs(x = 'Date' , y='Infections per day') +
      scale_x_date(date_breaks = "1 weeks",
                   date_labels = '%b\n%d',
                   limits = c(start_date, plot_end_date)) +
      #scale_y_continuous(trans = 'log') +
      scale_colour_manual(values = c("#0D4A70", "#67916E")) +
      theme_minimal() +
      theme(
        strip.text.y= element_text(size=20),
        strip.text.x= element_text(size=20),
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)
      )
    
    ggsave(paste0(plot_dir, '/', 'B117_deconv_', data_source, '.pdf'), height = 11, width = 16)
    # ggsave(paste0(plot_dir, '/', 'B117_deconv_', data_source, '.png'), height = 11, width = 16)
    
    ## Estimate Re ####
    intervals = list('CHE' = c(start_date-1, "2021-01-17"))
    
    orig_Re <- getReBootstrap(orig_deconv, interval_ends = intervals)
    b117_Re <- getReBootstrap(b117_deconv, interval_ends = intervals)
    
    all_Re <- bind_rows(orig_Re, b117_Re) %>%
      select(date, data_type, estimate_type, median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
      mutate(data_type = ifelse(data_type == 'infection_b117', 'B117', 'Non-B117'))
    
    write_csv(all_Re, paste0(plot_dir, '/', 'all_Re_', data_source, '.csv'))
    
    ## Plot Re ####
    method_names <- setNames(c("Sliding Window","Piecewise constant"),
                             c('Cori_slidingWindow', 'Cori_step'))
    
    ggplot(all_Re) +
      geom_ribbon(aes(x = date, ymin = median_R_lowHPD,
                      ymax = median_R_highHPD, fill = data_type), alpha = 0.5, show.legend = F) +
      geom_line(aes(x = date, y = median_R_mean, colour = data_type), alpha = 1, show.legend = T) +
      facet_wrap(vars(estimate_type), ncol = 1,
                 scale = 'free_y', labeller = labeller(estimate_type = method_names)) +
      geom_hline(yintercept = 1) +
      #coord_cartesian(ylim = c(0, 4)) +
      scale_colour_manual(values = c("#0D4A70", "#67916E")) +
      scale_fill_manual(values = c("#0D4A70", "#67916E")) +
      scale_x_date(date_breaks = "1 weeks",
                   date_labels = '%b\n%d',
                   limits = c(start_date, plot_end_date)) +
      labs( x = 'Date', y = 'Estimated reproductive number R',
            colour = 'Observation Type', fill = 'Observation Type') +
      theme_minimal() +
      guides(color = guide_legend(override.aes = list(size=5))) +
      theme(
        strip.text= element_text(size = 25),
        #strip.text.x= element_text(size=20),
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'bottom'
      )
    
    ggsave(paste0(plot_dir, '/', 'B117_Re_', data_source, '.pdf'), height = 13, width = 14)
    # ggsave(paste0(plot_dir, '/', 'B117_Re_', data_source, '.png'), height = 12, width = 14)
    
    ## Plot Re Difference / Ratio ####
    
    compare_Re <- orig_Re %>%
      select(date, data_type, estimate_type, median_R_mean, median_R_highHPD, median_R_lowHPD, sd_mean) %>%
      right_join(b117_Re, by = c('date', 'estimate_type')) %>%
      group_by(estimate_type) %>%
      mutate(diff_mean = median_R_mean.y - median_R_mean.x,
             diff_sd = sqrt(sd_mean.y^2 + sd_mean.x^2),
             ratio_mean = median_R_mean.y/median_R_mean.x,
             ratio_sd = sqrt((1/median_R_mean.x^2)*sd_mean.y^2 +
                               (median_R_mean.y^2/median_R_mean.x^4)*sd_mean.x^2)
      ) %>%
      ungroup()
    
    ggplot(compare_Re) +
      geom_ribbon(aes(x = date, ymin = diff_mean - diff_sd,
                      ymax = diff_mean + diff_sd), alpha = 0.5, show.legend = F) +
      geom_line(aes(x = date, y = diff_mean), alpha = 1, show.legend = T) +
      facet_wrap(vars(estimate_type), ncol = 1,
                 scale = 'free_y', labeller = labeller(estimate_type = method_names)) +
      coord_cartesian(ylim = c(-0.5, 1.5)) +
      scale_x_date(date_breaks = "1 weeks",
                   date_labels = '%b\n%d',
                   limits = c(start_date, plot_end_date)) +
      labs( x = 'Date', y = 'R(B117) - R(non-B117)',
            colour = 'Observation Type', fill = 'Observation Type') +
      theme_minimal() +
      guides(color = guide_legend(override.aes = list(size=5))) +
      theme(
        #strip.text.y= element_blank(),
        #strip.text.x= element_text(size=20),
        strip.text= element_text(size = 25),
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'bottom'
      )
    
    ggsave(paste0(plot_dir, '/', 'Re_diff_', data_source, '.pdf'), height = 12, width = 14)
    # ggsave(paste0(plot_dir, '/', 'Re_diff_', data_source, '.png'), height = 12, width = 14)
    
    
    ggplot(compare_Re) +
      geom_ribbon(aes(x = date, ymin = ratio_mean - ratio_sd,
                      ymax = ratio_mean + ratio_sd), alpha = 0.5, show.legend = F) +
      geom_line(aes(x = date, y = ratio_mean), alpha = 1, show.legend = T) +
      facet_wrap(vars(estimate_type), ncol = 1,
                 scale = 'free_y', labeller = labeller(estimate_type = method_names)) +
      coord_cartesian(ylim = c(1, 2.5)) +
      scale_x_date(date_breaks = "1 weeks",
                   date_labels = '%b\n%d',
                   limits = c(start_date, plot_end_date)) +
      labs( x = 'Date', y = 'Ratio of R(B117) to R(non-B117)',
            colour = 'Observation Type', fill = 'Observation Type') +
      theme_minimal() +
      guides(color = guide_legend(override.aes = list(size=5))) +
      theme(
        #strip.text.y= element_blank(),
        #strip.text.x= element_text(size=20),
        strip.text= element_text(size = 25),
        axis.text.y= element_text(size=20),
        axis.text.x= element_text(size=20),
        axis.title.y =  element_text(size=25),
        axis.title.x =  element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = 'bottom'
      )
    
    ggsave(paste0(plot_dir, '/', 'Re_ratio_', data_source, '.pdf'), height = 12, width = 14)
    # ggsave(paste0(plot_dir, '/', 'Re_ratio_', data_source, '.png'), height = 12, width = 14)
    
  }
}
