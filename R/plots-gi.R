# Genetic interaction calculation step related plots can go here

library(tidyverse)

## This plot is meant to be functionally equivalent to Fig S5K (for HeLa, equivalent of Fig 3a for PC9).
#From the figure's caption:
#Scatter plot of target-level observed versus expected CRISPR scores in the HeLa screen.
#The solid line is the linear regression line for the negative control (single KO) pgRNAs, 
#while dashed lines indicate ± 2 residuals.

## Input is the old results from GI Mapping (old_gi_results <- readRDS("d.HeLa_GI_scores_target"))
## column names include broad_target_type, mean_observed_CS, and mean_expected_CS (as well as others, but those are the relevant ones for this plot)


plot_main_scatter <- function(old_gi_results){

  regression_data <- old_gi_results %>%
    filter(broad_target_type == "single_targeting")

  model <- lm(mean_observed_CS ~ mean_expected_CS, data = regression_data)

  old_gi_results %>% 
    ggplot(aes(x=mean_expected_CS, 
               y=mean_observed_CS, 
               color=broad_target_type)) + 
    geom_point(alpha=0.8) + 
    scale_color_manual(values = c("single_targeting" = "gray50", 
                                  "double_targeting" = "mediumpurple3"), 
                       labels = c("single_targeting" = "control", 
                                  "double_targeting" = "DKO")) + 
    theme(panel.background = element_blank(), 
          panel.grid = element_blank()) + 
    theme_classic() +
    scale_x_continuous(breaks=seq(-5, 1), 
                       labels = seq(-5, 1)) + 
    theme(legend.position = "inside", 
          legend.position.inside = c(0.2, 0.8)) + 
    theme(legend.title = element_blank()) + 
    xlab("Expected CRISPR score\n(paralog 1 KO + paralog 2 KO)") + 
    ylab("Observed CRISPR score\n(paralog 1 & 2 DKO)") +
    geom_abline(slope = model$coefficients[["mean_expected_CS"]], 
                intercept = model$coefficients[["(Intercept)"]]) +
    geom_abline(slope = model$coefficients[["mean_expected_CS"]], 
                intercept = model$coefficients[["(Intercept)"]] + quantile(model$residuals)["75%"], 
                linetype=3) +
    geom_abline(slope = model$coefficients[["mean_expected_CS"]], 
                intercept = model$coefficients[["(Intercept)"]] + quantile(model$residuals)["25%"], 
                linetype=3)
}

# Not sure if this is handling the +- residuals lines correctly or not

# Want to rewrite this to work with the new gimap results as input instead
