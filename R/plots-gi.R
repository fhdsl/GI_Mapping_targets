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
  #NEED the mean_expected_CS, mean_observed_CS, broad_target_type (single_targeting, double_targeting) columns for this dataset
  regression_data <- old_gi_results %>%
    filter(broad_target_type == "single_targeting")

  model <- lm(mean_observed_CS ~ mean_expected_CS, data = regression_data)

  old_gi_results %>% 
    ggplot(aes(x=mean_expected_CS, 
               y=mean_observed_CS, 
               color=broad_target_type)) + 
    geom_point(alpha=0.7) + 
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

plot_rank_scatter <- function(old_gi_results){
  ## Only need the mean_GI_score column for this graph (they said they didn't need the inset)
  old_gi_results %>%
    mutate(Rank = dense_rank(mean_GI_score)) %>%
    ggplot(aes(x=Rank, 
               y=mean_GI_score)) + 
    geom_point(alpha=0.7) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank()) + 
    theme_classic() +
    ylab("GI score") +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -0.5, linetype = "dashed") +
    geom_hline(yintercept = 0.25, linetype = "dashed")
}

plot_volcano <- function(old_gi_results){
  ## NEED THE FDR and mean_GI_score columns for this graph.
  old_gi_results %>%
    mutate(logfdr = -log10(fdr),
           pointColor = case_when(logfdr < 1 ~ "darkgrey",
                                  ((mean_GI_score < -0.5) & (logfdr > 1)) ~ "dodgerblue3",
                                  ((mean_GI_score > 0.25) & (logfdr > 1)) ~ "darkred",
                                  .default="black")) %>%
    ggplot(aes(x = mean_GI_score,
               y = logfdr,
               color = pointColor)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.25, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_color_manual(values = c("darkgrey" = "darkgrey", 
                                  "dodgerblue3" = "dodgerblue3", 
                                  "darkred" = "darkred",
                                  "black" = "black")) +
    ylab("-log10(FDR)")
}


