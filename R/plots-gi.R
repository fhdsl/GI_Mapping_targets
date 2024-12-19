# Genetic interaction calculation step related plots can go here

library(tidyverse)

## This plot is meant to be functionally equivalent to Fig S5K (for HeLa, equivalent of Fig 3a for PC9).
#From the figure's caption:
#Scatter plot of target-level observed versus expected CRISPR scores in the HeLa screen.
#The solid line is the linear regression line for the negative control (single KO) pgRNAs, 
#while dashed lines indicate ± 2 residuals.

plot_exp_v_obs_scatter <- function(gimap_dataset, facet_rep = TRUE, reps_to_drop = c("Day05_RepA_early")){
  regression_data <- gimap_dataset$gi_scores %>%
    filter(target_type != "gene_gene") %>% #get only single targeting
    filter(!(rep %in% reps_to_drop))
  
  model <- lm(mean_observed_cs ~ mean_expected_cs, data = regression_data)
  
  gplot <- gimap_dataset$gi_scores %>%
    filter(!(rep %in% reps_to_drop)) %>%
    mutate(broad_target_type = case_when(target_type == "gene_gene" ~ "DKO",
                                         target_type == "ctrl_gene" ~ "control",
                                         target_type == "gene_ctrl" ~ "control")) %>%
    ggplot(aes(x=mean_expected_cs,
               y=mean_observed_cs,
               color = broad_target_type)) +
    geom_point(size=1, alpha=0.7) +
    scale_color_manual(values = c("control" = "gray50",
                                  "DKO" = "mediumpurple3")) + 
    theme_classic() +
    theme(legend.position = "inside", 
          legend.position.inside = c(0.1, 0.9),
          legend.title = element_blank()) +
    xlab("Expected CRISPR score\n(paralog 1 KO + paralog 2 KO)") + 
    ylab("Observed CRISPR score\n(paralog 1 & 2 DKO)") +
    geom_abline(slope = model$coefficients[["mean_expected_cs"]], 
                intercept = model$coefficients[["(Intercept)"]]) +
    geom_abline(slope = model$coefficients[["mean_expected_cs"]], 
                intercept = model$coefficients[["(Intercept)"]] + quantile(model$residuals)["75%"], 
                linetype=3) +
    geom_abline(slope = model$coefficients[["mean_expected_cs"]], 
                intercept = model$coefficients[["(Intercept)"]] + quantile(model$residuals)["25%"], 
                linetype=3)
    
  if (facet_rep){
    return(gplot + facet_wrap(~rep))
  } else{ return(gplot) }  
}

plot_rank_scatter <- function(gimap_dataset, reps_to_drop = c("Day05_RepA_early")){
  return(
    gimap_dataset$gi_scores %>%
      filter(target_type == "gene_gene") %>% #get only double targeting
      filter(!(rep %in% reps_to_drop)) %>%
      mutate(Rank = dense_rank(mean_gi_score)) %>%
      ggplot(aes(x=Rank, 
                 y=mean_gi_score,
                 color = rep)) + 
      geom_point(size=1, alpha=0.7) +
      theme_classic() +
      theme(legend.position = "inside", 
            legend.position.inside = c(0.9, 0.1),
            legend.title = element_blank()) +
      ylab("GI score") +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = -0.5, linetype = "dashed") +
      geom_hline(yintercept = 0.25, linetype = "dashed")
  )
}

plot_volcano <- function(gimap_dataset, facet_rep = TRUE, reps_to_drop = c("Day05_RepA_early")){
  gplot <- gimap_dataset$gi_scores %>%
    filter(target_type == "gene_gene") %>% #get only double targeting
    filter(!(rep %in% reps_to_drop)) %>%
    mutate(logfdr = -log10(fdr),
           pointColor = case_when(logfdr < 1 ~ "darkgrey",
                                  ((mean_gi_score < -0.5) & (logfdr > 1)) ~ "dodgerblue3",
                                  ((mean_gi_score > 0.25) & (logfdr > 1)) ~ "darkred",
                                  .default="black")) %>%
    ggplot(aes(x = mean_gi_score,
               y = logfdr,
               color = pointColor)) +
    geom_point(size=1, alpha=0.7) +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.25, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_color_manual(values = c("darkgrey" = "darkgrey", 
                                  "dodgerblue3" = "dodgerblue3", 
                                  "darkred" = "darkred",
                                  "black" = "black")) +
    ylab("-log10(FDR)")
  
  if (facet_rep){
    return(gplot + facet_wrap(~rep))
  } else{ return(gplot) }
}

plot_targets_bar <- function(gimap_dataset, target1, target2, reps_to_drop = c("Day05_RepA_early")){
  gimap_dataset$gi_scores %>%
    filter(!(rep %in% reps_to_drop)) %>%
    filter((grepl(target1,pgRNA_target)) | (grepl(target2, pgRNA_target))) %>%
    ggplot(aes(y = mean_observed_cs, #This is not the right column
               x = target_type,
               fill = target_type)) +
    geom_bar(stat = "identity") +
    geom_point(pch = 21, size=3) +
    theme_bw() +
    ylab("CRISPR score") +
    xlab("") +
    ggtitle(paste0(target1, "/", target2)) +
    geom_hline(yintercept = 0) +
    scale_x_discrete(labels = c("ctrl_gene" = paste0(target2, " KO"), #this assumes that target2 is the ctrl_{target2} gene
                                "gene_ctrl" = paste0(target1, " KO"), #this assumes that target1 is the {target1}_ctrl gene
                                "gene_gene" = "DKO")) +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1))
}

