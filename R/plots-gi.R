#' Expected vs Observed CRISPR Scatterplot
#' @description This plot is meant to be functionally equivalent to Fig S5K (for HeLa, equivalent of Fig 3a for PC9).
#' Scatter plot of target-level observed versus expected CRISPR scores in the screen.
#' The solid line is the linear regression line for the negative control (single KO) pgRNAs,
#' while dashed lines indicate the lower and upper quartile residuals.
#' @param gimap_dataset A special dataset structure that is originally setup using `setup_data()` and has had gi scores calculated with `calc_gi()`.
#' @param facet_rep Should the replicates be wrapped with facet_wrap()?
#' @param reps_to_drop Names of replicates that should be not plotted (Optional)
#' @import dplyr
#' @import ggplot2
#' @export
#' @return A ggplot2 scatterplot of the target level observed vs expected
#' CRISPR scores.
#' @examples \donttest{
#'
#' gimap_dataset <- get_example_data("gimap") %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA") %>%
#'   gimap_normalize(
#'     timepoints = "day",
#'      missing_ids_file =  tempfile()
#'   ) %>%
#'   calc_gi()
#'
#' # To plot results
#' plot_exp_v_obs_scatter(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' plot_rank_scatter(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' plot_volcano(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#'}
plot_exp_v_obs_scatter <- function(gimap_dataset, facet_rep = TRUE, reps_to_drop = "") {
  if (!("gimap_dataset" %in% class(gimap_dataset))) {
    stop(
      "This function only works",
      "with gimap_dataset objects which can be made with the setup_data() function."
    )
  }

  if (is.null(gimap_dataset$overall_results)) {
    stop(
      "This function only works with",
      "gimap_dataset objects which have had gi calculated with calc_gi()"
    )
  }


  regression_data <- gimap_dataset$gi_scores %>%
    filter(target_type != "gene_gene") %>% # get only single targeting
    filter(!(rep %in% reps_to_drop))

  model <- lm(mean_observed_cs ~ mean_expected_cs, data = regression_data)

  gplot <- gimap_dataset$gi_scores %>%
    filter(!(rep %in% reps_to_drop)) %>%
    mutate(broad_target_type = case_when(
      target_type == "gene_gene" ~ "DKO",
      target_type == "ctrl_gene" ~ "control",
      target_type == "gene_ctrl" ~ "control"
    )) %>%
    ggplot(aes(
      x = mean_expected_cs,
      y = mean_observed_cs,
      color = broad_target_type
    )) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual(values = c(
      "control" = "gray50",
      "DKO" = "mediumpurple3"
    )) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    xlab("Expected CRISPR score\n(paralog 1 KO + paralog 2 KO)") +
    ylab("Observed CRISPR score\n(paralog 1 & 2 DKO)") +
    geom_abline(
      slope = model$coefficients[["mean_expected_cs"]],
      intercept = model$coefficients[["(Intercept)"]]
    ) +
    geom_abline(
      slope = model$coefficients[["mean_expected_cs"]],
      intercept = model$coefficients[["(Intercept)"]] + quantile(model$residuals)["75%"],
      linetype = 3
    ) +
    geom_abline(
      slope = model$coefficients[["mean_expected_cs"]],
      intercept = model$coefficients[["(Intercept)"]] + quantile(model$residuals)["25%"],
      linetype = 3
    )

  if (facet_rep) {
    return(gplot + facet_wrap(~rep))
  } else {
    return(gplot)
  }
}
#' Rank plot for target-level GI scores
#' @description This plot is meant to be functionally equivalent to Fig 5a (for HeLa, equivalent of Fig 3c for PC9).
#' Rank plot of target-level GI scores.
#' Dashed horizontal lines are for GI scores of 0.25 and -0.5
#' @param gimap_dataset A special dataset structure that is originally setup using `setup_data()` and has had gi scores calculated with `calc_gi()`.
#' @param reps_to_drop Names of replicates that should be not plotted (Optional)
#' @import dplyr
#' @import ggplot2
#' @return A ggplot2 rankplot of the target level genetic interaction scores.
#' @export
#' @examples \donttest{
#'
#' gimap_dataset <- get_example_data("gimap") %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA") %>%
#'   gimap_normalize(
#'     timepoints = "day"
#'   ) %>%
#'   calc_gi()
#'
#' # To plot results
#' plot_exp_v_obs_scatter(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' plot_rank_scatter(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' plot_volcano(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' }
plot_rank_scatter <- function(gimap_dataset, reps_to_drop = "") {
  if (!("gimap_dataset" %in% class(gimap_dataset))) {
    stop(
      "This function only works",
      "with gimap_dataset objects which can be made with the setup_data() function."
    )
  }

  if (is.null(gimap_dataset$overall_results)) {
    stop(
      "This function only works with",
      "gimap_dataset objects which have had gi calculated with calc_gi()"
    )
  }

  return(
    gimap_dataset$gi_scores %>%
      filter(target_type == "gene_gene") %>% # get only double targeting
      filter(!(rep %in% reps_to_drop)) %>%
      mutate(Rank = dense_rank(mean_gi_score)) %>%
      ggplot(aes(
        x = Rank,
        y = mean_gi_score,
        color = rep
      )) +
      geom_point(size = 1, alpha = 0.7) +
      theme_classic() +
      theme(legend.title = element_blank()) +
      ylab("GI score") +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = -0.5, linetype = "dashed") +
      geom_hline(yintercept = 0.25, linetype = "dashed")
  )
}

#' Volcano plot for GI scores
#' @description This plot is meant to be functionally equivalent to Fig 5b (for HeLa, equivalent of Fig 3d for PC9).
#' Volcano plot of target-level GI scores
#' Blue points are synthetic lethal paralog GIs with GI < 0.5 and FDR < 0.1; red points are buffering paralog GIs with GI > 0.25 and FDR < 0.1.
#' @param gimap_dataset A special dataset structure that is originally setup using `setup_data()` and has had gi scores calculated with `calc_gi()`.
#' @param facet_rep Should the replicates be wrapped with facet_wrap()?
#' @param reps_to_drop Names of replicates that should be not plotted (Optional)
#' @import dplyr
#' @import ggplot2
#' @return A ggplot2 volcano plot of the target level genetic interaction scores.
#' @export
#' @examples \donttest{
#'
#' gimap_dataset <- get_example_data("gimap") %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA") %>%
#'   gimap_normalize(
#'     timepoints = "day"
#'   ) %>%
#'   calc_gi()
#'
#' # To plot results
#' plot_exp_v_obs_scatter(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' plot_rank_scatter(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' plot_volcano(gimap_dataset, reps_to_drop = "Day05_RepA_early")
#' }
plot_volcano <- function(gimap_dataset, facet_rep = TRUE, reps_to_drop = "") {
  if (!("gimap_dataset" %in% class(gimap_dataset))) {
    stop(
      "This function only works",
      "with gimap_dataset objects which can be made with the setup_data() function."
    )
  }

  if (is.null(gimap_dataset$overall_results)) {
    stop(
      "This function only works with",
      "gimap_dataset objects which have had gi calculated with calc_gi()"
    )
  }

  gplot <- gimap_dataset$gi_scores %>%
    filter(target_type == "gene_gene") %>% # get only double targeting
    filter(!(rep %in% reps_to_drop)) %>%
    mutate(
      logfdr = -log10(fdr),
      pointColor = case_when(logfdr < 1 ~ "darkgrey",
        ((mean_gi_score < -0.5) & (logfdr > 1)) ~ "dodgerblue3",
        ((mean_gi_score > 0.25) & (logfdr > 1)) ~ "darkred",
        .default = "black"
      )
    ) %>%
    ggplot(aes(
      x = mean_gi_score,
      y = logfdr,
      color = pointColor
    )) +
    geom_point(size = 1, alpha = 0.7) +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.25, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_color_manual(values = c(
      "darkgrey" = "darkgrey",
      "dodgerblue3" = "dodgerblue3",
      "darkred" = "darkred",
      "black" = "black"
    )) +
    ylab("-log10(FDR)") +
    xlab("Mean GI score")

  if (facet_rep) {
    return(gplot + facet_wrap(~rep))
  } else {
    return(gplot)
  }
}
#' Target bar plot for CRISPR scores
#' @description This plot is for when you'd like to examine a target pair
#' specifically -- meant to be functionally equivalent to Fig 3b
#' CRISPR scores for representative synthetic lethal paralog pairs.
#' Data shown are the mean CRISPR score for each single KO or DKO target across
#'  three biological replicates with replicate data shown in overlaid points.
#' @param gimap_dataset A special dataset structure that is originally setup using
#'  `setup_data()` and has had gi scores calculated with `calc_gi()`.
#' @param target1 Name of the first target to be plotted e.g.
#' @param target2 Name of the second target to be plotted e.g.
#' @param reps_to_drop Names of replicates that should be not plotted (Optional)
#' @import dplyr
#' @import ggplot2
#' @export
#' @return A ggplot2 bar plot of the specific target's genetic interaction scores.
#' @examples \donttest{
#'
#' gimap_dataset <- get_example_data("gimap") %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA") %>%
#'   gimap_normalize(
#'     timepoints = "day"
#'   ) %>%
#'   calc_gi()
#'
#' # To plot results, pick out two targets from the gi_score table
#' head(dplyr::arrange(gimap_dataset$gi_score, fdr))
#'
#' # "NDEL1_NDE1" is top result so let's plot that
#' plot_targets_bar(gimap_dataset, target1 = "NDEL1", target2 = "NDE1")
#' }
plot_targets_bar <- function(gimap_dataset, target1, target2, reps_to_drop = "") {
  if (!("gimap_dataset" %in% class(gimap_dataset))) {
    stop(
      "This function only works",
      "with gimap_dataset objects which can be made with the setup_data() function."
    )
  }

  if (is.null(gimap_dataset$overall_results)) {
    stop(
      "This function only works with",
      "gimap_dataset objects which have had gi calculated with calc_gi()"
    )
  }

  return(
    gimap_dataset$gi_scores %>%
      filter(!(rep %in% reps_to_drop)) %>%
      filter((grepl(target1, pgRNA_target)) | (grepl(target2, pgRNA_target))) %>%
      ggplot(aes(
        y = mean_observed_cs, # This is not the right column
        x = target_type,
        fill = target_type
      )) +
      geom_bar(stat = "summary", fun.y = "mean") +
      geom_point(pch = 21, size = 3) +
      theme_bw() +
      ylab("CRISPR score") +
      xlab("") +
      ggtitle(paste0(target1, "/", target2)) +
      geom_hline(yintercept = 0) +
      scale_x_discrete(labels = c(
        "ctrl_gene" = paste0(target2, " KO"),
        # this assumes that target2 is the ctrl_{target2} gene
        "gene_ctrl" = paste0(target1, " KO"),
        # this assumes that target1 is the {target1}_ctrl gene
        "gene_gene" = "DKO"
      )) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  )
}
