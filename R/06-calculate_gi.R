#' Calculate Genetic Interaction scores
#' @description Create results table that has CRISPR scores, Wilcoxon rank-sum test and t tests.
#' @param .data Data can be piped in with tidyverse pipes from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @export
#' @examples {
#'   gimap_dataset <- get_example_data("gimap")
#'
#'   # Highly recommended but not required
#'   run_qc(gimap_dataset)
#'
#'   gimap_dataset <- gimap_dataset %>%
#'     gimap_filter() %>%
#'     gimap_annotate(cell_line = "HELA") %>%
#'     gimap_normalize(
#'       timepoints = "day",
#'       replicates = "rep"
#'     ) %>%
#'     calc_crispr() %>%
#'     calc_gi()
#'
#'   saveRDS(gimap_dataset, "gimap_dataset_final.RDS")
#' }
calc_gi <- function(.data = NULL,
                    gimap_dataset) {
  # Code adapted from
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/04-calculate_GI_scores.Rmd

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  ## calculate expected double-targeting GI score by summing the two mean single-targeting
  ## CRISPR scores for that paralog pair
  gi_calc_df <- gimap_dataset$crispr_score %>%
    dplyr::mutate(expected_crispr = dplyr::case_when(
      target_type == "gene_gene" ~ mean_single_target_crispr_1 + mean_single_target_crispr_2,
      target_type == "gene_ctrl" ~ mean_single_target_crispr_1 + mean_double_control_crispr_2,
      target_type == "ctrl_gene" ~ mean_double_control_crispr_1 + mean_single_target_crispr_2
    ))
  
  # Calculating mean crisprs
  overall_results <- gi_calc_df %>%
    dplyr::group_by(rep, pgRNA_target_double) %>%
    dplyr::summarize(
      mean_expected_crispr = mean(expected_crispr, na.rm = TRUE),
      mean_observed_crispr = mean(double_crispr_score, na.rm = TRUE)
    ) %>%
    group_modify(~ broom::tidy(lm(mean_observed_crispr ~ mean_expected_crispr, data = .x)))

  # Run the overall linear model
  stats <- overall_results %>%
    dplyr::ungroup() %>%
    dplyr::select(term, estimate, rep) %>%
    pivot_wider(
      names_from = term,
      values_from = estimate
    ) %>%
    rename(intercept = "(Intercept)", slope = mean_expected_crispr)

  message("Calculating Genetic Interaction scores")
  # Do the linear model adjustments
  gi_calc_adj <- gi_calc_df %>%
    dplyr::left_join(stats, by = "rep") %>%
    dplyr::mutate(
      double_target_gi_score = double_crispr_score - (intercept + slope * expected_crispr),
      single_target_gi_score_1 = single_crispr_score_1 - (intercept + slope * expected_crispr),
      single_target_gi_score_2 = single_crispr_score_2 - (intercept + slope * expected_crispr)
    )

  # Which replicates we have?
  replicates <- unique(gi_calc_adj$rep)

  # TODO: this is an lapply inside an lapply at the end of the day. Not ideal
  target_results <- lapply(replicates,
    gimap_rep_stats,
    gi_calc_adj = gi_calc_adj
  )

  # Bring over the replicate names
  names(target_results) <- replicates

  # Turn into a data.frame
  target_results_df <- dplyr::bind_rows(target_results, .id = "replicate") %>%
    tidyr::pivot_wider(
      names_from = replicate,
      values_from = c(
        p_val_ttest,
        p_val_wil,
        fdr_vals_ttest,
        fdr_vals_wil
      )
    )

  # Store the useful bits
  gimap_dataset$gi_scores <- gi_calc_adj %>%
    dplyr::select(
      pgRNA_target_double,
      rep,
      double_target_gi_score,
      single_target_gi_score_1,
      single_target_gi_score_2,
      expected_crispr
    )

  # Store this
  gimap_dataset$results <-
    list(
      overall = overall_results,
      by_target = target_results_df
    )

  return(gimap_dataset)
}


#' Do tests for each replicate --an internal function used by calc_gi() function
#' @description Create results table that has t test p values
#' @param replicate a name of a replicate to filter out from gi_calc_adj
#' @param gi_calc_adj a data.frame with adjusted gi scores
gimap_rep_stats <- function(replicate, gi_calc_adj) {
  ## get a vector of GI scores for all single-targeting ("control") pgRNAs for each rep
  ## get double-targeting pgRNAs for this rep, do a t-test to compare the double-

  ## targeting GI scores for each paralog pair to the control vector

  ## adjust for multiple testing using the Benjamini-Hochberg method

  per_rep_stats <- gi_calc_adj %>%
    dplyr::filter(rep == replicate, pgRNA_target_double != "ctrl_ctrl")

  # Extract double scores
  double_scores <- per_rep_stats %>%
    dplyr::select(pgRNA_target_double, double_target_gi_score)

  # Extract single target scores
  single_scores <- per_rep_stats %>%
    dplyr::select(
      pgRNA_target_double, single_target_gi_score_1,
      single_target_gi_score_2
    )

  # What's the target list
  double_targets <- unique(double_scores$pgRNA_target_double)

  # Run the test for each target
  p_vals <- lapply(double_targets, function(target) {
    # Get the values for this particular target
    doubles <- dplyr::filter(double_scores, pgRNA_target_double == target) %>%
      dplyr::distinct()
    singles <- dplyr::filter(single_scores, pgRNA_target_double == target) %>%
      dplyr::distinct()

    p_val_ttest <- t.test(
      x = doubles$double_target_gi_score,
      y = c(singles$single_target_gi_score_1, singles$single_target_gi_score_2),
      paired = FALSE
    )$p.value

    p_val_wil <- wilcox.test(
      x = doubles$double_target_gi_score,
      y = c(singles$single_target_gi_score_1, singles$single_target_gi_score_2),
      paired = FALSE, exact = FALSE
    )$p.value

    p_vals <- data.frame(
      p_val_ttest,
      p_val_wil
    )

    return(p_vals)
  })

  # Put this together in a dataframe for this replicate
  p_vals_df <- data.frame(
    targets = double_targets,
    bind_rows(p_vals)
  )

  # Adjust for multiple testing using the Benjamini-Hochberg method
  p_vals_df$fdr_vals_ttest <- p.adjust(p_vals_df$p_val_ttest, method = "BH")
  p_vals_df$fdr_vals_wil <- p.adjust(p_vals_df$p_val_wil, method = "BH")

  return(p_vals_df)
}
