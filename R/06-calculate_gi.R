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
    dplyr::mutate(
      expected_crispr_double = mean_single_target_crispr_1 + mean_single_target_crispr_2,
      expected_crispr_single_1 = mean_single_target_crispr_1 + mean_double_control_crispr,
      expected_crispr_single_2 = mean_single_target_crispr_2 + mean_double_control_crispr
    )

  mean_double_expected_df <- gi_calc_df %>%
    dplyr::group_by(rep, pgRNA_target_double) %>%
    dplyr::summarize(
      mean_expected_double_crispr = mean(expected_crispr_double, na.rm = TRUE)
    )

  # Calculating mean crisprs
  reshaped_single_df <- gi_calc_df %>%
    # reshaping the data a bit so we can do the math later
    dplyr::select(
      rep,
      pgRNA_target_double,
      mean_single_target_crispr_1,
      mean_single_target_crispr_2,
      expected_crispr_single_1,
      expected_crispr_single_2,
      gene_symbol_1,
      gene_symbol_2
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        gene_symbol_1,
        gene_symbol_2
      ),
      values_to = "gene_symbol",
      names_to = "which_gene"
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        expected_crispr_single_1,
        expected_crispr_single_2
      ),
      values_to = "expected_crispr_single",
      names_to = "which_expected"
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        mean_single_target_crispr_1,
        mean_single_target_crispr_2
      ),
      values_to = "observed_crispr_single",
      names_to = "which_obs"
    )


  mean_single_expected_df <- reshaped_single_df %>%
    dplyr::group_by(rep, gene_symbol) %>%
    dplyr::summarize(
      mean_expected_single_crispr = mean(expected_crispr_single, na.rm = TRUE)
    )
  # Calculate the linear model from this
  single_lm <- reshaped_single_df %>%
    dplyr::left_join(mean_single_expected_df, by = c("rep" = "rep", "gene_symbol" = "gene_symbol")) %>%
    dplyr::group_by(rep) %>%
    group_modify(~ broom::tidy(lm(observed_crispr_single ~ mean_expected_single_crispr, data = .x)))

  # Run the overall linear model
  single_stats <- single_lm %>%
    dplyr::ungroup() %>%
    dplyr::select(term, estimate, rep) %>%
    pivot_wider(
      names_from = term,
      values_from = estimate
    ) %>%
    rename(intercept = "(Intercept)", slope = mean_expected_single_crispr)

  message("Calculating Genetic Interaction scores")

  # Do the linear model adjustments
  gi_calc_adj <- gi_calc_df %>%
    dplyr::left_join(single_stats, by = "rep") %>%
    dplyr::mutate(
      double_target_gi_score = double_crispr_score - (intercept + slope * expected_crispr_double),
      single_target_gi_score_1 = single_crispr_score_1 - (intercept + slope * expected_crispr_single_1),
      single_target_gi_score_2 = single_crispr_score_2 - (intercept + slope * expected_crispr_single_2)
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
      expected_crispr_single_1,
      expected_crispr_single_2,
      expected_crispr_double
    ) %>%
    dplyr::distinct()

  # Store this
  gimap_dataset$results <-
    list(
      overall = single_lm,
      by_target = target_results_df
    )

  return(gimap_dataset)
}


#' Do tests for each replicate --an internal function used by calc_gi() function
#' @description Create results table that has t test p values
#' @param replicate a name of a replicate to filter out from gi_calc_adj
#' @param gi_calc_adj a data.frame with adjusted gi scores
#' @importFrom stats p.adjust t.test wilcox.test
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
