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
#'       timepoints = "day"
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

  # Calculate means by targets, collapsing it if the other gene is a control
  double_mean_df <- gi_calc_df %>%
    dplyr::group_by(rep, pgRNA_target_double) %>%
    dplyr::summarize(
      mean_expected_double_crispr = mean(expected_crispr_double, na.rm = TRUE),
      mean_observed_double_crispr = mean(double_crispr_score, na.rm = TRUE),
    )

  # Calculating mean crisprs
  reshaped_single_df <- gi_calc_df %>%
    # reshaping the data a bit so we can do the math later
    dplyr::select(
      rep,
      pgRNA_target_double,
      pgRNA_target_single,
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
    ) %>%
    # Because of the pivoting we have some duplicates
    # So we want to get down to non-redundant data
    dplyr::distinct()


  ## calculate mean CRISPR score of single-targeting pgRNAs containing the same targeting
  ## sgRNA sequence but different control sgRNA sequences
  single_mean_df <- reshaped_single_df %>%
    dplyr::group_by(rep, pgRNA_target_single) %>%
    dplyr::summarize(
      mean_expected_single_crispr = mean(expected_crispr_single, na.rm = TRUE),
      mean_observed_single_crispr = mean(observed_crispr_single, na.rm = TRUE)
    )

  # Calculate the linear models from this
  single_lm_df <- single_mean_df %>%
    group_modify(~ broom::tidy(lm(mean_observed_single_crispr ~ mean_expected_single_crispr, data = .x)))

  double_lm_df <- double_mean_df %>%
    group_modify(~ broom::tidy(lm(mean_observed_double_crispr ~ mean_expected_double_crispr, data = .x)))

  # Run the overall linear model for single targets
  single_stats <- single_lm_df %>%
    dplyr::select(term, estimate, rep) %>%
    pivot_wider(
      names_from = term,
      values_from = estimate
    ) %>%
    rename(intercept = "(Intercept)", slope = mean_expected_single_crispr)

  # Run the overall linear model for double targets
  double_stats <- double_lm_df %>%
    dplyr::select(term, estimate, rep) %>%
    pivot_wider(
      names_from = term,
      values_from = estimate
    ) %>%
    rename(intercept = "(Intercept)", slope = mean_expected_double_crispr)

  message("Calculating Genetic Interaction scores")

  # Do the linear model adjustments
  gi_calc_single <- single_mean_df %>%
    dplyr::left_join(single_stats, by = "rep") %>%
    dplyr::mutate(
      single_target_gi_score = mean_observed_single_crispr - (intercept + slope * mean_expected_single_crispr),
    )

  # Do the linear model adjustments but don't collapse double
  gi_calc_double <- gi_calc_df %>%
    # Tack on the mean expected double crispr
    dplyr::left_join(dplyr::select(double_mean_df, mean_expected_double_crispr, rep, pgRNA_target_double),
                     by = c("rep", "pgRNA_target_double")) %>%
    # Using the single target's linear model here
    dplyr::left_join(single_stats, by = "rep") %>%
    dplyr::mutate(
      double_target_gi_score = double_crispr_score - (intercept + slope * mean_expected_double_crispr)
    )

  # Which replicates we have?
  replicates <- unique(gi_calc_df$rep)

  # TODO: this is an lapply inside an lapply at the end of the day. Not ideal
  target_results <- lapply(replicates, function(replicate) {
    gimap_rep_stats(replicate,
                    gi_calc_single = gi_calc_single,
                    gi_calc_double = gi_calc_double)
  }
  )

  # Bring over the replicate names
  names(target_results) <- replicates

  # Turn into a data.frame
  target_results_df <- dplyr::bind_rows(target_results, .id = "replicate") %>%
    tidyr::pivot_wider(
      names_from = replicate,
      values_from = c(
        p_val,
        fdr
      )
    )

  # Store the useful bits
  gimap_dataset$gi_scores <- gi_calc_double

  # Store this
  gimap_dataset$results <-
    list(
      overall = single_lm_df,
      by_target = target_results_df
    )

  return(gimap_dataset)
}


#' Do tests for each replicate --an internal function used by calc_gi() function
#' @description Create results table that has t test p values
#' @param replicate a name of a replicate to filter out from gi_calc_adj
#' @param gi_calc_adj a data.frame with adjusted gi scores
#' @importFrom stats p.adjust t.test wilcox.test
gimap_rep_stats <- function(replicate, gi_calc_double,  gi_calc_single) {
  ## get a vector of GI scores for all single-targeting ("control") pgRNAs for each rep
  ## get double-targeting pgRNAs for this rep, do a t-test to compare the double-

  ## targeting GI scores for each paralog pair to the control vector

  ## adjust for multiple testing using the Benjamini-Hochberg method

  per_rep_stats_double <- gi_calc_double %>%
    dplyr::filter(rep == replicate) %>%
    dplyr::ungroup()

  per_rep_stats_single <- gi_calc_single %>%
    dplyr::filter(rep == replicate)

  rep_gi_scores <- per_rep_stats_double %>%
    group_by(pgRNA_target_double) %>%
    mutate(p_val = t.test(x = per_rep_stats_single$mean_observed_single_crispr,
                          y = double_target_gi_score,
                          paired = FALSE)$p.value)


  ## adjust for multiple testing using the Benjamini-Hochberg method
  d_p_val <- rep_gi_scores %>%
    dplyr::select(pgRNA_target_double, p_val) %>%
    arrange(p_val) %>%
    distinct(p_val, .keep_all = TRUE)

  fdr_vals <- p.adjust(d_p_val$p_val, method = "BH")

  d_fdr <- tibble("fdr" = fdr_vals) %>%
     bind_cols(d_p_val) %>%
    dplyr::select(pgRNA_target_double, p_val, fdr)

  return(d_fdr)
}
