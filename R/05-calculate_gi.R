#' Calculate Genetic Interaction scores
#' @description Create results table that has CRISPR scores, Wilcoxon rank-sum
#' test and t tests.
#' The output of the `gimap` package is genetic interaction scores which _is the
#' distance between the observed CRISPR score and the expected CRISPR score._
#' The expected CRISPR scores are what we expect for the CRISPR values should
#' two genes be unrelated to each other. The further away an observed CRISPR
#' scoreis from its expected the more we suspect genetic interaction.
#' This can be true in a positive way (a CRISPR knockout pair caused more cell
#' proliferation than expected) or in a negative way (a CRISPR knockout pair
#' caused more cell lethality than expected).
#'
#' The genetic interaction scores are based on a linear model calculated for
#' each sample where `observed_crispr_single` is the outcome variable and
#' `expected_crispr_single` is the predictor variable.
#' For each sample: lm(observed_crispr_single ~ expected_crispr_single)
#'
#' Using `y = mx+b`, we can fill in the following values:
#' * `y` = observed CRISPR score
#' * `x` = expected CRISPR score
#' * `m` = slope from linear model for this sample
#' * `b` = intercept from linear model for this sample
#'
#' The intercept and slope from this linear model are used to adjust the CRISPR
#' scores for each sample:
#' single target gi score =
#'   observed single crispr - (intercept + slope * expected single crispr)
#' double_target_gi_score =
#'   double crispr score - (intercept + slope * expected double crispr)
#' These single and double target genetic interaction scores are calculated at
#' the construct level and are then summarized using a t-test to see if the the
#' distribution of the set of double targeting constructs is significantly
#' different than the overall distribution single targeting constructs.
#' After multiple testing correction, FDR values are reported.
#' Low FDR value for a double construct means high suspicion of paralogs.
#' @param .data Data can be piped in with tidyverse pipes from function to
#' function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the
#' `setup_data()` function.
#' @param use_lfc Should Log fold change be used to calculate GI scores instead
#' of CRISPR scores? If you do not have negative controls or CRISPR scores you
#' will need to set this to TRUE.
#' @param stats_by_rep Should statistics be calculated per rep or should replicates
#' be collapsed?
#' @return A gimap dataset with statistics and genetic interaction scores
#' calculated. Overall results in the returned object can be obtained using
#' gimap_dataset$overall_results Whereas target level genetic interaction
#' scores can be retrieved using `gimap_dataset$gi_scores`.
#' @import dplyr
#' @importFrom stats lm
#' @export
#' @examples \donttest{
#'
#' gimap_dataset <- get_example_data("gimap",
#'                                   data_dir = tempdir()) %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA",
#'                  annot_dir = tempdir()) %>%
#'   gimap_normalize(
#'     timepoints = "day",
#'     missing_ids_file = tempfile()
#'   ) %>%
#'   calc_gi()
#'
#' saveRDS(gimap_dataset, file.path(tempdir(), "gimap_dataset_final.RDS"))
#' }
calc_gi <- function(.data = NULL,
                    gimap_dataset,
                    use_lfc = FALSE,
                    stats_by_rep = FALSE) {
  # Summary the calculation
  # single_target_crispr_1 = geneA_nt1, geneA_nt2...
  # single_target_crispr_2 = nt1_geneB, nt2_geneB...
  # double_crispr_score = geneA_geneBpg1, geneA_geneBpg2...

  # mean_double_control_crispr = mean for the same control sequence

  # expected_crispr_double=single_target_crispr_1 + single_target_crispr_2
  # expected_crispr_single_1=single_target_crispr_1 + mean_double_control_crispr
  # expected_crispr_single_2=single_target_crispr_2 + mean_double_control_crispr

  # linear model is lm(observed_single_crispr ~ expected_single_crispr)

  # single_target_gi_score = observed_single_crispr -
  #      (intercept + slope * expected_single_crispr)
  # double_target_gi_score = double_crispr_score -
  #      (intercept + slope * expected_double_crispr)

  # Code adapted from
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/
  # 04-calculate_GI_scores.Rmd

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) {
    stop(
      "This function only works",
      " with gimap_dataset objects which can be made with the",
      " setup_data() function."
    )
  }

  if (is.null(gimap_dataset$normalized_log_fc)) {
    stop(
      "This function only works",
      "with already normalized gimap_dataset objects",
      "which can be done with the gimap_normalize() function."
    )
  }
  lfc_adj <- gimap_dataset$normalized_log_fc

  if (!"crispr_score" %in% colnames(lfc_adj)) {

   if (!use_lfc) {
      stop("No CRISPR scores found, you can set use_lfc = TRUE to calculate
         genetic interaction scores using log fold change")
    }
  }
  if (use_lfc) {
    lfc_adj$crispr_score <- lfc_adj$lfc
    message("Calculating Genetic Interaction scores using LFC")
  }
  # Get mean control target CRISPR scores -- they will be used for expected
  # calculations
  control_target_df <- lfc_adj %>%
    dplyr::filter(target_type == "ctrl_ctrl") %>%
    tidyr::pivot_longer(
      cols = c(gRNA1_seq, gRNA2_seq),
      names_to = "position",
      values_to = "control_gRNA_seq"
    ) %>%
    # If there's the same control sequence, and rep
    dplyr::group_by(rep, pgRNA_target, control_gRNA_seq, norm_ctrl_flag) %>%
    # Then take the mean for when controls have the same sequence
    dplyr::summarize(
      mean_double_control_crispr =
        mean(crispr_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(
      rep,
      pgRNA_target,
      control_gRNA_seq,
      mean_double_control_crispr,
      norm_ctrl_flag,
    )
  # This means we have a mean double control crispr for each rep and
  # control sequence

  # Calculate CRISPR scores for single targets
  single_crispr_df <- lfc_adj %>%
    dplyr::filter(target_type %in% c("ctrl_gene", "gene_ctrl")) %>%
    # We will be joining things based on the gRNA sequences so
    # we do some recoding here
    mutate(
      targeting_gRNA_seq = case_when(
        target_type == "gene_ctrl" ~ gRNA1_seq,
        target_type == "ctrl_gene" ~ gRNA2_seq
      ),
      control_gRNA_seq = case_when(
        target_type == "gene_ctrl" ~ gRNA2_seq,
        target_type == "ctrl_gene" ~ gRNA1_seq
      ),
      gene_symbol = dplyr::case_when(
        target_type == "gene_ctrl" ~ gene1_symbol,
        target_type == "ctrl_gene" ~ gene2_symbol
      ),
    ) %>%
    dplyr::left_join(control_target_df,
      by = c("rep" = "rep", "control_gRNA_seq" = "control_gRNA_seq"),
      suffix = c("", "_control")
    ) %>%
    group_by(rep, pgRNA_target, targeting_gRNA_seq) %>%
    # Taking the mean of the single target crisprs that have the same
    # targeting sequence
    mutate(mean_single_crispr = mean(crispr_score, na.rm = TRUE)) %>%
    dplyr::select(rep,
      pgRNA_target,
      targeting_gRNA_seq,
      control_gRNA_seq,
      single_crispr = crispr_score,
      mean_single_crispr,
      mean_double_control_crispr,
      norm_ctrl_flag
    ) %>%
    dplyr::distinct()



  if (nrow(control_target_df) > 0) {
    single_crispr_df <- single_crispr_df %>%
      ## calculate expected double-targeting GI score by summing the two mean
      ## single-targeting
      ## CRISPR scores for that paralog pair
      dplyr::mutate(
        expected_single_crispr = single_crispr + mean_double_control_crispr,
      )
  } else {
    single_crispr_df <- single_crispr_df %>%
      ## calculate expected double-targeting GI score by summing the two mean
      ## single-targeting
      ## CRISPR scores for that paralog pair
      dplyr::mutate(
        expected_single_crispr = single_crispr,
      )
  }

  # Calculate expected
  expected_single_crispr_df <- single_crispr_df %>%
    dplyr::select(rep, pgRNA_target, targeting_gRNA_seq, mean_single_crispr) %>%
    dplyr::distinct()

  # Now put it all together into one df
  double_crispr_df <- lfc_adj %>%
    dplyr::filter(target_type == "gene_gene") %>%
    dplyr::select(
      pg_ids,
      rep,
      crispr_score,
      gRNA1_seq,
      gRNA2_seq,
      pgRNA_target,
      norm_ctrl_flag
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(expected_single_crispr_df,
      by = c("rep" = "rep", "gRNA1_seq" = "targeting_gRNA_seq"),
      suffix = c("", "_1")
    ) %>%
    dplyr::left_join(expected_single_crispr_df,
      by = c("rep" = "rep", "gRNA2_seq" = "targeting_gRNA_seq"),
      suffix = c("", "_2")
    ) %>%
    dplyr::select(pg_ids,
      rep,
      double_crispr = crispr_score,
      norm_ctrl_flag,
      gRNA1_seq,
      gRNA2_seq,
      pgRNA_target,
      mean_single_crispr_1 = mean_single_crispr,
      mean_single_crispr_2
    ) %>%
    dplyr::mutate(
      expected_double_crispr = mean_single_crispr_1 + mean_single_crispr_2
    )

  # Calculate the linear models from this
  single_lm_df <- single_crispr_df %>%
    dplyr::filter(!is.na(single_crispr), !is.na(expected_single_crispr)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(rep) %>%
    group_modify(~ broom::tidy(
      lm(single_crispr ~ expected_single_crispr, data = .x)
    )) %>%
    dplyr::select(term, estimate, rep) %>%
    pivot_wider(
      names_from = term,
      values_from = estimate
    ) %>%
    rename(intercept = "(Intercept)", slope = expected_single_crispr)

  message("Calculating Genetic Interaction scores")

  # Do the linear model adjustments
  gi_calc_single <- single_crispr_df %>%
    dplyr::left_join(single_lm_df, by = "rep") %>%
    dplyr::mutate(
      single_gi_score = single_crispr -
        (intercept + slope * expected_single_crispr)
    )

  # Do the linear model adjustments but don't collapse double
  gi_calc_double <- double_crispr_df %>%
    # Using the single target's linear model here
    dplyr::left_join(single_lm_df, by = "rep") %>%
    dplyr::mutate(
      double_gi_score = double_crispr -
        (intercept + slope * expected_double_crispr)
    )

  # Which replicates we have?
  replicates <- unique(gi_calc_double$rep)

  if (stats_by_rep) {
    target_results <- lapply(replicates, function(replicate) {
      gimap_stats(replicate,
        gi_calc_single = gi_calc_single,
        gi_calc_double = gi_calc_double
      )
    })

    # Bring over the replicate names
    names(target_results) <- replicates

    # Turn into a data.frame
    target_results_df <- dplyr::bind_rows(target_results, .id = "rep")

    ## Clean up the data and get some means
    gi_calc_double <- gi_calc_double %>%
      dplyr::group_by(
        rep,
        pgRNA_target
      )
  } else {
    target_results_df <- gimap_stats(
      gi_calc_single = gi_calc_single,
      gi_calc_double = gi_calc_double
    )
    ## Clean up the data and get some means
    gi_calc_double <- gi_calc_double %>%
      dplyr::group_by(
        pgRNA_target
      )
  }

  ## Clean up the data and get some means
  gi_calc_double <- gi_calc_double %>%
    dplyr::summarize(
      mean_expected_cs = mean(expected_double_crispr, na.rm = TRUE),
      mean_observed_cs = mean(double_crispr, na.rm = TRUE),
      mean_gi_score = mean(double_gi_score, na.rm = TRUE)
    )

  if (stats_by_rep) {
    gi_calc_double <- gi_calc_double %>%
      # Collapse to just stats and don't care about pg_ids anymore
      dplyr::select(
        rep,
        pgRNA_target,
        mean_expected_cs,
        mean_observed_cs,
        mean_gi_score
      )
  } else {
    gi_calc_double <- gi_calc_double %>%
      # Collapse to just stats and don't care about pg_ids anymore
      dplyr::select(
        pgRNA_target,
        mean_expected_cs,
        mean_observed_cs,
        mean_gi_score
      )
  }
  gi_calc_double <- gi_calc_double %>%
    dplyr::mutate(target_type = "gene_gene") %>%
    dplyr::distinct()

  if (stats_by_rep) {
    # Same kind of reformatting for single
    gi_calc_single <- gi_calc_single %>%
      dplyr::group_by(
        rep,
        pgRNA_target
      )
  } else {
    # Same kind of reformatting for single
    gi_calc_single <- gi_calc_single %>%
      dplyr::group_by(
        pgRNA_target
      )
  }
  gi_calc_single <- gi_calc_single %>%
    dplyr::summarize(
      mean_expected_cs = mean(expected_single_crispr, na.rm = TRUE),
      mean_observed_cs = mean(single_crispr, na.rm = TRUE),
      mean_gi_score = mean(single_gi_score, na.rm = TRUE)
    ) %>%
    dplyr::mutate(target_type = dplyr::case_when(
      grepl("^ctrl_*", pgRNA_target) ~ "ctrl_gene",
      grepl("*_ctrl$", pgRNA_target) ~ "gene_ctrl"
    ))

    if (stats_by_rep) {
      gi_calc_single <- gi_calc_single %>%
        dplyr::select(
          rep,
          target_type,
          pgRNA_target,
          mean_expected_cs,
          mean_observed_cs,
          mean_gi_score
        ) %>%
        dplyr::distinct()
    } else {
      gi_calc_single <- gi_calc_single %>%
        dplyr::select(
          target_type,
          pgRNA_target,
          mean_expected_cs,
          mean_observed_cs,
          mean_gi_score
        ) %>%
        dplyr::distinct()
    }

  all_gi_scores <- dplyr::bind_rows(gi_calc_double, gi_calc_single)


  if (stats_by_rep) {
  # Add on test results
  all_gi_scores <- all_gi_scores %>%
    dplyr::left_join(target_results_df,
      by = c("pgRNA_target", "rep")
    )
  } else {
    # Add on test results
    all_gi_scores <- all_gi_scores %>%
      dplyr::left_join(target_results_df,
                       by = c("pgRNA_target")
      )
  }

  # Store the useful bits
  gimap_dataset$gi_scores <- all_gi_scores

  # Store this
  gimap_dataset$overall_results <- single_lm_df

  means_by_rep <- dplyr::bind_rows(
    dplyr::select(single_crispr_df,
                  rep, pgRNA_target, norm_ctrl_flag, mean_score = mean_single_crispr,
                  seq = targeting_gRNA_seq),
    dplyr::select(double_crispr_df,
                  rep, pgRNA_target, norm_ctrl_flag, mean_score = double_crispr),
    dplyr::select(control_target_df,
                  rep, pgRNA_target, norm_ctrl_flag, mean_score = mean_double_control_crispr,
                  seq = control_gRNA_seq))

  if (!use_lfc) {
    # Save at the target level but call them crispr scores
    gimap_dataset$crispr_score$single_crispr_score <- single_crispr_df
    gimap_dataset$crispr_score$double_crispr_score <- double_crispr_df
    gimap_dataset$crispr_score$neg_control_crispr <- control_target_df
    gimap_dataset$crispr_score$means_by_rep <- means_by_rep
  } else {
    # Save at the target level but they aren't crispr scores
    colnames(single_crispr_df) <- gsub("crispr", "lfc", colnames(single_crispr_df))
    colnames(double_crispr_df) <- gsub("crispr", "lfc", colnames(double_crispr_df))
    colnames(control_target_df) <- gsub("crispr", "lfc", colnames(control_target_df))

    gimap_dataset$lfc$single_lfc <- single_crispr_df
    gimap_dataset$lfc$double_lfc <- double_crispr_df
    gimap_dataset$lfc$neg_control_crispr <- control_target_df

    # The final results also aren't crispr scores
    gimap_dataset$gi_scores <- gimap_dataset$gi_scores %>%
      dplyr::rename(
        mean_expected_lfc = mean_expected_cs,
        mean_observed_lfc = mean_observed_cs
      )
    gimap_dataset$crispr_score$means_by_rep <- means_by_rep
  }

  return(gimap_dataset)
}


#' Do tests --an internal function used by calc_gi() function
#' @description Create results table that has t test p values
#' @param gi_calc_single a data.frame with adjusted single gi scores
#' @param gi_calc_double a data.frame with adjusted double gi scores
#' @param replicate a name of a replicate to filter out from gi_calc_adj Optional
#' @importFrom stats p.adjust t.test wilcox.test
gimap_stats <- function(gi_calc_double, gi_calc_single, replicate = NULL) {
  ## get a vector of GI scores for all single-targeting ("control") pgRNAs
  ## for each rep
  ## get double-targeting pgRNAs for this rep, do a t-test to compare the
  ## double-targeting GI scores for each paralog pair to the control vector
  ## adjust for multiple testing using the Benjamini-Hochberg method

  if (!is.null(replicate)) {
    gi_calc_double <- gi_calc_double %>%
      dplyr::filter(rep == replicate)

    gi_calc_single <- gi_calc_single %>%
      dplyr::filter(rep == replicate)
  }

  gi_scores <- gi_calc_double %>%
    group_by(pgRNA_target) %>%
    mutate(p_val = t.test(
      x = gi_calc_single$single_gi_score,
      y = double_gi_score, # all 16 construct guides here
      paired = FALSE
    )$p.value)


  ## adjust for multiple testing using the Benjamini-Hochberg method
  d_p_val <- gi_scores %>%
    dplyr::select(pgRNA_target, p_val) %>%
    arrange(p_val) %>%
    distinct(p_val, .keep_all = TRUE)

  fdr_vals <- p.adjust(d_p_val$p_val, method = "BH")

  d_fdr <- tibble("fdr" = fdr_vals) %>%
    bind_cols(d_p_val) %>%
    dplyr::select(pgRNA_target, p_val, fdr)

  return(d_fdr)
}
