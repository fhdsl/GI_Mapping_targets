#' Calculate CRISPR scores
#' @description This calculates the log fold change for a gimap dataset based on the annotation and metadata provided.
#' @param .data Data can be piped in with tidyverse pipes from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param normalized Default is TRUE meaning that we should expect to look for normalized data in the gimap_dataset.
#' @export
#' @examples \dontrun{
#'
#' gimap_dataset <- get_example_data("gimap")
#'
#' # Highly recommended but not required
#' run_qc(gimap_dataset)
#'
#' gimap_dataset <- gimap_dataset %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA") %>%
#'   gimap_normalize(
#'     timepoints = "day"
#'   ) %>%
#'   calc_crispr()
#'
#' # To see results
#' gimap_dataset$crispr_score
#' }
calc_crispr <- function(.data = NULL,
                        gimap_dataset,
                        normalized = TRUE) {
  # Code adapted from
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/03-filter_and_calculate_LFC.Rmd

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  if (is.null(gimap_dataset$normalized_log_fc) && normalized) {
    stop("No normalized data found in this gimap_dataset. Make sure you have run the gimap_normalize() function or set normalized = FALSE")
  }

  # If data is normalized and normalized is TRUE, we set this to the source data
  if (!is.null(gimap_dataset$normalized_log_fc) && normalized) {
    source_data <- gimap_dataset$normalized_log_fc
  }
  if (gimap_dataset$filtered_data$filter_step_run) {
    pg_ids <- gimap_dataset$filtered_data$metadata_pg_ids$id
  } else {
    pg_ids <- gimap_dataset$metadata$pg_ids
  }

  # If normalized is set to FALSE, we use rawish data but attach annotation
  if (!normalized) {
    source_data <- gimap_dataset$transformed_data$log2_cpm %>%
      as.data.frame() %>%
      dplyr::mutate(pg_ids = pg_ids) %>%
      tidyr::pivot_longer(-pg_ids) %>%
      dplyr::left_join(gimap_dataset$metadata$sample_metadata, by = c("name" = "col_names")) %>%
      dplyr::group_by(timepoints, pg_ids) %>%
      dplyr::summarize(timepoint_avg = mean(value)) %>%
      tidyr::pivot_wider(
        values_from = timepoint_avg,
        names_from = timepoints
      )
  }

  # Calculate medians based on single, double targeting as well as if they are unexpressed control genes
  medians_df <- source_data %>%
    dplyr::group_by(target_type, unexpressed_ctrl_flag) %>%
    dplyr::summarize(median = median(lfc_adj))

  message("Calculating CRISPR score")

  lfc_df <- source_data %>%
    dplyr::left_join(medians_df, by = c("target_type", "unexpressed_ctrl_flag")) %>%
    dplyr::mutate(
      # Since the pgPEN library uses non-targeting controls, we adjusted for the
      # fact that single-targeting pgRNAs generate only two double-strand breaks
      # (1 per allele), whereas the double-targeting pgRNAs generate four DSBs.
      # To do this, we set the median (adjusted) LFC for unexpressed genes of each group to zero.
      crispr_score = lfc_adj - median
    )

  # Get mean control target CRISPR scores -- they will be used for expected calculations
  control_target_df <- lfc_df %>%
    dplyr::filter(target_type == "ctrl_ctrl") %>%
    tidyr::pivot_longer(
      cols = c(gRNA1_seq, gRNA2_seq),
      names_to = "position",
      values_to = "control_gRNA_seq"
    ) %>%
    dplyr::group_by(rep, control_gRNA_seq) %>%
    dplyr::summarize(mean_double_control_crispr = mean(crispr_score, na.rm = TRUE)) %>%
    dplyr::select(rep, control_gRNA_seq, mean_double_control_crispr)

  # Calculate CRISPR scores for single targets
  single_target_df <- lfc_df %>%
    dplyr::filter(target_type %in% c("ctrl_gene", "gene_ctrl")) %>%
    # We will be joining things based on the gRNA sequences so we do some recoding here
    mutate(
      targeting_gRNA_seq = case_when(
        target_type == "gene_ctrl" ~ gRNA1_seq,
        target_type == "ctrl_gene" ~ gRNA2_seq
      ),
      gene_symbol = dplyr::case_when(
        target_type == "gene_ctrl" ~ gene1_symbol,
        target_type == "ctrl_gene" ~ gene2_symbol
      ),
      control_gRNA_seq = case_when(
        target_type == "gene_ctrl" ~ gRNA2_seq,
        target_type == "ctrl_gene" ~ gRNA1_seq
      )
    ) %>%
    group_by(rep, pgRNA_target, targeting_gRNA_seq) %>%
    # Taking the mean of the single target crisprs
    mutate(mean_single_target_crispr = mean(crispr_score)) %>%
    dplyr::select(rep,
      pgRNA_target,
      targeting_gRNA_seq,
      mean_single_target_crispr,
      single_crispr_score = crispr_score
    )

  # Now put it all together into one df
  crispr_df <- lfc_df %>%
    dplyr::filter(target_type %in% c("gene_gene", "ctrl_ctrl")) %>%
    dplyr::select(pg_ids,
                  rep,
                  double_crispr_score = crispr_score,
                  target_type,
                  pgRNA_target,
                  gRNA1_seq,
                  gRNA2_seq,
                  pgRNA_target_double = pgRNA_target
                 ) %>%

    # Join on single target crispr scores
    dplyr::left_join(single_target_df,
      by = c("rep" = "rep", "gRNA1_seq" = "targeting_gRNA_seq"),
      relationship = "many-to-many",
      suffix = c("_double", "_1")
    ) %>%
    dplyr::left_join(single_target_df,
      by = c("rep" = "rep", "gRNA2_seq" = "targeting_gRNA_seq"),
      relationship = "many-to-many",
      suffix = c("_1", "_2")
    ) %>%
    # Join on the control CRISPR
    dplyr::left_join(control_target_df,
      by = c("rep" = "rep", "gRNA1_seq" = "control_gRNA_seq"),
      relationship = "many-to-many",
      suffix = c("", "_1")
    ) %>%
    dplyr::left_join(control_target_df,
      by = c("rep" = "rep", "gRNA2_seq" = "control_gRNA_seq"),
      relationship = "many-to-many",
      suffix = c("", "_2")
    ) %>%
    dplyr::select(
      pg_ids,
      target_type,
      rep,
      double_crispr_score,
      single_crispr_score_1,
      single_crispr_score_2,
      pgRNA_target_double,
      mean_single_target_crispr_1,
      mean_single_target_crispr_2,
      pgRNA1_seq = gRNA1_seq,
      pgRNA2_seq = gRNA2_seq,
      mean_double_control_crispr_1 = mean_double_control_crispr,
      mean_double_control_crispr_2
    )

  # Save at the target level
  gimap_dataset$crispr_score <- crispr_df

  return(gimap_dataset)
}
