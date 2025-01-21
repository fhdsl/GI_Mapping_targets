utils::globalVariables(c(
  "pg_ids", "plot_theme()", "negative_control", "positive_control", "mean_observed_cs", "timepoints", "value", "timepoint_avg", "target_type",
  "unexpressed_ctrl_flag", "median", "lfc_adj", "median", "gRNA1_seq", "gRNA2_seq",
  "control_gRNA_seq", "crispr_score", "pgRNA_target", "mean_double_control_crispr",
  "pgRNA_target", "targeting_gRNA_seq", "mean_single_crispr", "double_crispr_score",
  "single_crispr_score_1", "single_crispr_score_2", "pgRNA_target_double", "mean_single_crispr_1",
  "mean_single_crispr_2", "mean_double_control_crispr_2",
  "expected_crispr", "term", "estimate", "mean_expected_crispr", "intercept", "slope",
  "p_val_ttest", "p_val_wil", "fde_vals_ttest", "fdr_vals_wil", "double_gi_score",
  "single_gi_score_1", "single_gi_score_2", "gene", "DepMap_ID",
  "gene1_symbol", "gene2_symbol", "expressed_flag", "norm_ctrl_flag", "bool_vals",
  "filter_name", "counts", "numzero", "name", "value", "lfc_plasmid_vs_late", "lfc_adj",
  "double_gi_score", "count_normalized", "construct",
  "filterFlag", "plasmid_log2_cpm", "log2_cpm", "gene_symbol", "gene_symbol_1", "gene_symbol_2",
  "mean_double_control_crispr_1", "expected_crispr_double", "expected_crispr_single_1",
  "expected_crispr_single_2", "fdr_vals_ttest", "read_table", "stripped_cell_line_name",
  "comparison", ".", "col_names", "lfc_adj1", "t.test", "wilcox.test", "p.adjust",
  "cor", "quantile", "var", "browseURL", "single_crispr", "mean_single_crispr_2",
  "mean_single_crispr_1", "expected_single_crispr", "double_crispr", "double_gi_score",
  "fdr", "lfc", "mean_expected_cs", "mean_gi_score", "mean_single_crispr",
  "expected_double_crispr", "p_val", "single_gi_score", "Rank", "broad_target_type",
  "logfdr", "pointColor"
))


#' Returns example data for package
#' @description This function loads and returns example data for the packagae. Which dataset is returned must be specified
#' @param which_data options are "count" or "meta"; specifies which example dataset should be returned
#' @export
#' @examples \dontrun{
#'
#' gimap_dataset <- get_example_data("gimap")
#' }
get_example_data <- function(which_data) {
  if (which_data == "count") {
    file <- list.files(
      pattern = "PP_pgPEN_HeLa_counts.txt",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_tsv(file))
  } else if (which_data == "count_treatment") {
    file <- list.files(
      pattern = "counts_pgPEN_PC9_example.tsv",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_tsv(file))
  } else if (which_data == "meta") {
    file <- list.files(
      pattern = "pgRNA_ID_pgPEN_library_comp.csv",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_csv(file, skip = 1))
  } else if (which_data == "gimap") {
    file <- list.files(
      pattern = "gimap_dataset.RDS",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_rds(file))
  } else if (which_data == "annotation") {
    file <- list.files(
      pattern = "pgPEN_annotations.txt",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_tsv(file, show_col_types = FALSE))
  } else if (which_data == "final_treatment") {
    file <- list.files(
      pattern = "gimap_dataset_final_treatment.RDS",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_rds(file))
  } else if (which_data == "final_treatment") {
    file <- list.files(
      pattern = "gimap_dataset_final_treatment.RDS",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_rds(file))
  } else {
    stop("Specification for `which_data` not understood; Need to use 'gimap', 'count', 'meta', or 'annotation' ")
  }
}


#' Get file path to an default credentials RDS
#' @export
#' @return Returns the file path to folder where the example data is stored
example_data_folder <- function() {
  file <- list.files(
    pattern = "example_data.md",
    recursive = TRUE,
    system.file("extdata", package = "gimap"),
    full.names = TRUE
  )
  dirname(file)
}

# This function sets up the example count data
save_example_data <- function() {
  example_data <- get_example_data("count")

  example_pg_metadata <- get_example_data("meta")

  example_counts <- example_data %>%
    dplyr::select(c("Day00_RepA", "Day05_RepA", "Day22_RepA", "Day22_RepB", "Day22_RepC")) %>%
    as.matrix()

  example_pg_id <- example_data %>%
    dplyr::select("id")

  example_pg_metadata <- example_data %>%
    dplyr::select(c("id", "seq_1", "seq_2"))

  example_sample_metadata <- data.frame(
    col_names = c("Day00_RepA", "Day05_RepA", "Day22_RepA", "Day22_RepB", "Day22_RepC"),
    day = as.numeric(c("0", "5", "22", "22", "22")),
    rep = as.factor(c("RepA", "RepA", "RepA", "RepB", "RepC"))
  )

  gimap_dataset <- setup_data(
    counts = example_counts,
    pg_ids = example_pg_id,
    sample_metadata = example_sample_metadata
  )

  example_folder <- list.files(
    pattern = "PP_pgPEN_HeLa_counts.txt",
    recursive = TRUE,
    system.file("extdata", package = "gimap"),
    full.names = TRUE
  )

  saveRDS(gimap_dataset, file.path(dirname(example_folder), "gimap_dataset.RDS"))
}

plot_options <- function() {
  list(theme_bw(base_size = 12))
}
