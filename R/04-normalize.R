#' Normalize Log fold changes
#' @description This calculates the log fold change for a gimap dataset based on the annotation and metadata provided.
#' @param .data Data can be piped in with a tidyverse pipe from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param control_name A name that specifies the data either in the treatments column that should be used as the control. This could be the Day 0 of treatment or an untreated sample.
#'  For timepoints testing it will be assumed that the mininmum timepoint is the control.
#' @param timepoints Specifies the column name of the metadata set up in `$metadata$sample_metadata`
#' that has a factor that represents the timepoints. Timepoints will be made into three categories:
#' plasmid for the earliest time point, early for all middle timepoints and late for the latest timepoints.
#' The late timepoints will be the focus for the calculations. The column used for timepoints must be numeric or at least ordinal.
#' @param treatments Specifies the column name of the metadata set up in `$metadata$sample_metadata`
#' that has a factor that represents column that specifies the treatment applied to each. The replicates will be kept collapsed to an average.
#' @param num_ids_wo_annot default is 20; the number of pgRNA IDs to display to console if they don't have corresponding annotation data;
#' ff there are more IDs without annotation data than this number, the output will be sent to a file rather than the console.
#' @param rm_ids_wo_annot default is TRUE; whether or not to filter out pgRNA IDs from the input dataset that don't have corresponding annotation data available
#' @exports
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
#'     timepoints = "day",
#'   )
#'
#' # To see results
#' gimap_dataset$normalized_log_fc
#' }
gimap_normalize <- function(.data = NULL,
                            gimap_dataset,
                            timepoints = NULL,
                            treatments = NULL,
                            control = NULL,
                            num_ids_wo_annot = 20,
                            rm_ids_wo_annot = TRUE) {

  # Code adapted from
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/03-filter_and_calculate_LFC.Rmd

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  # Based on log fold change calculations and other handling will go based on the code in:
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/03-filter_and_calculate_LFC.Rmd

  if (gimap_dataset$filtered_data$filter_step_run) {
    dataset <- gimap_dataset$filtered_data$transformed_log2_cpm
    pg_ids <- gimap_dataset$filtered_data$metadata_pg_ids$id
  } else {
    dataset <- gimap_dataset$transformed_data$log2_cpm
    pg_ids <- gimap_dataset$metadata$pg_ids
  }

  # Doing some reshaping to get one handy data frame
  lfc_df <- dataset %>%
    as.data.frame() %>%
    dplyr::mutate(pg_ids = pg_ids) %>%
    tidyr::pivot_longer(-pg_ids, values_to = "log2_cpm") %>%
    # Adding on metadata
    dplyr::left_join(gimap_dataset$metadata$sample_metadata, by = c("name" = "col_names")) %>%
    dplyr::select(-name)


  ### IF WE HAVE TREATMENTS
  if (!is.null(treatments)) {
    if (!(treatments %in% colnames(gimap_dataset$metadata$sample_metadata))) {
      stop("The column name specified for 'treatments' does not exist in gimap_dataset$metadata$sample_metadata")
    }
    # Rename and recode the timepoints variable
    gimap_dataset$metadata$sample_metadata <- gimap_dataset$metadata$sample_metadata %>%
      dplyr::rename(treatments = all_of(treatments))

    # Collapse reps
    lfc_df <- lfc_df %>%
      tidyr::pivot_wider(values_from = "log2_cpm",
                         names_from = treatments,
                         values_fn = mean)

    if (!is.null(control_name)) {
      if (!(control_name %in% colnames(lfc_df))) {
        stop("There are no samples with the label of 'control' in the treatments column")
      }
      gimap_dataset$metadata$sample_metadata <- gimap_dataset$metadata$sample_metadata %>%
        dplyr::rename(control = control_name)
    }
  }

  ### IF WE HAVE TIMEPOINTS
  if (!is.null(timepoints)) {
    if (!(timepoints %in% colnames(gimap_dataset$metadata$sample_metadata))) {
      stop("The column name specified for 'timepoints' does not exist in gimap_dataset$metadata$sample_metadata")
    }

    # Rename and recode the timepoints variable
    gimap_dataset$metadata$sample_metadata <- gimap_dataset$metadata$sample_metadata %>%
      dplyr::rename(timepoints = all_of(timepoints)) %>%
      # Note that timepoints are extablished as three categories: control, early, or late.
      dplyr::mutate(timepoints = dplyr::case_when(
        timepoints == min(timepoints) ~ "control",
        timepoints == max(timepoints) ~ "late",
        TRUE ~ "early"
      ))

    # Collapse reps
    lfc_df <- lfc_df %>%
      tidyr::pivot_wider(values_from = "log2_cpm",
                         names_from = treatments,
                         values_fn = mean)
  }

  ### Stop if no annotations
  if (is.null(gimap_dataset$annotation)) {
    stop(
      "No annotations are stored in this gimap_dataset, annotations are needed so we know what genes should be used as controls.",
      "Please run gimap_annotate() function on your gimap_dataset and then retry this function."
    )
  }

  message("Normalizing Log Fold Change")


  ##### Checking for missing ids
  missing_ids <- data.frame(
    missing_ids = setdiff(lfc_df$pg_ids, gimap_dataset$annotation$pgRNA_id)
  )

  if ((nrow(missing_ids) > 0) & (nrow(missing_ids) < num_ids_wo_annot)){
    message("The following ", nrow(missing_ids), " IDs were not found in the annotation data: \n", paste0(missing_ids, collapse = ", "))
  } else {
    missing_ids_file <- file.path("missing_ids_file.csv")
    readr::write_csv(missing_ids, missing_ids_file)
  }

  if ((nrow(missing_ids) > 0) & (rm_ids_wo_annot == TRUE)){
    lfc_df <- lfc_df %>%
      filter(!pg_ids %in% missing_ids)
    message("The input data for the IDs which were not found in the annotation data has been filtered out and will not be included in the analysis output.")
  } else{
    message("The input data for the IDs which were not found in the annotation data will be kept throughout the analysis, but any data from the annotation won't be available for them.")
  }

  ###################### Subtract the control column (so either day 0 or pretreatment)


  comparison_df <-  lfc_df %>%
      dplyr::mutate_at(dplyr::vars(!c(pg_ids, control)), ~.x - control) %>%
      dplyr::select(-control)  %>%
      dplyr::left_join(gimap_dataset$annotation, by = c("pg_ids" = "pgRNA_id"))

  ########################### Perform adjustments #############################

  ### Calculate medians
  neg_control_median_df <- comparison_df %>%
    dplyr::filter(norm_ctrl_flag == "negative_control") %>%
    dplyr::select(pg_ids, dplyr::starts_with("late"))

  # Find a median for each rep, so apply across columns
  neg_control_median <- apply(neg_control_median_df[, -1], 2, median)

  # First and second adjustments to LFC
  lfc_df_adj <- comparison_df %>%
    #subtract the correct replicate negative control median from the late vs plasmid difference
    mutate(across(names(neg_control_median), ~ . - neg_control_median[cur_column()])) %>%
    tidyr::pivot_longer(dplyr::starts_with("late"),
                        names_to = "rep",
                        values_to = "lfc_adj1")  %>%
    group_by(rep) %>%
    dplyr::mutate(
      # Then, divide by the median of negative controls (double non-targeting) minus
      # median of positive controls (targeting 1 essential gene).
      # This will effectively set the median of the positive controls (essential genes) to -1.
      lfc_adj = lfc_adj1 / (median(lfc_adj1[norm_ctrl_flag == "negative_control"], na.rm = TRUE) - median(lfc_adj1[norm_ctrl_flag == "positive_control"]))
    ) %>%
    ungroup()

  # Save this at the construct level
  gimap_dataset$normalized_log_fc <- lfc_df_adj

  return(gimap_dataset)
}
