test_that("Test normalization", {
  testthat::skip_on_cran()
  gimap_dataset <- get_example_data("gimap") %>%
    gimap_filter() %>%
    gimap_annotate(cell_line = "HELA") %>%
    gimap_normalize(
      timepoints = "day"
    )

  # make sure the important columns are there
  testthat::expect_true(
    all(c("target_type", "lfc", "rep", "crispr_score", "unexpressed_ctrl_flag")
    %in% colnames(gimap_dataset$normalized_log_fc))
  )

  neg_controls <- gimap_dataset$normalized_log_fc %>%
    dplyr::filter(norm_ctrl_flag == "negative_control") %>%
    dplyr::group_by(rep) %>%
    dplyr::summarize(neg_ctrl_med = median(crispr_score)) %>%
    dplyr::pull(neg_ctrl_med)

  # We expect negative controls to be now equal to 0
  testthat::expect_equal(neg_controls[2:4], c(0, 0, 0))

  pos_controls <- gimap_dataset$normalized_log_fc %>%
    dplyr::filter(norm_ctrl_flag == "positive_control") %>%
    dplyr::group_by(rep) %>%
    dplyr::summarize(pos_ctrl_med = median(crispr_score)) %>%
    dplyr::pull(pos_ctrl_med)

  # We expect positive controls to be now equal to -1
  testthat::expect_equal(
    round(pos_controls[-1]),
    round(c(-1, -1, -1))
  )
})
