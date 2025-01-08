test_that("Test normalization", {

   gimap_dataset <- get_example_data("gimap") %>%
     gimap_filter() %>%
     gimap_annotate(cell_line = "HELA") %>%
     gimap_normalize(
       timepoints = "day"
     )

  # make sure the important columns are there
  testthat::expect_true(all(c("target_type", "lfc_adj", "rep") %in% colnames(gimap_dataset$normalized_log_fc)))

  neg_controls <- gimap_dataset$normalized_log_fc %>%
    dplyr::filter(target_type == "ctrl_ctrl") %>%
    dplyr::group_by(rep) %>%
    dplyr::summarize(neg_ctrl_med = median(lfc_adj)) %>%
    dplyr::pull(neg_ctrl_med)

  # We expect negative controls to be now equal to 0
  testthat::expect_equal(neg_controls[2:4], c(0,0,0))

  pos_controls <- gimap_dataset$normalized_log_fc %>%
    dplyr::filter(target_type %in% c("gene_ctrl", "ctrl_gene")) %>%
    dplyr::group_by(rep) %>%
    dplyr::summarize(pos_ctrl_med = median(lfc_adj)) %>%
    dplyr::pull(pos_ctrl_med)

  # We expect positive controls to be now equal to -1
  testthat::expect_equal(round(sum(pos_controls[-1])), round(sum(-0.3979652 -0.3998103 -0.3862449)))

})
