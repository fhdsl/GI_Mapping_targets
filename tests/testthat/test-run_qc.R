test_that("HTML file is created and content is correct", {
  testthat::skip_on_cran()
  gimap_dataset <- get_example_data("gimap")

  html_file <- run_qc(gimap_dataset, open_results = FALSE, overwrite = TRUE)

  expect_true(file.exists(html_file))

  qc_cdf(gimap_dataset)

  qc_sample_hist(gimap_dataset)

  qc_variance_hist(gimap_dataset)

  qc_constructs_countzero_bar(gimap_dataset)

  qc_cor_heatmap(gimap_dataset)

  qc_plasmid_histogram(gimap_dataset)

  gimap_dataset <- get_example_data("gimap_treatment")

  html_file <- run_qc(gimap_dataset, open_results = FALSE, overwrite = TRUE)

  expect_true(file.exists(html_file))
})
