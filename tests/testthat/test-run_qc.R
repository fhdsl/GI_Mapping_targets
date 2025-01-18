test_that("HTML file is created and content is correct", {
  gimap_dataset <- get_example_data("gimap")

  html_file <- run_qc(gimap_dataset, open_results = FALSE)

  expect_true(file.exists(html_file))
})
