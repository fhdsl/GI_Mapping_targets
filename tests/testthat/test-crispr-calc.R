test_that("Test CRISPR calculations", {
   gimap_dataset <- get_example_data("gimap") %>%
     gimap_filter() %>%
     gimap_annotate(cell_line = "HELA") %>%
     gimap_normalize(
       timepoints = "day"
     ) %>%
     calc_crispr()

  # Test that single crisprs are turning out the same
  testthat::expect_true(round(gimap_dataset$single_crispr_score$single_crispr[1], 3) == round(4.762139, 3))
  testthat::expect_true(round(gimap_dataset$single_crispr_score$mean_single_crispr[1], 3) == round(3.141807, 3))
  testthat::expect_true(round(gimap_dataset$single_crispr_score$expected_single_crispr[1], 3) == round(3.976161, 3))

  # Test that double crisprs are turning out the same
  testthat::expect_true(round(gimap_dataset$double_crispr_score$expected_double_crispr[1], 3)  == round(3.376975, 3))
  testthat::expect_true(round(gimap_dataset$double_crispr_score$mean_single_crispr_1[1], 3) == round(0.2351685, 3))
  testthat::expect_true(round(gimap_dataset$double_crispr_score$double_crispr[1], 3) == round(-8.034962, 3))

})
