test_that("Annotation options", {


   gimap_dataset <- get_example_data("gimap") %>%
     gimap_filter() %>%
     gimap_annotate(cell_line = "HELA")


   testthat::expect_warning(
     gimap_dataset <- get_example_data("gimap") %>%
       gimap_filter() %>%
       gimap_annotate(cell_line_annotate = FALSE) %>%
       gimap_normalize(timepoints = "day"))

   gimap_dataset <- gimap_dataset %>%
     gimap_filter() %>%
     gimap_annotate(cell_line_annotate = FALSE) %>%
     gimap_normalize(timepoints = "day",
                     normalize_by_unexpressed = FALSE)


})
