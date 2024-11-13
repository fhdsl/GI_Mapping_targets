test_that("Test annotation step", {

  # Test with annotations
  gimap_dataset <- get_example_data("gimap") %>%
    gimap_filter() %>%
    gimap_annotate(cell_line = "HELA")

  testthat::expect_named(gimap_dataset$annotation,
                         c("pgRNA_id", "gRNA1_seq", "gRNA2_seq", "target_type",
                           "paralog_pair", "paralog_pair_id", "gene1_ensembl_id", "gene1_symbol",
                           "gene2_symbol", "gene2_ensembl_id", "log2_cn_gene1", "log2_cn_gene2",
                           "gene1_essential_flag", "gene2_essential_flag",  "pgRNA_target", "log2_tpm_gene1",
                           "gene1_expressed_flag",  "log2_tpm_gene2", "gene2_expressed_flag", "norm_ctrl_flag",
                           "unexpressed_ctrl_flag")
  )


  # Test without annotations
  gimap_dataset <- get_example_data("gimap") %>%
    gimap_filter() %>%
    gimap_annotate(depmap_annotate = FALSE)

})
