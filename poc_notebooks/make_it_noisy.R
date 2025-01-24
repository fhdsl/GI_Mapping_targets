
set.seed(1234)
library(magrittr)
example_data <- readr::read_tsv("counts_pgPEN.PC.correct.v3.txt")


noised_up_data <-
  example_data %>%
  dplyr::mutate_at(dplyr::vars(pretreatment:osimB), ~.x + runif(nrow(example_data), min =-2, max = 5)) %>%
  dplyr::mutate_at(dplyr::vars(pretreatment:osimB), ~.x + rnorm(nrow(example_data), mean = 0, sd = 2)) %>%
  dplyr::mutate_at(dplyr::vars(pretreatment:osimB), round) %>%
  dplyr::mutate_at(dplyr::vars(pretreatment:osimB), abs) %>%
  dplyr::rename(id = pgRNA, drug1A = osimA, drug1B = osimB) %>%
  readr::write_tsv("counts_pgPEN_PC9_example.tsv")
