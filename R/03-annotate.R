#' Annotate gimap data
#' @description In this function, a `gimap_dataset` is annotated as far as which
#' genes should be used as controls.
#' @param .data Data can be piped in with tidyverse pipes from function to function.
#'  But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param annotation_file If no file is given, will attempt to use the design file
#' from
#' https://media.addgene.org/cms/filer_public/a9/9a/a99a9328-324b-42ff-8ccc-30c544b899e4/pgrna_library.xlsx
#' @param control_genes A vector of gene symbols (e.g. AAMP) that should be labeled
#' as control genes. These will be used for log fold change calculations. If no list is given
#' then DepMap Public 23Q4 Achilles_common_essentials.csv is used https://depmap.org/portal/download/all/
#' @param cell_line_annotate (Optional) TRUE or FALSE you'd also like to have cell_line_annotation
#' from DepMap.
#' @param cell_line which cell line are you using? (e.g., HELA, PC9, etc.).
#' Required argument if cell_line_annotate is TRUE.
#' @param custom_tpm (Optional) You may supply your own data frame of transcript
#' per million expression to be used for this calculation if you can't or don't
#' want to use DepMap data annotation for your cell_line. This data frame needs to have
#' two columns: 'log2_tpm' that has the log2 tpm expression data for this cell line and
#' and 'genes' which needs to be gene symbols that match those in the data. eg. "NDL1".
#' Note that you can use custom_tpm with cell_line_annotate but your custom_tpm will be used
#' instead of the tpm data from DepMap. However other data from DepMap like CN will be added.
#' @importFrom stringr word
#' @import dplyr
#' @export
#' @examples \dontrun{
#'
#' gimap_dataset <- get_example_data("gimap")
#'
#' # Highly recommended but not required
#' run_qc(gimap_dataset)
#'
#' # By default DepMap annotation will be used to determine genes which are
#' # unexpressed. In the `gimap_normalize` this will by default be used to normalize to.
#' gimap_dataset <- gimap_dataset %>%
#'   gimap_filter() %>%
#'   gimap_annotate(cell_line = "HELA")
#'
#' # To see anotations you can do this:
#' gimap_dataset$annotation
#'
#'
#' # You can also say cell_line_annotate = false if you don't want to use DepMap
#' # annotation BUT if you don't also specify that you say you are
#' # `normalize_by_unexpressed = FALSE` in the normalize step you will get a warning.
#' gimap_dataset <- get_example_data("gimap") %>%
#'  gimap_filter() %>%
#'  gimap_annotate(cell_line_annotate = FALSE) %>%
#'  gimap_normalize(timepoints = "day",
#'                  normalize_by_unexpressed = FALSE)
#'
#' ### CUSTOM TPM example
#' # Lastly, this is also an option:
#' where custom data is provided to `custom_tpm` is a data frame with
#' `genes` and `log2_tpm` as the columns.
#' gimap_dataset <- get_example_data("gimap") %>%
#'  gimap_filter() %>%
#'  gimap_annotate(custom_tpm = custom_tpm)) %>%
#'  gimap_normalize(timepoints = "day")
#'
#' }
gimap_annotate <- function(.data = NULL,
                           gimap_dataset,
                           annotation_file = NULL,
                           control_genes = NULL,
                           cell_line_annotate = TRUE,
                           custom_tpm = NULL,
                           cell_line = NULL) {
  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) {
    stop("This function only works with gimap_dataset objects which can be made",
    "with the setup_data() function.")
  }

  if (cell_line_annotate & is.null(cell_line)) {
    stop("By default `cell_line_annotate` = TRUE which means if you want DepMap annotation",
         " you must specify which cell line you are using with the `cell_line` argument",
         " OR set `cell_line_annotate` = FALSE.",
         " However this means you cannot use expression data for normalization",
         " Aka in gimap_normalize(): `normalize_by_unexpressed` must be set to FALSE",
         " unless you supply your own expression data using the `custom_tpm` argument")
  }

  # Get the annotation data based on the pg construct design
  if (!is.null(annotation_file)) {
    if (!file.exists(annotation_file)) {
      stop("The annotation_file specified cannot be found. Please double check the file path")
    }
    annotation_df <- read_table(annotation_file)
  } else {
    annotation_df <- get_example_data("annotation")
  }

  message("Annotating Data")

  ############################ CONTROL GENE ANNOTATION #########################
  # If control genes aren't provided then we get some from DepMap
  if (!is.null(control_genes)) {
    if (!file.exists(control_genes)) {
      stop("The annotation_file specified cannot be found. Please double check the file path")
    }
    control_genes <- read_table(control_genes)[, 1]
  } else {
    # This file is from https://depmap.org/portal/download/all/ and from DepMap Public 19Q3 All Files
    # Essential gene labeling is from inst/extdata/Achilles_common_essentials.csv
    control_genes <- crtl_genes()
  }

  ############################ Get TPM data ####################################
  if (cell_line_annotate) {
    # This is used to flag things
    ## get TPM and CN information (w/ option for user to upload their own info)
    depmap_metadata <- readr::read_csv(
      "https://figshare.com/ndownloader/files/35020903", show_col_types = FALSE)

    my_depmap_id <- depmap_metadata %>%
      dplyr::filter(stripped_cell_line_name == toupper(cell_line)) %>%
      dplyr::pull(DepMap_ID)

    if (length(my_depmap_id) == 0) {
      stop("The cell line specified, ",
           cell_line,
           "was not found in the DepMap data. Run supported_cell_lines() to see the full list")
    }

    tpm_file <- file.path(
      system.file("extdata", package = "gimap"),
      "CCLE_expression.csv")

    if (!file.exists(tpm_file)) tpm_setup()

    tpm <- readr::read_csv(tpm_file,
      show_col_types = FALSE,
      col_select = c("genes", dplyr::all_of(my_depmap_id))
    ) %>%
      dplyr::rename(log2_tpm = my_depmap_id)

    ############################ COPY NUMBER ANNOTATION ##########################
    cn_file <- file.path(system.file("extdata", package = "gimap"), "CCLE_gene_cn.csv")
    if (!file.exists(cn_file)) cn_setup()

    # Read in the CN data
    depmap_cn <- readr::read_csv(cn_file,
      show_col_types = FALSE,
      col_select = c("genes", dplyr::all_of(my_depmap_id))
    ) %>%
      dplyr::rename(log2_cn = dplyr::all_of(my_depmap_id))

    annotation_df <- annotation_df %>%
      dplyr::left_join(depmap_cn, by = c("gene1_symbol" = "genes")) %>%
      dplyr::left_join(depmap_cn, by = c("gene2_symbol" = "genes"),
                       suffix = c("_gene1", "_gene2"))
  }

  # If people supply their own tpm file we need to check it for stuff
  if (!is.null(custom_tpm)) {
    stopifnot("custom_tpm must be a data.frame or tibble" =
                is.data.frame(custom_tpm),
              "custom_tpm must contain a column called 'log2_tpm'" = "log2_tpm" %in%
                colnames(custom_tpm),
              "custom_tpm must contain a column called 'genes'" = "genes" %in%
                colnames(custom_tpm))

    gene_matches <- sum(match(custom_tpm$genes, annotation_df$genes))
    percent <- gene_matches/length(custom_tpm$genes)*100

    message(gene_matches,
            " :number of genes have matches in the custom_tpm data \n",
            percent, "% :percent of genes with matches in the custom_tpm data")

    stopifnot("less than half of the genes have custom_tpm matches" =
                percent > 50)
  }
  if (!is.null(custom_tpm) | cell_line_annotate) {
  tpm <- tpm %>%
    dplyr::mutate(expressed_flag = dplyr::case_when(
      log2_tpm < 1 ~ FALSE,
      log2_tpm >= 1 ~ TRUE,
      is.na(log2_tpm) ~ NA
    ))
  }
  ############################ ANNOTATION COMBINING ############################
  # This set up is more or less the same as the original
  annotation_df <- annotation_df %>%
    dplyr::mutate(
      gene1_essential_flag = gene1_symbol %in% control_genes,
      gene2_essential_flag = gene2_symbol %in% control_genes,
      pgRNA_target = dplyr::case_when(
        target_type == "gene_gene" ~ paste(gene1_symbol, gene2_symbol, sep = "_"),
        target_type == "gene_ctrl" ~ paste(gene1_symbol, "ctrl", sep = "_"),
        target_type == "ctrl_gene" ~ paste("ctrl", gene2_symbol, sep = "_"),
        TRUE ~ target_type
      )
    ) %>%
    dplyr::mutate(norm_ctrl_flag = dplyr::case_when(
      target_type == "gene_gene" ~ "double_targeting",
      target_type == "gene_ctrl" & gene1_essential_flag == TRUE ~ "positive_control",
      target_type == "ctrl_gene" & gene2_essential_flag == TRUE ~ "positive_control",
      target_type == "gene_ctrl" & gene1_essential_flag != TRUE ~ "single_targeting",
      target_type == "ctrl_gene" & gene2_essential_flag != TRUE ~ "single_targeting",
      target_type == "ctrl_ctrl" ~ "negative_control"
    )) %>%
    dplyr::mutate(
      norm_ctrl_flag = factor(norm_ctrl_flag, levels = c(
        "negative_control",
        "positive_control",
        "single_targeting",
        "double_targeting"
      ))
    )
  if (cell_line_annotate | !is.null(custom_tpm)) {
    annotation_df <- annotation_df %>%
      dplyr::left_join(tpm, by = c("gene1_symbol" = "genes")) %>%
      dplyr::rename(gene1_expressed_flag = expressed_flag) %>%
      dplyr::left_join(tpm, by = c("gene2_symbol" = "genes"), suffix = c("_gene1", "_gene2")) %>%
      dplyr::rename(gene2_expressed_flag = expressed_flag) %>%
      dplyr::mutate(
        unexpressed_ctrl_flag = dplyr::case_when(
          norm_ctrl_flag == "double_targeting" & gene1_expressed_flag == FALSE &
            gene2_expressed_flag == FALSE ~ TRUE,
          norm_ctrl_flag == "single_targeting" & (gene1_expressed_flag == FALSE |
                                                    gene2_expressed_flag == FALSE) ~ TRUE,
          TRUE ~ FALSE
        )
      )
  }
  ################################ STORE IT ####################################

  if (gimap_dataset$filtered_data$filter_step_run) {
    keep_for_annotdf <- annotation_df$pgRNA_id %in%
      unlist(gimap_dataset$filtered_data$metadata_pg_ids)
    annotation_df <- annotation_df[keep_for_annotdf, ]
  }

  gimap_dataset$annotation <- annotation_df

  return(gimap_dataset)
}


# This function sets up the tpm data from DepMap is called by the `gimap_annotate()` function
tpm_setup <- function() {
  tpm_file <- file.path(
    system.file("extdata", package = "gimap"),
    "CCLE_expression.csv"
  )

  download.file("https://figshare.com/ndownloader/files/34989919",
    destfile = tpm_file)

  data_df <- readr::read_csv(tpm_file,
    show_col_types = FALSE,
    name_repair = make.names
  )

  cell_line_ids <- data_df$X

  genes <- stringr::word(colnames(data_df)[-1], sep = "\\.\\.", 1)

  colnames(data_df) <- c("cell_line_ids", genes)

  data_df <- as.data.frame(t(data_df[, -1]))
  colnames(data_df) <- cell_line_ids
  data_df$genes <- genes

  data_df %>%
    dplyr::select(genes, dplyr::everything()) %>%
    readr::write_csv(tpm_file)

  return(tpm_file)
}

# This function sets up the tpm data from DepMap is called by the `gimap_annotate()`
# function if the cn_annotate = TRUE
cn_setup <- function() {
  options(timeout = 1000)

  cn_file <- file.path(
    system.file("extdata", package = "gimap"),
    "CCLE_gene_cn.csv"
  )

  download.file("https://figshare.com/ndownloader/files/34989937",
    destfile = cn_file)

  data_df <- readr::read_csv(cn_file,
    show_col_types = FALSE,
    name_repair = make.names
  )

  cell_line_ids <- data_df$X

  genes <- stringr::word(colnames(data_df)[-1], sep = "\\.\\.", 1)

  colnames(data_df) <- c("cell_line_ids", genes)

  data_df <- as.data.frame(t(data_df[, -1]))
  colnames(data_df) <- cell_line_ids
  data_df$genes <- genes

  data_df %>%
    dplyr::select(genes, dplyr::everything()) %>%
    readr::write_csv(cn_file)

  return(cn_file)
}

# This function sets up the control genes file from DepMap is called by the `gimap_annotate()`
crtl_genes <- function() {
  crtl_genes_file <- file.path(
    system.file("extdata", package = "gimap"),
    "Achilles_common_essentials.csv"
  )

  if (!file.exists(crtl_genes_file)) {
    download.file("https://figshare.com/ndownloader/files/34989871",
      destfile = crtl_genes_file)

    crtl_genes <- readr::read_csv(crtl_genes_file, show_col_types = FALSE) %>%
      tidyr::separate(col = gene, into = c("gene_symbol", "entrez_id"), remove = FALSE, extra = "drop")

    readr::write_csv(crtl_genes, crtl_genes_file)
  } else {
    crtl_genes <- readr::read_csv(crtl_genes_file)
  }

  return(crtl_genes$gene_symbol)
}


#' List the supported cell lines
#' @description This function downloads the metadata for DepMap and lists which cell
#' lines are supported.
#' @export
#' @examples \dontrun{
#'
#' cell_lines <- supported_cell_lines()
#' }
supported_cell_lines <- function() {
  depmap_metadata <- readr::read_csv("https://figshare.com/ndownloader/files/35020903",
                                     show_col_types = FALSE)

  return(sort(depmap_metadata$stripped_cell_line_name))
}
