# This script was supposed to run on phoenix, but the lack of internet access means I need to run it locally.
# The generated file was moved to the /resources folder on phoenix

library(magrittr)
library(stringr)
library(tibble)
library(AnnotationHub)

  ## Setup to allow user to run interactively or as part of snakemake workflow
  # if (exists("snakemake")) {
  #   ens_species <- snakemake@params[["species"]] %>%
  #     str_replace("_", " ")
  #   ens_release <- snakemake@params[["ensembl_release"]] %>%
  #     as.character()
  #   intervalFile <- "resources/exons.intervals"
  # } else {
  #   ## Edit these variables manually if running interactively
  #   ens_species <- "Mus musculus"
  #   ens_release <- "101"
  #   intervalFile <- dirname(rstudioapi::getSourceEditorContext()$path) %>%
  #     str_remove("workflow/scripts") %>%
  #     paste0("resources/exons.intervals")
  # }
  #
  # ## Restrict sequences to primary assembly
  # if (ens_species == "Mus musculus") {
  #   primary_chrs <- c(seq(1:19), "X", "Y")
  # } else if (ens_species == "Danio rerio") {
  #   primary_chrs <- c(seq(1:25))
  # } else if (ens_species == "Homo sapiens") {
  #   primary_chrs <- c(seq(1:22), "X", "Y")
  # }

primary_chrs <- c(seq(1:25))

## Get exon intervals via AnnotationHub
ah <- AnnotationHub() %>%
  subset(species == "Danio rerio") %>%
  subset(rdataclass == "EnsDb")


ahId <- ah$ah_id[str_detect(ah$title, "101")]
ensDb <- ah[[ahId]]
exons <- exonsBy(ensDb, by = "gene")
exonRanges <- exons %>%
    unlist() %>%
    GenomicRanges::reduce() %>%
    as.data.frame() %>%
    dplyr::filter(seqnames %in% primary_chrs) %>%
    dplyr::select(chromosome = seqnames, start, end)

  ## Save intervals in correct file format
  if (!dir.exists(dirname(intervalFile))) {
    dir.create(dirname(intervalFile), recursive = TRUE)
  }
  tibble(
    interval = paste0(
      exonRanges$chromosome,
      ":",
      exonRanges$start,
      "-",
      exonRanges$end
    )
  ) %>%
    write.table(
      file = "exons.intervals",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
