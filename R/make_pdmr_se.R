# make summarized experiment objects

library(SummarizedExperiment)
library(tidyverse)

annot <- read_csv("data/raw/pdmr_sample_info_has_rnaseq_20240208.csv")
annot <- mutate(annot, sample = paste(`Patient ID`, `Specimen ID`, `Sample ID`, sep = "_")) |>
  filter(`PDM Type` == "PDX") |>
  column_to_rownames("sample")

samples <- paste0("data/raw/", rownames(annot), ".RSEM.genes.txt")
samples <- samples[samples %in% list.files("data/raw", pattern = "genes\\.txt", full.names = TRUE)]

read_rsem <- function(f, value = "TPM"){
  x <- readr::read_tsv(f, show_col_types = FALSE)
  x <- x[,c("gene_id", value)]
  colnames(x) <- c(
    "gene_id", 
    stringr::str_remove(basename(f), "\\.RSEM\\.genes\\.txt")
  )
  
  x |>
    tibble::column_to_rownames(var = "gene_id") |>
    t()
}


pdmr_tpm <- sapply(samples, read_rsem, simplify = F)
pdmr_tpm <- do.call(rbind, pdmr_tpm)
pdmr_tpm <- pdmr_tpm[,colSums(pdmr_tpm)>0]

# convert to ensg
genes <- gprofiler2::gconvert(colnames(pdmr_tpm), target = "ENSG") |>
  dplyr::filter(! target == "nan") |> 
  dplyr::distinct(input, .keep_all = T)

pdmr_tpm <- pdmr_tpm[,genes$input]
colnames(pdmr_tpm) <- genes$target


pdmr_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    TPM = t(pdmr_tpm)
  ),
  colData = annot[rownames(pdmr_tpm), ],
  metadata = list(
    date_compiled = Sys.Date(),
    creator = "Sam Lee",
    description = "PDMR PDX models with matched primary tumours in the database."
  )
)

pdmr_se <- scuttle::logNormCounts(pdmr_se, assay.type = "TPM", name = "logTPM")

write_rds(pdmr_se, "data/clean/pdmr_paired_pdx_se.rds")
