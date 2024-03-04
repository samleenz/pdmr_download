# Mon Mar  4 11:40:37 2024 ------------------------------
# make summarized experiment objects for the SOC samples

library(SummarizedExperiment)
library(tidyverse)

annot <- read_csv("data/clean/pdmr_soc_annotations.csv") |>
  as.data.frame() |>
  column_to_rownames("sample")

samples <- paste0("data/raw/tmp/", rownames(annot), ".RSEM.genes.txt")
samples <- samples[samples %in% list.files("data/raw/tmp", pattern = "genes\\.txt", full.names = TRUE)]

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


pdmr_count <- sapply(samples, read_rsem, value = "expected_count", simplify = F)
pdmr_count <- do.call(rbind, pdmr_count)
pdmr_count <- pdmr_count[,colSums(pdmr_count)>0]

pdmr_tpm <- sapply(samples, read_rsem, simplify = F)
pdmr_tpm <- do.call(rbind, pdmr_tpm)
pdmr_tpm <- pdmr_tpm[, colnames(pdmr_count)]

# convert to ensg
genes <- gprofiler2::gconvert(colnames(pdmr_tpm), target = "ENSG")


rData <- genes |>
  dplyr::filter(! target == "nan") |> 
  dplyr::distinct(input, .keep_all = T) |>
  dplyr::distinct(target, .keep_all = T) |>
  select(HGNC = input, ENSGID = target, description) |>
  mutate(rname = ENSGID) |>
  as.data.frame() |>
  column_to_rownames("rname")

lost_genes <- setdiff(colnames(pdmr_tpm), rData$HGNC)
write_lines(lost_genes, "data/clean/pdmr_soc_hgnc_lost_in_mapping.txt")
    
pdmr_tpm <- pdmr_tpm[,rData$HGNC]
colnames(pdmr_tpm) <- rData$ENSGID

pdmr_count <- pdmr_count[,rData$HGNC]
colnames(pdmr_count) <- rData$ENSGID


pdmr_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    counts = t(pdmr_count),
    TPM = t(pdmr_tpm)
  ),
  colData = annot[rownames(pdmr_tpm), ],
  rowData = rData,
  metadata = list(
    date_compiled = Sys.Date(),
    creator = "Sam Lee",
    description = "PDMR samples annotated as SOC/HGSOC."
  )
)

pdmr_se

write_rds(pdmr_se, "data/clean/pdmr_soc_se.rds")
