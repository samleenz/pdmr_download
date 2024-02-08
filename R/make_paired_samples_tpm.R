# 20220412
# sam lee

library(tidyverse)
annot_clean <- read_csv("data/clean/pdmr_paired_annotations.csv")

samples <- paste0("data/raw/", annot_clean$sample, ".RSEM.genes.txt")
samples <- samples[samples %in% list.files("data/raw", pattern = "genes\\.txt", full.names = TRUE)]

# samples <- list.files("data/raw", pattern = "*genes.txt", full.names = TRUE)
# samples <- samples[1:10]



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

pdmrPaired <- list(
  annot = column_to_rownames(annot_clean, var = "sample")[rownames(pdmr_tpm), ],
  TPM = pdmr_tpm
)

save(pdmrPaired, file = "/stornext/HPCScratch/lee.sa/data/clean/pdmrPaired.RData")


# 
# pdmr_tpm[1:5,1:5]
# 
# 
# pc <- pdmr_tpm2 |>
#   ncla::transformData("log2") |>
#   ncla::transformData("pca10")
# 
# um <- pdmr_tpm2 |>
#   ncla::transformData("log2") |>
#   umap::umap()
# 
# ncla::plotGG_UMAP_eda(um, "one")
# 
# plot(pc[,1], pc[,2])