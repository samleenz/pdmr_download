library(tidyverse)
library(foreach)
library(doParallel)

####
# Download rnaseq data (RSEM gene level) for all PDMR samples that have data available for the primary tumour
# 
####

# variables ---------------------------------------------------------------
outdir <- "data/raw"
pdmr_table <- read_csv("data/clean/pdmr_sample_table_rnaseq.csv")
cl <- 4

registerDoParallel(4)

# functions ---------------------------------------------------------------
smoosh <- function(x){
  #' paste idiom but make it a function
  #' combine the first three fields of the pdmr URL
  #' [patient]_[specimin]_[sample]
  #'
  paste0(x[1, 1:3], collapse = "_")
}

# main --------------------------------------------------------------------
have_primary <- pdmr_table |>
  dplyr::filter(`PDM Type` == "Patient/Originator Specimen") |>
  dplyr::pull(`Patient ID`)

have_pdx <- pdmr_table |>
  dplyr::filter(`PDM Type` == "PDX") |>
  dplyr::pull(`Patient ID`)


pdmr_paired_primary <- pdmr_table |>
  dplyr::filter(`Patient ID` %in% intersect(have_primary, have_pdx))

write_csv(pdmr_paired_primary, file = "data/raw/paired_samples.csv")

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
  }



saved <- foreach(nme = pdmr_paired_primary$`RSEM(genes)`, .errorhandling = "remove") %dopar% {
  ## create a nice name for saving
  clean_nme <- nme |>
    basename() |>
    stringr::str_split(pattern = "~", simplify = T) |>
    smoosh()
  
  ## save the file
  download.file(
    url = paste0("https://pdmdb.cancer.gov", nme),
    destfile = file.path(outdir, paste0(clean_nme, ".RSEM.genes.txt"))
  )
  
  nme
}

## see if any files were missed
difference <- setdiff(pdmr_paired_primary$`RSEM(genes)`, saved)
if(length(difference) > 0){
  write_lines(difference, file = "data/raw/paired_samples_missing.txt")
}
  
