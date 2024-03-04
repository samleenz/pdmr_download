library(tidyverse)
library(foreach)
library(doParallel)

####
# Download rnaseq data (RSEM gene level) for all PDMR samples of SOC
# 
####

# variables ---------------------------------------------------------------
outdir <- "data/raw/tmp"
pdmr_table <- read_csv("data/raw/pdmr_soc_rnaseq.csv")
cl <- 4

registerDoParallel(4)

# functions ---------------------------------------------------------------
smoosh <- function(x){
  #' paste idiom but make it a function
  #' combine the first three fields of the pdmr URL
  #' [patient]_[specimin]_[sample]
  #'
  paste0(x[1, 1:3], collapse = "_") |>
    str_remove_all("-")
}

# main --------------------------------------------------------------------
# have_primary <- pdmr_table |>
#   dplyr::filter(`PDM Type` == "Patient/Originator Specimen") |>
#   dplyr::pull(`Patient ID`)
# 
# have_pdx <- pdmr_table |>
#   dplyr::filter(`PDM Type` == "PDX") |>
#   dplyr::pull(`Patient ID`)
# 
# 
# pdmr_paired_primary <- pdmr_table |>
#   dplyr::filter(`Patient ID` %in% intersect(have_primary, have_pdx))
# 
# write_csv(pdmr_paired_primary, file = "data/raw/paired_samples.csv")

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}



saved <- foreach(nme = pdmr_table$`RSEM(genes)`, .errorhandling = "remove") %dopar% {
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
difference <- setdiff(pdmr_table$`RSEM(genes)`, saved)
if(length(difference) > 0){
  write_lines(difference, file = "data/raw/soc_samples_missing.txt")
}

