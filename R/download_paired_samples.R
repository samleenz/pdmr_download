library(tidyverse)

pdmr_table <- read_csv("~/data/clean/pdmr_sample_table_rnaseq.csv")


have_primary <- pdmr_table |>
  filter(`PDM Type` == "Patient/Originator Specimen") |>
  pull(`Patient ID`)


pdmr_paired_primary <- pdmr_table |>
  filter(`Patient ID` %in% have_primary)


outdir <- "~/Desktop/20220411_test"


smoosh <- function(x){
  #' paste idiom but make it a function
  #' combine the first three fields of the pdmr URL
  #' [patient]_[specimin]_[sample]
  #'
  paste0(x[1, 1:3], collapse = "_")
}

for(nme in pdmr_paired_primary$`RSEM(genes)`[1:10]){
  
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
  
  
  # print(clean_nme)
  
}

x <- read_tsv(file.path(test, paste0(clean_nme, ".RSEM.genes.txt")))


x
