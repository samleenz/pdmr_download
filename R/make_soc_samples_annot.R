# Mon Mar  4 11:05:14 2024 ------------------------------
# rnaseq table downloaded from PDMR website aftering querying genomic
# analysis:rnaseq for all samples with the OncoTree code 'SOC'
# sample table from the sample tab
library(tidyverse)

make_sample_name <- function(x){
  x |>
    mutate(`Specimen ID` = str_remove(`Specimen ID`, "-")) |>
    mutate(`Sample ID` = str_remove(`Sample ID`, "-")) |>
    mutate(sample = paste(`Patient ID`, `Specimen ID`, `Sample ID`, sep = "_"))
}

soc_samples <- read_csv("data/raw/pdmr_soc_samples.csv") |>
  make_sample_name()
dat <- read_csv("data/raw/pdmr_soc_rnaseq.csv") |>
  make_sample_name()

dat 

annot <- dat |>
  left_join(
    soc_samples |>
      select(- `Patient ID`, - `Specimen ID`, - `Sample ID`, - `CTEP SDCDescription`, - `DiagnosisSubtype`, - `OncoTreeCode`, - `Disease BodyLocation`)
    ) |>
  select(
    sample, 
    `Patient ID`,
    `Specimen ID`,
    `Sample ID`,
    `PDM Type`,
    Passage, 
    DiagnosisSubtype,
    OncoTreeCode,
    PathologyAvail,
    `OncoKB Cancer GenePanel Data Avail`,
    `Whole ExomeSequence Avail`,
    RNASeqAvail
  ) 

write_csv(annot, file = "data/clean/pdmr_soc_annotations.csv")

