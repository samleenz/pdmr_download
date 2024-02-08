raw <- read_csv("data/raw/pdm_samples.csv")
ctep <- read_csv("data/raw/SDCv10_M10_ctep.txt", skip = 1)
pdmr_paired_primary <- read_csv("data/raw/paired_samples.csv")

annot_paired <- dplyr::filter(raw, `Sample ID` %in% pdmr_paired_primary$`Sample ID`) |>
  mutate(sample = paste(annot_paired$`Patient ID`, annot_paired$`Specimen ID`, annot_paired$`Sample ID`, sep = "_")) |>
  select(
    sample, 
    `Patient ID`,
    `Specimen ID`,
    `Sample ID`,
    `PDM Type`,
    Passage, 
    `Disease BodyLocation`,
    `CTEP SDCDescription`,
    `CTEP SDCCode`
    )

filter(annot_paired, sample %in% rownames(pdmr_tpm)) |>
  View()

write_csv(annot_paired, file = "data/clean/pdmr_paired_annotations.csv")
