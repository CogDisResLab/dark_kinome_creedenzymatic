# Creedenzymatic Analysis

library(tidyverse)
library(creedenzymatic)

process_creedenzymatic <-
  function(krsa_path, uka_path, peptide_path) {
    krsa_data <- read_tsv(krsa_path, show_col_types = FALSE) |>
      select(Kinase, Score = AvgZ) |>
      read_krsa(trns = "abs", sort = "desc")

    uka_data <- read_tsv(uka_path, show_col_types = FALSE) |>
      select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
      read_uka(trns = "abs", sort = "desc")

    peptide_data <-
      read_tsv(peptide_path, show_col_types = FALSE) |>
      select(Peptide, Score = totalMeanLFC)

    kea3_data <-
      read_kea(
        peptide_data,
        sort = "asc",
        trns = "abs",
        method = "MeanRank",
        lib = "kinase-substrate"
      )

    # ptmsea_data <-
    #   read_ptmsea(peptide_data)

    combined <- combine_tools(
      KRSA_df = krsa_data,
      UKA_df = uka_data,
      KEA3_df = kea3_data
      # PTM_SEA_df = ptmsea_data
    )

    combined
  }

krsa_files <- c(
  "results/acrossChip_KRSA_Table_females.txt",
  "results/acrossChip_KRSA_Table_males.txt"
)

uka_files <- c(
  "results/UKA_Femalestxt.txt",
  "results/UKA_Males.txt"
)

peptide_files <- c(
  "results/females_LFC_df.txt",
  "results/males_LFC_df.txt"
)

result <-
  list(
    krsa_path = krsa_files,
    uka_path = uka_files,
    peptide_path = peptide_files
  ) |>
  pmap(process_creedenzymatic) |>
  set_names(c("Females", "Males")) |>
  imap_dfr(~ write_csv(.x, str_glue("results/{.y}_creedenzymatic.csv")), .id = "Comparison")
