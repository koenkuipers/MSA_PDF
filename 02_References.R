# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# libraries
library(dplyr)
library(stringr)

# Data
## All external data used is publicly available:
### PREDICTS: https://doi.org/10.5519/j4sh7e0w
### BioTIME: https://biotime.st-andrews.ac.uk/download_form.php?dl=csv_download 

## We have downloaded the data and stored the csv datasets in our 'Data' folder with the following (original) names:
### PREDICTS: PREDICTS_references.csv
### BioTIME metadata: BioTIMECitations_24_06_2021.csv

pred.df <- read.csv("./PREDICTS_references.csv")
kuip.df <- read.csv("./Data/Agridiv_data.csv")
mido.df <- read.csv("./Data/Midolo_etal_2018_data_N-application_FINAL_IA.csv")
biot.df <- read.csv("./Data/BioTIMECitations_24_06_2021.csv")
dat.df <- read.csv("./Data/01_PDF_MSA.csv")

# modify data
## PREDICTS
pred.df <- rename(pred.df, Authors_year = Reference)
pred.df <- mutate(pred.df, Database = "PREDICTS",
                  Source_ID = as.character(Source_ID),
                  Year = as.integer(Year))
pred.df <- select(pred.df, Database, Source_ID, Authors_year, Year, Journal, DOI)

## Kuipers
kuip.df <- rename(kuip.df, Authors_year = Reference, Year = Reference_yr)
kuip.df <- mutate(kuip.df,
                  Authors_year = str_replace_all(Authors_year, "_", " "),
                  Database = "Kuipers_GCB_2023")
kuip.df <- mutate(kuip.df, Journal = word(Authors_year, str_count(Authors_year, " ") + 1, sep = " "))
kuip.df <- mutate(kuip.df, Authors_year = word(Authors_year, 1, str_count(Authors_year, " "), sep = " "))
kuip.df <- mutate(kuip.df, Authors_year = str_replace_all(Authors_year, "et al ", "et al. "),
                  Source_ID = as.character(Source_ID),
                  Year = as.integer(Year))
kuip.df <- select(kuip.df, Database, Source_ID, Authors_year, Year, Journal, DOI)
kuip.df <- distinct(kuip.df)

## Midolo
mido.df <- mutate(mido.df,
                  Database = "Midolog_GEB_2019",
                  Source_ID = Experiment,
                  Authors_year = paste(word(Experiment, 2, str_count(Experiment, " ") + 1, sep = " "), word(Experiment, 1, sep = " "), sep = " "),
                  Year = word(Experiment, 1, sep = " "),
                  Journal = NA,
                  DOI = NA)
mido.df <- mutate(mido.df, Authors_year = str_replace_all(Authors_year, " A", " "))
mido.df <- mutate(mido.df, Authors_year = str_replace_all(Authors_year, " B", " "))
mido.df <- mutate(mido.df, Authors_year = str_replace_all(Authors_year, " C", " "),
                  Source_ID = as.character(Source_ID),
                  Year = as.integer(Year))
mido.df <- select(mido.df, Database, Source_ID, Authors_year, Year, Journal, DOI)
mido.df <- distinct(mido.df)

## BioTIME
biot.df <- mutate(biot.df,
                  Database = "BioTIME",
                  Source_ID = as.character(STUDY_ID),
                  Year = word(word(CITATION_LINE, 2, sep = "\\("), 1, sep = "\\)"),
                  Journal = word(word(word(CITATION_LINE, 2, sep = "\\) "), 2, sep = "\\. "), 1, sep = ", "),
                  DOI = NA)
biot.df <- mutate(biot.df, Authors_year = paste(word(CITATION_LINE, 1, sep = " \\("), Year, sep = " "))
biot.df$DOI[str_which(biot.df$Journal, "doi")] <- biot.df$Journal[str_which(biot.df$Journal, "doi")]
biot.df <- select(biot.df, Database, Source_ID, Authors_year, Journal, DOI)
biot.df <- distinct(biot.df)

## Reference dataframe
ref.df <- bind_rows(pred.df, kuip.df, mido.df, biot.df)
ref.df <- rename(ref.df, Source_ID_original = Source_ID)

rm(pred.df, kuip.df, mido.df, biot.df)

ref.df <- left_join(ref.df, distinct(select(dat.df, Database, Source_ID, Source_ID_original)))
ref.df <- filter(ref.df, !is.na(Source_ID))
ref.df <- select(ref.df,
                 Database, Source_ID, Source_ID_original,
                 Authors_year, Year, Journal, DOI)
ref.df <- arrange(ref.df, Source_ID)

# save data
write.csv(ref.df, "./Data/02_References.csv", row.names = FALSE)
