# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# libraries
library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(igraph)

sf_use_s2(use_s2 = FALSE)

# Data
## All data used is publicly available:
### PREDICTS: https://doi.org/10.5519/j4sh7e0w
### Kuipers GCB 2023: http://doi.org/10.5281/zenodo.8199982
### Midolo GEB 2019: https://www.globio.info/impacts-of-nitrogen-addition-on-plant-species-richness-and-abundance-a-global-meta%e2%80%90analysis 
### BioTIME: https://biotime.st-andrews.ac.uk/download_form.php?dl=csv_download 

## We have downloaded the data and stored the csv datasets in our 'Data' folder with the following (original) names:
### PREDICTS: PREDICTS_database.csv
### Kuipers GCB 2023: Agridiv_data.csv
### Midolo GEB 2019: Midolo_etal_2018_data_N-application_FINAL_IA.csv
### BioTIME: BioTIMEQuery_24_06_2021.csv
### BioTIME metadata: BioTIMEMetadata_24_06_2021.csv

# Functions
link_dist <- function(dat.df, dist_m, pressure = TRUE) {
  dat.sf <- sf::st_as_sf(distinct(select(dat.df, Longitude, Latitude)),
                         coords = c("Longitude", "Latitude"),
                         crs = sf::st_crs(4326),
                         remove = FALSE)
  dat.sf <- sf::st_transform(dat.sf, crs = "ESRI:54012")
  dat.sf <- sf::st_buffer(dat.sf, dist = dist_m)
  if (pressure == TRUE) {
    dat.sf <- dplyr::left_join(dat.sf,
                             distinct(select(dat.df,
                                             SS_ID,
                                             Longitude, Latitude,
                                             Pressure_level)))
  } else {
    dat.sf <- dplyr::left_join(dat.sf,
                               distinct(select(dat.df,
                               SS_ID,
                               Longitude, Latitude)))
  }
  for (i in unique(dat.df$SS_ID)) {
    if (which(unique(dat.df$SS_ID) == i) %% 10 == 0) {
      print(paste(which(unique(dat.df$SS_ID) == i),
                  length(unique(dat.df$SS_ID)),
                  sep = "/"))
    }
    int.l <- sf::st_intersects(filter(dat.sf, SS_ID == i),
                               filter(dat.sf, SS_ID == i),
                               sparse = TRUE)
    clust.l <- igraph::graph_from_adj_list(int.l)
    dat.sf$Site_ID[dat.sf$SS_ID == i] <- igraph::components(clust.l)$membership
    rm(int.l, clust.l)
    gc()
  }
  dat.sf
}

pdf_msa <- function(dat.df) {
  dat.df <- mutate(dat.df, RIA = Mean_measurement_t / Mean_measurement_c)
  dat.df <- filter(dat.df, !is.na(RIA) | is.finite(RIA), Mean_measurement_c > 0)
  dat.df$RIA[dat.df$RIA > 1] <- 1
  if("Taxon_group" %in% colnames(dat.df)) {
    dat.df <- group_by(dat.df,
                       SS_ID,
                       Longitude,
                       Latitude,
                       Taxon_group,
                       Pressure_level)
  } else {
     dat.df <- group_by(dat.df,
                       SS_ID,
                       Longitude,
                       Latitude,
                       Pressure_level)
    }
}

# PREDICTS ----------------------------------------------------------------
dat.df <- read.csv("./Data/PREDICTS_database.csv")

# modify data
## select columns and filter rows
dat.df <- select(dat.df,
                 Source_ID, SS,
                 Longitude, Latitude,
                 Predominant_land_use, Use_intensity,
                 Rank, Taxon_number, Kingdom, Phylum, Class,
                 Diversity_metric_type, Diversity_metric,
                 Effort_corrected_measurement)
dat.df <- filter(dat.df,
                 Diversity_metric_type == "Abundance",
                 Diversity_metric %in% c("abundance", "density", "effort-corrected abundance", "sign abundance", "group abundance", "biomass", "percent cover"),
                 Rank == "Species",
                 Predominant_land_use != "Cannot decide",
                 Kingdom != "Fungi")

## rename columns
dat.df <- rename(dat.df,
                 SS_ID = SS,
                 Measurement = Effort_corrected_measurement,
                 Pressure_level = Predominant_land_use)

## rename attributes
dat.df$Pressure_level[str_which(dat.df$Pressure_level, "Secondary")] <- "Secondary vegetation"
dat.df$Pressure_level[dat.df$Pressure_level == "Primary vegetation" & dat.df$Use_intensity != "Minimal use"] <- "Primary vegetation - Light/intense use"
dat.df$Pressure_level[dat.df$Pressure_level == "Primary vegetation" & dat.df$Use_intensity == "Minimal use"] <- "Primary vegetation - Minimal use"
dat.df$Pressure_level[dat.df$Pressure_level == "Plantation forest"] <- "Forest_plantation"

dat.df$Taxon_group[dat.df$Kingdom == "Plantae"] <- "Plants"
dat.df$Taxon_group[dat.df$Phylum == "Chordata"] <- "Vertebrates" # Chordata in PREDICTS contains Aves, Mammalia, Reptilia, and Amphibia
dat.df$Taxon_group[dat.df$Kingdom == "Animalia" & dat.df$Phylum != "Chordata"] <- "Invertebrates"

## link sites within each other's radius of 5km
dat.df <- filter(dat.df, !is.na(Longitude), !is.na(Latitude))
dat.sf <- link_dist(dat.df, 5000)
dat.df <- left_join(dat.df, st_drop_geometry(dat.sf))
rm(dat.sf)
gc()

## aggregate site by Pressure_level
dat.df <- group_by(dat.df,
                   Source_ID, SS_ID, Site_ID,
                   Pressure_level,
                   Taxon_group, Taxon_number)
dat.df <- mutate(dat.df,
                 Longitude = mean(Longitude, na.rm = TRUE), Latitude = mean(Latitude, na.rm = TRUE),
                 Mean_measurement_t = mean(Measurement, na.rm = TRUE))
dat.df <- ungroup(dat.df)
dat.df <- select(dat.df,
                 Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_group, Taxon_number,
                 Pressure_level,
                 Mean_measurement_t)
dat.df <- distinct(dat.df)

## identify control measurements in natural conditions
control.df <- filter(dat.df, Pressure_level == "Primary vegetation - Minimal use")
control.df <- select(control.df, SS_ID, Taxon_number, Mean_measurement_t)
control.df <- rename(control.df, Mean_measurement_c = Mean_measurement_t)
control.df <- group_by(control.df, SS_ID, Taxon_number)
control.df <- summarise(control.df,
                        Mean_measurement_c = mean(Mean_measurement_c, na.rm = T))
control.df <- ungroup(control.df)
SS.v <- unique(pull(control.df, SS_ID))
dat.df <- filter(dat.df, SS_ID %in% SS.v)
rm(SS.v)

## complete data for species without treatment measurements
cols.v <- colnames(dat.df)
dat.df <- group_by(dat.df, Source_ID, SS_ID, Longitude, Latitude, Taxon_group)
dat.df <- complete(dat.df, Taxon_number, nesting(Pressure_level), fill = list(Mean_measurement_t = 0))
dat.df <- ungroup(dat.df)
dat.df <- relocate(dat.df, any_of(cols.v))
rm(cols.v)

## merge treatment and control data
dat.df <- filter(dat.df, Pressure_level != "Primary vegetation - Minimal use")
dat.df <- left_join(dat.df, control.df)
rm(control.df)
gc()


# calculate PDF and MSA
dat.df <- pdf_msa(dat.df)

## filter for measurements of more than 1 species
dat.df <- filter(dat.df, n_species > 1)

## filter, select, and sort relevant data
dat.df <- mutate(dat.df,
                 Database = "PREDICTS",
                 Pressure = "Land_use",
                 Pressure_unit = "Land_use_category",
                 MSA_loss = 1 - MSA)
dat.df <- select(dat.df,
                 Database, Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_group,
                 Pressure, Pressure_level, Pressure_unit,
                 n_species,
                 PDF,
                 MSA_loss)
dat.df <- distinct(dat.df)
dat.df <- arrange(dat.df, Source_ID, SS_ID, Taxon_group, Pressure_level)

# save data
write.csv(dat.df, "./Data/01_PREDICTS.csv", row.names = FALSE)


###### Kuipers_GCB_2023
# read data
dat.df <- read.csv("./Data/Agridiv_data.csv")

# modify data
## select columns and filter rows
dat.df <- select(dat.df,
                 Source_ID, SS_ID,
                 Longitude, Latitude,
                 Land_use, Group,
                 ITIS_TSN, Binomial,
                 Diversity_metric_type, Diversity_metric,
                 Diversity_value, Diversity_value_ref)

dat.df <- filter(dat.df,
                 Diversity_metric_type == "Abundance",
                 !is.na(Diversity_value_ref), Diversity_value_ref > 0,
                 !is.na(Binomial))

## rename columns
dat.df <- rename(dat.df,
                 Pressure_level = Land_use,
                 Taxon_number = ITIS_TSN,
                 Measurement_t = Diversity_value, Measurement_c = Diversity_value_ref)

## rename attributes
dat.df$Pressure_level[dat.df$Pressure_level == "Primary_vegetation" & dat.df$Group == "Reference"] <- "Primary vegetation - Minimal use"
dat.df$Pressure_level[dat.df$Pressure_level == "Primary_vegetation" & dat.df$Group != "Reference"] <- "Primary vegetation - Light/intense use"
dat.df$Pressure_level <- str_replace_all(dat.df$Pressure_level, "_", " ")

## link sites within each other's radius of 5km
dat.df <- filter(dat.df, !is.na(Longitude), !is.na(Latitude))
dat.sf <- link_dist(dat.df, 5000)
dat.df <- left_join(dat.df, st_drop_geometry(dat.sf))
rm(dat.sf)
gc()

## aggregate site by Pressure_level
dat.df <- group_by(dat.df,
                   Source_ID, SS_ID, Site_ID,
                   Pressure_level,
                   Taxon_number)
dat.df <- mutate(dat.df,
                 Longitude = mean(Longitude, na.rm = TRUE), Latitude = mean(Latitude, na.rm = TRUE),
                 Mean_measurement_t = mean(Measurement_t, na.rm = TRUE), Mean_measurement_c = mean(Measurement_c, na.rm = T))
dat.df <- ungroup(dat.df)
dat.df <- select(dat.df,
                 Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_number,
                 Pressure_level,
                 Mean_measurement_t, Mean_measurement_c)
dat.df <- distinct(dat.df)

## complete data for species without treatment measurements
cols.v <- colnames(dat.df)
dat.df <- group_by(dat.df, Source_ID, SS_ID, Longitude, Latitude)
dat.df <- complete(dat.df, Taxon_number, nesting(Pressure_level), fill = list(Mean_measurement_t = 0, Mean_measurement_c = 0))
dat.df <- ungroup(dat.df)
dat.df <- relocate(dat.df, any_of(cols.v))
rm(cols.v)

# calculate PDF and MSA
dat.df <- filter(dat.df, Pressure_level != "Primary vegetation - Minimal use")

dat.df <- pdf_msa(dat.df)

## filter for measurements of more than 1 species
dat.df <- filter(dat.df, n_species > 1)

## filter, select, and sort relevant data
dat.df <- mutate(dat.df,
                 Database = "Kuipers_GCB_2023",
                 Taxon_group = "Vertebrates",
                 Pressure = "Land_use",
                 Pressure_unit = "Land_use_category",
                 MSA_loss = 1 - MSA)
dat.df <- select(dat.df,
                 Database, Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_group,
                 Pressure, Pressure_level, Pressure_unit,
                 n_species,
                 PDF,
                 MSA_loss)
dat.df <- distinct(dat.df)
dat.df <- arrange(dat.df, Source_ID, SS_ID, Taxon_group, Pressure_level)

# save data
write.csv(dat.df, "./Data/01_Kuipers_GCB_2023.csv", row.names = FALSE)


###### Midolo_GEB_2019
# read data
dat.df <- read.csv("./Data/Midolo_etal_2018_data_N-application_FINAL_IA.csv")

# modify data
## select columns and filter rows
dat.df <- select(dat.df,
                 Experiment,
                 LONG, LAT,
                 Nadd,
                 Species,
                 Treatment, Control)

dat.df <- filter(dat.df,
                 !is.na(Control), Control > 0,
                 !is.na(Species)) # no effect for the dataset of midolo

## rename columns
dat.df <- rename(dat.df,
                 Source_ID = Experiment,
                 Longitude = LONG, Latitude = LAT,
                 Pressure_level = Nadd,
                 Taxon_number = Species,
                 Measurement_t = Treatment, Measurement_c = Control)

## add SS_ID column
dat.df <- mutate(dat.df, SS_ID = Source_ID)

## link sites within each other's radius of 5km [no effect for the dataset of midolo]
dat.df <- filter(dat.df, !is.na(Longitude), !is.na(Latitude))
dat.sf <- link_dist(dat.df, 5000)
dat.df <- left_join(dat.df, st_drop_geometry(dat.sf))
rm(dat.sf)
gc()

## aggregate site by Pressure_level
dat.df <- group_by(dat.df,
                   Source_ID, Site_ID,
                   Pressure_level,
                   Taxon_number)
dat.df <- mutate(dat.df,
                 Longitude = mean(Longitude, na.rm = TRUE), Latitude = mean(Latitude, na.rm = TRUE),
                 Mean_measurement_t = mean(Measurement_t, na.rm = TRUE), Mean_measurement_c = mean(Measurement_c, na.rm = T))
dat.df <- ungroup(dat.df)
dat.df <- select(dat.df,
                 Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_number,
                 Pressure_level,
                 Mean_measurement_t, Mean_measurement_c)
dat.df <- distinct(dat.df)

## complete data for species without treatment measurements [no effect for the dataset of Midolo]
cols.v <- colnames(dat.df)
dat.df <- group_by(dat.df, Source_ID, Longitude, Latitude)
dat.df <- complete(dat.df, Taxon_number, nesting(Pressure_level), fill = list(Mean_measurement_t = 0, Mean_measurement_c = 0))
dat.df <- ungroup(dat.df)
dat.df <- relocate(dat.df, any_of(cols.v))
rm(cols.v)

n.df <- read.csv("./Data/01_n.csv")
n.df$n_species_comparisons[n.df$Database == "Midolo_GEB_2019"] <- nrow(dat.df)

# calculate PDF and MSA
dat.df <- pdf_msa(dat.df)

## filter for measurements of more than 1 species
dat.df <- filter(dat.df, n_species > 1)

## filter, select, and sort relevant data
dat.df <- mutate(dat.df,
                 Database = "Midolo_GEB_2019",
                 Taxon_group = "Plants",
                 Pressure = "N_addition",
                 Pressure_unit = "kg_N_ha_yr",
                 MSA_loss = 1 - MSA)
dat.df <- select(dat.df,
                 Database, Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_group,
                 Pressure, Pressure_level, Pressure_unit,
                 n_species,
                 PDF,
                 MSA_loss)
dat.df <- distinct(dat.df)
dat.df <- arrange(dat.df, Source_ID, SS_ID, Taxon_group, Pressure_level)

# save data
write.csv(dat.df, "./Data/01_Midolo_GEB_2019.csv", row.names = FALSE)


###### BioTIME
met.df <- read.csv("./Data/BioTIMEMetadata_24_06_2021.csv")
dat.df <- read.csv("./Data/BioTIMEQuery_24_06_2021.csv")

# modify data
## select columns and filter rows
met.df <- filter(met.df,
                 ABUNDANCE_TYPE != "Presence/Absence",
                 !is.na(ABUNDANCE_TYPE), REALM == "Terrestrial",
                 !(TAXA %in% c("All", "Fungi")))
met.df <- select(met.df, STUDY_ID, TAXA)

dat.df <- select(dat.df,
                 STUDY_ID, SAMPLE_DESC,
                 LATITUDE, LONGITUDE,
                 YEAR,
                 ID_SPECIES, GENUS_SPECIES,
                 sum.allrawdata.ABUNDANCE)
dat.df <- filter(dat.df,
                 STUDY_ID %in% pull(met.df, STUDY_ID),
                 !is.na(GENUS_SPECIES))
dat.df <- left_join(dat.df, met.df)

dat.df <- mutate(dat.df,
                 SS_ID = paste(STUDY_ID,
                               word(SAMPLE_DESC, start = 3, end = str_count(SAMPLE_DESC, "_") + 1, sep = "_"),
                               sep = "_"))
dat.df <- select(dat.df, -SAMPLE_DESC)
rm(met.df)

## rename columns
dat.df <- rename(dat.df,
                 Source_ID = STUDY_ID,
                 Latitude = LATITUDE, Longitude = LONGITUDE,
                 Pressure_level = YEAR,
                 Taxon_group = TAXA, Taxon_number = ID_SPECIES, Species = GENUS_SPECIES,
                 Mean_measurement = sum.allrawdata.ABUNDANCE)

## rename attributes
dat.df$Taxon_group[dat.df$Taxon_group == "Terrestrial plants"] <- "Plants"
dat.df$Taxon_group[dat.df$Taxon_group == "Terrestrial invertebrates"] <- "Invertebrates"
dat.df$Taxon_group[dat.df$Taxon_group %in% c("Birds", "Mammals", "Reptiles", "Amphibians")] <- "Vertebrates"

dat.df <- select(dat.df, Source_ID, SS_ID, Latitude, Longitude, Taxon_group, Taxon_number, Species, Pressure_level, Mean_measurement)
dat.df <- arrange(dat.df, SS_ID, Taxon_group, Taxon_number, Pressure_level)

## link sites within each other's radius of 5km
dat.df <- filter(dat.df, !is.na(Longitude), !is.na(Latitude))
dat.df <- mutate(dat.df,
                 Longitude = round(Longitude, 2),
                 Latitude = round(Latitude, 2))
dat.df <- group_by(dat.df,
                   SS_ID,
                   Latitude, Longitude,
                   Taxon_number,
                   Pressure_level)
dat.df <- mutate(dat.df, Mean_measurement = mean(Mean_measurement, na.rm = T))
dat.df <- ungroup(dat.df)
dat.df <- distinct(dat.df)

dat.sf <- link_dist2(dat.df, 5000)
dat.df <- left_join(dat.df, st_drop_geometry(dat.sf))
rm(dat.sf)
gc()

## aggregate site by Pressure_level
dat.df <- group_by(dat.df,
                   Source_ID, SS_ID, Site_ID,
                   Pressure_level,
                   Taxon_group, Taxon_number)
dat.df <- mutate(dat.df,
                 Longitude = mean(Longitude, na.rm = TRUE), Latitude = mean(Latitude, na.rm = TRUE),
                 Mean_measurement_t = mean(Mean_measurement, na.rm = TRUE))
dat.df <- ungroup(dat.df)
dat.df <- select(dat.df, -Mean_measurement)
dat.df <- distinct(dat.df)

dat.df <- mutate(dat.df, Year = Pressure_level)
dat.df <- group_by(dat.df,
                   Source_ID, SS_ID, Site_ID,
                   Taxon_group, Taxon_number)
dat.df <- mutate(dat.df,
                 Pressure_level = Year - min(Year),
                 Mean_measurement_c = Mean_measurement_t[which.min(Year)])
dat.df <- ungroup(dat.df)

dat.df <- select(dat.df,
                 Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_group, Taxon_number,
                 Pressure_level,
                 Mean_measurement_t, Mean_measurement_c)

## complete data for species without treatment measurements
cols.v <- colnames(dat.df)
dat.df <- group_by(dat.df, Source_ID, SS_ID, Longitude, Latitude, Taxon_group)
dat.df <- complete(dat.df, Taxon_number, nesting(Pressure_level), fill = list(Mean_measurement_t = 0, Mean_measurement_c = 0))
dat.df <- ungroup(dat.df)
dat.df <- relocate(dat.df, any_of(cols.v))
rm(cols.v)

# calculate PDF and MSA
dat.df <- pdf_msa(dat.df)

## filter for measurements of more than 1 species
dat.df <- filter(dat.df,
                 n_species > 1,
                 Mean_measurement_c > 0,
                 Pressure_level > 0)

## filter, select, and sort relevant data
dat.df <- mutate(dat.df,
                 Database = "BioTIME",
                 Pressure = "Time_as_proxy",
                 Pressure_unit = "Years_since_reference_measure",
                 MSA_loss = 1 - MSA)
dat.df <- select(dat.df,
                 Database, Source_ID, SS_ID,
                 Longitude, Latitude,
                 Taxon_group,
                 Pressure, Pressure_level, Pressure_unit,
                 n_species,
                 PDF,
                 MSA_loss)
dat.df <- distinct(dat.df)
dat.df <- arrange(dat.df, Source_ID, SS_ID, Taxon_group, Pressure_level)

# save data
write.csv(dat.df, "./Data/01_BioTIME.csv", row.names = FALSE)

rm(list=ls())
gc()


###### Combine data
# read data
pred.df <- read.csv("./Data/01_PREDICTS.csv")
kuip.df <- read.csv("./Data/01_Kuipers_GCB_2023.csv")
mido.df <- read.csv("./Data/01_Midolo_GEB_2019.csv")
biot.df <- read.csv("./Data/01_BioTIME.csv")

pred.df <- mutate(pred.df,
                  Source_ID = as.character(Source_ID),
                  Pressure_level = as.character(Pressure_level))
kuip.df <- mutate(kuip.df,
                  Source_ID = as.character(Source_ID),
                  Pressure_level = as.character(Pressure_level))
mido.df <- mutate(mido.df,
                  Source_ID = as.character(Source_ID),
                  Pressure_level = as.character(Pressure_level))
biot.df <- mutate(biot.df,
                  Source_ID = as.character(Source_ID),
                  Pressure_level = as.character(Pressure_level))

dat.df <- bind_rows(pred.df, kuip.df, mido.df, biot.df)
rm(pred.df, kuip.df, mido.df, biot.df)

dat.df <- mutate(dat.df, Source_ID_original = Source_ID, SS_ID_original = SS_ID)
dat.df <- mutate(dat.df, Source_ID = match(paste(Database, Source_ID_original), unique(paste(Database, Source_ID_original))))
dat.df <- group_by(dat.df, Database, Source_ID)
dat.df <- mutate(dat.df, SS_ID = match(SS_ID_original, unique(SS_ID_original)))
dat.df <- mutate(dat.df, SS_ID = paste(Source_ID, SS_ID, sep = "_"))
dat.df <- ungroup(dat.df)
dat.df <- mutate(dat.df, as.numeric(Source_ID))

dat.df <- select(dat.df,
                 Database, Source_ID, SS_ID, Source_ID_original, SS_ID_original,
                 Longitude, Latitude,
                 Taxon_group,
                 Pressure, Pressure_level, Pressure_unit,
                 n_species,
                 PDF, MSA_loss)

write.csv(dat.df, "./Data/01_PDF_MSA.csv", row.names = FALSE)