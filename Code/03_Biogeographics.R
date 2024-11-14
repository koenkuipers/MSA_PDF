# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# libraries
library(dplyr)
library(sf)

sf_use_s2(use_s2 = F)

# Data
## All external data used is publicly available:
### Ecoregions: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world

## We have downloaded the data and stored the csv datasets in our 'Data' folder with the following (original) names:
### Ecoregions: wwf_terr_ecos.shp

dat.df <- read.csv("./Data/01_PDF_MSA.csv")
eco.sf <- st_read("./Data/wwf_terr_ecos.shp")
rea.df <- read.csv("./Data/Realm_legend.csv")
bio.df <- read.csv("./Data/Biomes_legend.csv")

# modify data
dat.sf <- st_as_sf(dat.df, coords = c("Longitude", "Latitude"), remove = FALSE, crs = st_crs(4326))
eco.sf <- select(eco.sf, REALM, BIOME)
eco.sf <- rename(eco.sf, Realm_ID = REALM, Biome_ID = BIOME)
rea.df$Realm_ID[is.na(rea.df$Realm_ID)] <- "NA"
bio.df <- select(bio.df, Biome_ID, Biome_name)

# match ecoregion data to dat.sf
dat.sf <- st_intersection(dat.sf, eco.sf)

# add realm and biome to dat.df
dat.df <- left_join(dat.df, st_drop_geometry(dat.sf))

# add realm and biome names
dat.df <- left_join(dat.df, rea.df)
dat.df <- left_join(dat.df, bio.df)
dat.df <- rename(dat.df,
                 Realm = Realm_name,
                 Biome = Biome_name)
dat.df <- select(dat.df,
                 Database, Source_ID, SS_ID, Source_ID_original, SS_ID_original,
                 Longitude, Latitude,
                 Realm, Biome,
                 Taxon_group,
                 Pressure, Pressure_level, Pressure_unit,
                 n_species,
                 PDF, MSA_loss)

write.csv(dat.df, "./Data/03_PDF_MSA.csv", row.names = FALSE)
