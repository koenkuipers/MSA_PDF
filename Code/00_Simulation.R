# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# Libraries
library(ggplot2)
library(egg)
library(forcats)
library(dplyr)
library(tidyr)
library(stringr)

# Create data
## Population data
dat.df <- data.frame(population_id = rep(rep(1:1000, each = 101), 4),
                     pressure_intensity = rep(seq(0, 1, 0.01), 1000 * 4),
                     intraspecies_diversity = rep(c("min", "max"), each = (1000 * 101 * 4) / 2),
                     interspecies_diversity = rep(rep(c("min", "max"), each = (1000 * 101 * 4) / 4), 2),
                     abundance = NA)

dat.df <- data.frame(population_id = rep(rep(1:1000, each = 1001), 4),
                     pressure_intensity = rep(seq(0, 1, 0.001), 1000 * 4),
                     intraspecies_diversity = rep(c("min", "max"), each = (1000 * 1001 * 4) / 2),
                     interspecies_diversity = rep(rep(c("min", "max"), each = (1000 * 1001 * 4) / 4), 2),
                     abundance = NA)

dat.df[which(dat.df$intraspecies_diversity == "min" & dat.df$interspecies_diversity == "min" & dat.df$pressure_intensity < 0.5), "abundance"] <- 1
dat.df[which(dat.df$intraspecies_diversity == "min" & dat.df$interspecies_diversity == "min" & dat.df$pressure_intensity >= 0.5), "abundance"] <- 0
dat.df[which(dat.df$intraspecies_diversity == "max" & dat.df$interspecies_diversity == "min"), "abundance"] <- 1 - dat.df[which(dat.df$intraspecies_diversity == "max" & dat.df$interspecies_diversity == "min"), "pressure_intensity"]
dat.df[which(dat.df$intraspecies_diversity == "min" & dat.df$interspecies_diversity == "max"), "abundance"] <- 1

dat.df[which(dat.df$intraspecies_diversity == "min"), "intraspecies_diversity"] <- "No intraspecies variability"
dat.df[which(dat.df$intraspecies_diversity == "max"), "intraspecies_diversity"] <- "Intraspecies variability"
dat.df[which(dat.df$interspecies_diversity == "min"), "interspecies_diversity"] <- "No interspecies variability"
dat.df[which(dat.df$interspecies_diversity == "max"), "interspecies_diversity"] <- "Interspecies variability"

for(i in 1:1000){
  dat.df[which(dat.df$intraspecies_diversity == "min" & dat.df$interspecies_diversity == "max" & dat.df$population_id == i & dat.df$pressure_intensity >= i / 1000), "abundance"] <- 0
}
rm(i)

for(i in 1:1000){
  dat.df[which(dat.df$intraspecies_diversity == "max" & dat.df$interspecies_diversity == "max" & dat.df$population_id == i), "abundance"] <- 1 + -(1 / (1 - dat.df[which(dat.df$intraspecies_diversity == "max" & dat.df$interspecies_diversity == "max" & dat.df$population_id == i), "pressure_intensity"][i])) * dat.df[which(dat.df$intraspecies_diversity == "max" & dat.df$interspecies_diversity == "max" & dat.df$population_id == i), "pressure_intensity"]
}
rm(i)

dat.df[which(dat.df$abundance < 0), "abundance"] <- 0

## Diversity data
div.df <- dat.df[which(dat.df$population_id %in% c(1, seq(10, 1000, 10))), ]
div.df$intraspecies_diversity <- fct_rev(as.factor(div.df$intraspecies_diversity))
div.df$interspecies_diversity <- fct_rev(as.factor(div.df$interspecies_diversity))

wid.df <- pivot_wider(div.df, names_from = population_id, values_from = abundance)

div.df <- unique(div.df[, c("pressure_intensity", "intraspecies_diversity", "interspecies_diversity")])
div.df$abundance <- rowMeans(wid.df[, 4:ncol(wid.df)])
div.df$richness <- rowSums(wid.df[, 4:ncol(wid.df)] > 0) / length(4:ncol(wid.df))
div.df <- pivot_longer(div.df, abundance:richness, names_to = "diversity_metric", values_to = "diversity_value")
div.df$diversity_value <-  1 - div.df$diversity_value
div.df$diversity_metric <- str_to_sentence((div.df$diversity_metric))
div.df[which(div.df$diversity_metric == "Richness"), "diversity_metric"] <- "PDF"
div.df[which(div.df$diversity_metric == "Abundance"), "diversity_metric"] <- "MSA-loss"

# Write data
write.csv(dat.df, file = "./Data/00_Population_simulation.csv", row.names = F)
write.csv(div.df, file = "./Data/00_Diversity_simulation.csv", row.names = F)