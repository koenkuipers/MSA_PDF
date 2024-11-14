# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# Libraries
library(dplyr)
library(ggplot2)
library(lemon)
library(sf)
library(spData)
library(egg)
library(forcats)
library(grid)

# Data
pop.df <- read.csv("./Data/00_Population_simulation.csv")
div.df <- read.csv("./Data/00_Diversity_simulation.csv")

dat.df <- read.csv("./Data/03_PDF_MSA.csv")

pdf.df <- read.csv("./Data/05_Regression_PDF_figure_data.csv")
taxon.df <- read.csv("./Data/05_Regression_PDF_Taxon_figure_data.csv")
realm.df <- read.csv("./Data/05_Regression_PDF_Realm_figure_data.csv")

pdf10.df <- read.csv("./Data/05_Regression_PDF_n10_figure_data.csv")
taxon10.df <- read.csv("./Data/05_Regression_PDF_Taxon_n10_figure_data.csv")
realm10.df <- read.csv("./Data/05_Regression_PDF_Realm_n10_figure_data.csv")

taxon.df <- mutate(taxon.df, Taxon_group = fct_relevel(Taxon_group, "Plants", "Invertebrates"))
taxon10.df <- mutate(taxon10.df, Taxon_group = fct_relevel(Taxon_group, "Plants", "Invertebrates"))
dat.df <- mutate(dat.df, Taxon_group = fct_relevel(Taxon_group, "Plants", "Invertebrates"))

# Filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species > 1)
dat10.df <- filter(dat.df, !is.na(Realm), n_species >= 10)

# Manipulate data
dat.df <- mutate(dat.df, panel = "A)")
dat10.df <- mutate(dat10.df, panel = "B)")
dat3.df <- bind_rows(dat.df, dat10.df)

pdf.df <- mutate(pdf.df, panel = "A)")
pdf10.df <- mutate(pdf10.df, panel = "B)")
pdf3.df <- bind_rows(pdf.df, pdf10.df)

# colours
col.chr <- "#4969E1"
col2.chr <- colorspace::lighten("#4969E1", 0.5)


# Fig 1. Simulation -------------------------------------------------------

pop.df <- pop.df[which(pop.df$population_id %in% c(1, seq(10, 1000, 10))), ]
pop.df$intraspecies_diversity <- gsub(" va", "\nva", pop.df$intraspecies_diversity)
pop.df$interspecies_diversity <- gsub(" va", "\nva", pop.df$interspecies_diversity)
pop.df$intraspecies_diversity <- fct_rev(as.factor(pop.df$intraspecies_diversity))

div.df$intraspecies_diversity <- gsub("Minimum", "No", div.df$intraspecies_diversity)
div.df$intraspecies_diversity <- gsub("Maximum i", "I", div.df$intraspecies_diversity)
div.df$intraspecies_diversity <- gsub("s ", "s\n", div.df$intraspecies_diversity)
div.df$interspecies_diversity <- gsub("Minimum", "No", div.df$interspecies_diversity)
div.df$interspecies_diversity <- gsub("Maximum i", "I", div.df$interspecies_diversity)
div.df$interspecies_diversity <- gsub("s ", "s\n", div.df$interspecies_diversity)
div.df$intraspecies_diversity <- gsub("diversity", "variability", div.df$intraspecies_diversity)
div.df$interspecies_diversity <- gsub("diversity", "variability", div.df$interspecies_diversity)
div.df$intraspecies_diversity <- fct_rev(as.factor(div.df$intraspecies_diversity))

wid.df <- pivot_wider(div.df, names_from = diversity_metric, values_from = diversity_value)
wid.df$relationship <- wid.df$`PDF` / wid.df$`MSA-loss`
wid.df <- unique(wid.df[, c(1:3, ncol(wid.df))])
wid.df[which(is.na(wid.df$relationship)), "relationship"] <- 1

wid2.df <- pivot_wider(div.df, names_from = diversity_metric, values_from = diversity_value)
wid2.df$relationship <- wid2.df$`PDF` / wid2.df$`MSA-loss`
wid2.df[which(is.na(wid2.df$relationship)), "relationship"] <- 1
wid2.df$panel <- "A"
wid2.df$panel[1002:(1002+1000)] <- "B"
wid2.df$panel[2003:(2003+1000)] <- "C"
wid2.df$panel[3004:(3004+1000)] <- "D"

step.n <- 50
pop.v <- c(1, seq(step.n, 1000, step.n))

p0 <- ggplot(data = filter(pop.df, population_id %in% pop.v), mapping = aes(x = pressure_intensity, y = 1 - abundance, group = population_id)) +
  geom_line(alpha = 0.1) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  coord_fixed() +
  ylab("Relative population abundance loss") +
  xlab("Relative pressure intensity") +
  facet_rep_grid(interspecies_diversity ~ intraspecies_diversity) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        text = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        strip.text = element_text(size = 7, color = "black"),
        strip.background = element_blank())

p1 <- ggplot(data = wid.df, mapping = aes(x = `PDF`, y = `MSA-loss`)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  coord_fixed() +
  facet_rep_grid(interspecies_diversity ~ intraspecies_diversity) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        strip.text = element_text(size = 7, color = "black"),
        legend.position = c(0.15, 0.1))

p2 <- ggplot(data = wid2.df, mapping = aes(x = `PDF`, y = `MSA-loss`)) +
  geom_ribbon(aes(`PDF`, ymin = `MSA-loss`, ymax = 1), fill = "lightgrey") +
  geom_ribbon(data = wid2.df[wid.df$panel == "D", ], aes(`PDF`, ymin = `MSA-loss`, ymax = `PDF`), fill = "grey") +
  geom_line(aes(col = panel), linewidth = 1) +
  scale_colour_manual(values = c("#000000", "#000000", "#000000", "#000000")) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  coord_fixed() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 7, color = "black"),
        strip.text = element_text(size = 7, color = "black"),
        legend.position = "none")

f1 <- ggarrange(p0, p1, p2, ncol = 3, labels = c("A)", "B)", "C)"), label.args = list(gp = gpar(font = 1, cex = 0.6)))
ggsave(f1, filename = "./Figures/Figure_1.pdf", width = 63*3, height = 70, units = "mm")
ggsave(f1, filename = "./Figures/Figure_1.png", width = 63*3, height = 70, units = "mm")


# Fig 2. Map --------------------------------------------------------------

# spatial data
dat.sf <- st_as_sf(dat.df, coords = c("Longitude", "Latitude"), crs = st_crs(4326))
dat.sf <- st_transform(dat.sf, crs = st_crs("ESRI:54012"))
world.sf <- spData::world
world.sf <- st_transform(world.sf, crs = st_crs("ESRI:54012"))

# plot
f2 <- ggplot(data = dat.sf) +
  geom_sf(data = world.sf, col = "lightgrey", fill = "lightgrey") +
  geom_sf(col = "#4969E1", alpha = 0.25) +
  theme(panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"))

ggsave(f2, filename = "./Figures/Figure_2.pdf", width = 90, height = 50, units = "mm")
ggsave(f2, filename = "./Figures/Figure_2.png", width = 90, height = 50, units = "mm")

# Fig. 3 MSA~PDF ----------------------------------------------------------

f3 <- ggplot(data = pdf3.df, mapping = aes(x = PDF, y = MSA_loss_pred_med)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(data = dat3.df, aes(x = PDF, y = MSA_loss, size = sqrt(n_species)), alpha = 0.2, colour = col2.chr) +
  geom_ribbon(aes(ymin = MSA_loss_pred_lwr, ymax = MSA_loss_pred_upr), alpha = 0.5, fill = col.chr) +
  geom_line(lwd = 1, colour = "white") +
  facet_rep_wrap(panel ~ ., ncol = 3) +
  scale_size(range = c(1, 10)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1)) +
  ylab("MSA loss") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, colour = "black", hjust = 0),
        axis.line = element_line(),
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "none",
        aspect.ratio = 1)

ggsave(f3, filename = "./Figures/Figure_3.pdf", width = 63*2, height = 70, units = "mm")
ggsave(f3, filename = "./Figures/Figure_3.png", width = 63*2, height = 70, units = "mm")


# Fig. 4 MSA~PDF*Taxon ----------------------------------------------------

taxcol.v <- c(3, 1, 2)

f4 <- ggplot(data = taxon.df, mapping = aes(x = PDF, y = MSA_loss_pred_med, col = Taxon_group, fill = Taxon_group)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(data = dat.df, aes(x = PDF, y = MSA_loss, size = sqrt(n_species)), alpha = 0.2) +
  geom_ribbon(aes(ymin = MSA_loss_pred_lwr, ymax = MSA_loss_pred_upr), alpha = 0.5) +
  geom_line(lwd = 1, colour = "white") +
  facet_rep_wrap(Taxon_group ~ ., ncol = 3) +
  scale_size(range = c(1, 10)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), clip = "off") +
  ylab("MSA loss") +
  scale_color_manual(values = hcl.colors(6, "pastel 1")[taxcol.v]) +
  scale_fill_manual(values = hcl.colors(6, "dark 3")[taxcol.v]) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, colour = "black"),
        axis.line = element_line(),
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "none",
        aspect.ratio = 1)

ggsave(f4, filename = "./Figures/Figure_4.pdf", width = 63*3, height = 70, units = "mm")
ggsave(f4, filename = "./Figures/Figure_4.png", width = 63*3, height = 70, units = "mm")



# Fig. 5 MSA~PDF*Realm ----------------------------------------------------

f5 <- ggplot(data = filter(realm.df, Realm != "Oceania"), mapping = aes(x = PDF, y = MSA_loss_pred_med, col = Realm, fill = Realm)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(data = filter(dat.df, Realm != "Oceania"), aes(x = PDF, y = MSA_loss, size = sqrt(n_species)), alpha = 0.2) +
  geom_ribbon(aes(ymin = MSA_loss_pred_lwr, ymax = MSA_loss_pred_upr), alpha = 0.5) +
  geom_line(lwd = 1, colour = "white") +
  facet_rep_wrap(Realm ~ .) +
  scale_size(range = c(1, 10)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), clip = "off") +
  ylab("MSA loss") +
  scale_color_manual(values = hcl.colors(6, "pastel 1")) +
  scale_fill_manual(values = hcl.colors(6, "dark 3")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, colour = "black"),
        axis.line = element_line(),
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "none",
        aspect.ratio = 1)

ggsave(f5, filename = "./Figures/Figure_5.pdf", width = 63*3, height = 70*2, units = "mm")
ggsave(f5, filename = "./Figures/Figure_5.png", width = 63*3, height = 70*2, units = "mm")


# Fig. S1 MSA~PDF*Taxon (n>10) --------------------------------------------

taxcol.v <- c(3, 1, 2)

fs1 <- ggplot(data = taxon10.df, mapping = aes(x = PDF, y = MSA_loss_pred_med, col = Taxon_group, fill = Taxon_group)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(data = dat10.df, aes(x = PDF, y = MSA_loss, size = sqrt(n_species)), alpha = 0.2) +
  geom_ribbon(aes(ymin = MSA_loss_pred_lwr, ymax = MSA_loss_pred_upr), alpha = 0.5) +
  geom_line(lwd = 1, colour = "white") +
  facet_rep_wrap(Taxon_group ~ ., ncol = 3) +
  scale_size(range = c(1, 10)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), clip = "off") +
  ylab("MSA loss") +
  scale_color_manual(values = hcl.colors(6, "pastel 1")[taxcol.v]) +
  scale_fill_manual(values = hcl.colors(6, "dark 3")[taxcol.v]) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, colour = "black"),
        axis.line = element_line(),
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "none",
        aspect.ratio = 1)

ggsave(fs1, filename = "./Figures/Figure_S1.pdf", width = 63*3, height = 70, units = "mm")
ggsave(fs1, filename = "./Figures/Figure_S1.png", width = 63*3, height = 70, units = "mm")


# Fig. S2 MSA~PDF*Realm (n>10) --------------------------------------------

fs2 <- ggplot(data = filter(realm10.df, Realm != "Oceania"), mapping = aes(x = PDF, y = MSA_loss_pred_med, col = Realm, fill = Realm)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(data = filter(dat10.df, Realm != "Oceania"), aes(x = PDF, y = MSA_loss, size = sqrt(n_species)), alpha = 0.2) +
  geom_ribbon(aes(ymin = MSA_loss_pred_lwr, ymax = MSA_loss_pred_upr), alpha = 0.5) +
  geom_line(lwd = 1, colour = "white") +
  facet_rep_wrap(Realm ~ .) +
  scale_size(range = c(1, 10)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), clip = "off") +
  ylab("MSA loss") +
  scale_color_manual(values = hcl.colors(6, "pastel 1")) +
  scale_fill_manual(values = hcl.colors(6, "dark 3")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, colour = "black"),
        axis.line = element_line(),
        text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "none",
        aspect.ratio = 1)

ggsave(fs2, filename = "./Figures/Figure_S2.pdf", width = 63*3, height = 70*2, units = "mm")
ggsave(fs2, filename = "./Figures/Figure_S2.png", width = 63*3, height = 70*2, units = "mm")