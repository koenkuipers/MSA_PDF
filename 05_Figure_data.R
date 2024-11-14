# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# libraries
library(languageserver)
library(httpgd)
library(dplyr)
library(lme4)
library(glmmTMB)

# functions
pred.f <- function(pred.df) {
  # simulation
  b.bm <- lme4::bootMer(mod.glm, nsim = 1000, FUN = function(x){predict(x, newdata = pred.df, re.form = NA, allow.new.levels = T)})
  b.bm$t <- b.bm$t[complete.cases(b.bm$t), ]

  pred.df <- mutate(pred.df, logitMSA_loss_pred = predict(mod.glm, pred.df, re.form = NA, allow.new.levels = T),
                    logitMSA_loss_pred_lwr = apply(b.bm$t, 2, FUN = function(x){quantile(x, 0.025)}),
                    logitMSA_loss_pred_med = apply(b.bm$t, 2, FUN = function(x){median(x)}),
                    logitMSA_loss_pred_mea = apply(b.bm$t, 2, FUN = function(x){mean(x)}),
                    logitMSA_loss_pred_upr = apply(b.bm$t, 2, FUN = function(x){quantile(x, 0.975)}))
  pred.df <- mutate(pred.df, MSA_loss_pred = exp(logitMSA_loss_pred) / (1 + exp(logitMSA_loss_pred)),
                    MSA_loss_pred_lwr = exp(logitMSA_loss_pred_lwr) / (1 + exp(logitMSA_loss_pred_lwr)),
                    MSA_loss_pred_med = exp(logitMSA_loss_pred_med) / (1 + exp(logitMSA_loss_pred_med)),
                    MSA_loss_pred_mea = exp(logitMSA_loss_pred_mea) / (1 + exp(logitMSA_loss_pred_mea)),
                    MSA_loss_pred_upr = exp(logitMSA_loss_pred_upr) / (1 + exp(logitMSA_loss_pred_upr)))
  
  # replace MSA-loss>PDF
  pred.df$MSA_loss_pred[pred.df$MSA_loss_pred < pred.df$PDF] <- NA
  pred.df$MSA_loss_pred_lwr[pred.df$MSA_loss_pred_lwr < pred.df$PDF] <- pred.df$PDF[pred.df$MSA_loss_pred_lwr < pred.df$PDF]
  pred.df$MSA_loss_pred_med[pred.df$MSA_loss_pred_med < pred.df$PDF] <- NA
  pred.df$MSA_loss_pred_mea[pred.df$MSA_loss_pred_mea < pred.df$PDF] <- NA
  pred.df$MSA_loss_pred_upr[pred.df$MSA_loss_pred_upr < pred.df$PDF] <- NA
  pred.df$MSA_loss_pred_lwr[is.na(pred.df$MSA_loss_pred_upr)] <- NA
  
  pred.df
}

# MSA_loss ~ PDF ----------------------------------------------------------

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species > 1)

# regression models
mod.glm <- glmmTMB(MSA_loss ~ 1 + PDF + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# prepare dataframe of the regression curve
pred.df <- data.frame(Source_ID = 0,
                      PDF = seq(0, 1, 0.01),
                      n_species = 1)

pred.df <- pred.f(pred.df)

write.csv(pred.df, "./Data/05_Regression_PDF_figure_data.csv", row.names = FALSE)

rm(dat.df, mod.glm, pred.df)
gc()


# MSA-loss ~ PDF * Taxon --------------------------------------------------

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species > 1)

# regression models
mod.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Taxon_group + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# prepare dataframe of the regression curve
pred.df <- data.frame(Source_ID = 0,
                      PDF = rep(seq(0, 1, 0.01), length(unique(dat.df$Taxon_group))),
                      Taxon_group = rep(sort(unique(dat.df$Taxon_group)), each = length(seq(0, 1, 0.01))),
                      n_species = 1)

pred.df <- pred.f(pred.df)

write.csv(pred.df, "./Data/05_Regression_PDF_Taxon_figure_data.csv", row.names = FALSE)

rm(dat.df, mod.glm, pred.df)
gc()


# MSA-loss ~ PDF * Realm --------------------------------------------------

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species > 1)

# regression models
mod.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Realm + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# prepare dataframe of the regression curve
pred.df <- data.frame(Source_ID = 0,
                      PDF = rep(seq(0, 1, 0.01), length(unique(dat.df$Realm))),
                      Realm = rep(sort(unique(dat.df$Realm)), each = length(seq(0, 1, 0.01))),
                      n_species = 1)

pred.df <- pred.f(pred.df)

write.csv(pred.df, "./Data/05_Regression_PDF_Realm_figure_data.csv", row.names = FALSE)

rm(dat.df, mod.glm, pred.df)
gc()


# MSA_loss ~ PDF n>=10 ----------------------------------------------------

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species >= 10)

# regression models
mod.glm <- glmmTMB(MSA_loss ~ 1 + PDF + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# prepare dataframe of the regression curve
pred.df <- data.frame(Source_ID = 0,
                      PDF = seq(0, 1, 0.01),
                      n_species = 1)

pred.df <- pred.f(pred.df)

write.csv(pred.df, "./Data/05_Regression_PDF_n10_figure_data.csv", row.names = FALSE)
rm(dat.df, mod.glm, pred.df)
gc()


# MSA-loss ~ PDF * Taxon n>=10 --------------------------------------------

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species >= 10)

# regression models
mod.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Taxon_group + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# prepare dataframe of the regression curve
pred.df <- data.frame(Source_ID = 0,
                      PDF = rep(seq(0, 1, 0.01), length(unique(dat.df$Taxon_group))),
                      Taxon_group = rep(sort(unique(dat.df$Taxon_group)), each = length(seq(0, 1, 0.01))),
                      n_species = 1)

pred.df <- pred.f(pred.df)

write.csv(pred.df, "./Data/05_Regression_PDF_Taxon_n10_figure_data.csv", row.names = FALSE)

rm(dat.df, mod.glm, pred.df)
gc()


# MSA-loss ~ PDF * Realm n>=10 --------------------------------------------

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species >= 10, Realm != "Oceania")

# regression models
mod.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Realm + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# prepare dataframe of the regression curve
pred.df <- data.frame(Source_ID = 0,
                      PDF = rep(seq(0, 1, 0.01), length(unique(dat.df$Realm))),
                      Realm = rep(sort(unique(dat.df$Realm)), each = length(seq(0, 1, 0.01))),
                      n_species = 1)

pred.df <- pred.f(pred.df)

write.csv(pred.df, "./Data/05_Regression_PDF_Realm_n10_figure_data.csv", row.names = FALSE)

rm(dat.df, mod.glm, pred.df)
gc()
