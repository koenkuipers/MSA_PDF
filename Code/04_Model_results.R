# Author: KJJ Kuipers
# Date: 20/09/2024
# Code for: Kuipers, KJJ, A Melki, S Morel, AM Schipper (2024) Relationships between mean species abundance (MSA) and potentially disappeared fraction of species (PDF) are consistent butt also uncertain. Ecological Indicators.

# Directory
setwd("Drive:/folderpath/")

# functions
response.f <- function(model, data) {
  mod.df <- data.frame(Variable = row.names(coef(summary(model))$cond),
                       Element = "Intercept",
                       Estimate = round(coef(summary(model))$cond[, 1], 2),
                       lwr = round(confint(model)[1:(nrow(confint(model))-1), 1], 2),
                       upr = round(confint(model)[1:(nrow(confint(model))-1), 2], 2),
                       se = round(coef(summary(model))$cond[, 2], 2),
                       p = round(coef(summary(model))$cond[, 4], 4))
  mod.df$Variable <- stringr::str_replace_all(mod.df$Variable, "\\(", "")
  mod.df$Variable <- stringr::str_replace_all(mod.df$Variable, "\\)", "")
  mod.df$Element[str_which(mod.df$Variable, "PDF")] <- "Slope"
  
  if(any(stringr::str_detect(mod.df$Variable, "Taxon_group")) == TRUE) {
    
     mod.df$Variable <- stringr::str_replace_all(mod.df$Variable, ":Taxon_group", ":")
    
     mod.df$Estimate[str_which(mod.df$Variable, "Taxon_group")] <- mod.df$Estimate[str_which(mod.df$Variable, "Taxon_group")] + mod.df$Estimate[1]
     mod.df$lwr[str_which(mod.df$Variable, "Taxon_group")] <- mod.df$lwr[str_which(mod.df$Variable, "Taxon_group")] + mod.df$Estimate[1]
     mod.df$upr[str_which(mod.df$Variable, "Taxon_group")] <- mod.df$upr[str_which(mod.df$Variable, "Taxon_group")] + mod.df$Estimate[1]

     mod.df$Estimate[str_which(mod.df$Variable, "PDF:")] <- mod.df$Estimate[str_which(mod.df$Variable, "PDF:")] + mod.df$Estimate[2]
     mod.df$lwr[str_which(mod.df$Variable, "PDF:")] <- mod.df$lwr[str_which(mod.df$Variable, "PDF:")] + mod.df$Estimate[2]
     mod.df$upr[str_which(mod.df$Variable, "PDF:")] <- mod.df$upr[str_which(mod.df$Variable, "PDF:")] + mod.df$Estimate[2]

     mod.df$Variable <- stringr::str_replace_all(mod.df$Variable, "Taxon_group", "")

    if(any(stringr::str_detect(mod.df$Variable, "Invertebrates")) == FALSE) {
      mod.df$Variable[1:2] <- c("Invertebrates", "PDF:Invertebrates")
    } else if(any(stringr::str_detect(mod.df$Variable, "Plants")) == FALSE) {
      mod.df$Variable[1:2] <- c("Plants", "PDF:Plants")
    } else if(any(stringr::str_detect(mod.df$Variable, "Vertebrates")) == FALSE) {
      mod.df$Variable[1:2] <- c("Vertebrates", "PDF:Vertebrates")
    }
    n.df <- data.frame("Variable" = names(table(data$Taxon_group)),
                       "n" = as.vector(table(data$Taxon_group)))
    mod.df <- left_join(mod.df, n.df)
  }
  
  if(any(stringr::str_detect(mod.df$Variable, "Realm")) == TRUE) {
     
     mod.df$Variable <- stringr::str_replace_all(mod.df$Variable, ":Realm", ":")
    
     mod.df$Estimate[str_which(mod.df$Variable, "Realm")] <- mod.df$Estimate[str_which(mod.df$Variable, "Realm")] + mod.df$Estimate[1]
     mod.df$lwr[str_which(mod.df$Variable, "Realm")] <- mod.df$lwr[str_which(mod.df$Variable, "Realm")] + mod.df$Estimate[1]
     mod.df$upr[str_which(mod.df$Variable, "Realm")] <- mod.df$upr[str_which(mod.df$Variable, "Realm")] + mod.df$Estimate[1]

     mod.df$Estimate[str_which(mod.df$Variable, "PDF:")] <- mod.df$Estimate[str_which(mod.df$Variable, "PDF:")] + mod.df$Estimate[2]
     mod.df$lwr[str_which(mod.df$Variable, "PDF:")] <- mod.df$lwr[str_which(mod.df$Variable, "PDF:")] + mod.df$Estimate[2]
     mod.df$upr[str_which(mod.df$Variable, "PDF:")] <- mod.df$upr[str_which(mod.df$Variable, "PDF:")] + mod.df$Estimate[2]

     mod.df$Variable <- stringr::str_replace_all(mod.df$Variable, "Realm", "")
    
    if(any(stringr::str_detect(mod.df$Variable, "Afrotropic")) == FALSE) {
      mod.df$Variable[1:2] <- c("Afrotropic", "PDF:Afrotropic")
    } else if(any(stringr::str_detect(mod.df$Variable, "Australasia")) == FALSE) {
      mod.df$Variable[1:2] <- c("Australasia", "PDF:Australasia")
    } else if(any(stringr::str_detect(mod.df$Variable, "Indomalay")) == FALSE) {
      mod.df$Variable[1:2] <- c("Indomalay", "PDF:Indomalay")
    } else if(any(stringr::str_detect(mod.df$Variable, "Nearctic")) == FALSE) {
      mod.df$Variable[1:2] <- c("Nearctic", "PDF:Nearctic")
    } else if(any(stringr::str_detect(mod.df$Variable, "Neotropic")) == FALSE) {
      mod.df$Variable[1:2] <- c("Neotropic", "PDF:Neotropic")
    } else if(any(stringr::str_detect(mod.df$Variable, "Oceania")) == FALSE) {
      mod.df$Variable[1:2] <- c("Oceania", "PDF:Oceania")
    } else if(any(stringr::str_detect(mod.df$Variable, "Palearctic")) == FALSE) {
      mod.df$Variable[1:2] <- c("Palearctic", "PDF:Palearctic")
    }
    n.df <- data.frame("Variable" = names(table(data$Realm)),
                       "n" = as.vector(table(data$Realm)))
    mod.df <- left_join(mod.df, n.df)
  }
  
  if(!("n" %in% colnames(mod.df))) {
    mod.df$n <- nrow(dat.df)
  }
  
  row.names(mod.df) <- NULL
  mod.df <- dplyr::arrange(mod.df, Element, Variable)
  mod.df
}

# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species > 1)

# regression models
mod1.glm <- glmmTMB(MSA_loss ~ 1 + PDF + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)
taxon.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Taxon_group + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)
realm.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Realm + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# regression results
mod1.df <- response.f(mod1.glm, dat.df)
taxon.df <- response.f(taxon.glm, dat.df)
realm.df <- response.f(realm.glm, dat.df)

write.csv(mod1.df, "./Data/04_Model_results_PDF.csv", row.names = FALSE)
write.csv(taxon.df, "./Data/04_Model_results_taxon.csv", row.names = FALSE)
write.csv(realm.df, "./Data/04_Model_results_realm.csv", row.names = FALSE)


# read data
dat.df <- read.csv("./Data/03_PDF_MSA.csv")

# filter data
dat.df <- filter(dat.df, !is.na(Realm), n_species >= 10, Realm != "Oceania")

# regression models
mod1.glm <- glmmTMB(MSA_loss ~ 1 + PDF + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)
taxon.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Taxon_group + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)
realm.glm <- glmmTMB(MSA_loss ~ 1 + PDF * Realm + (1|Source_ID), weights = sqrt(n_species), data = dat.df, family = beta_family(link = "logit"), REML = TRUE, na.action = na.fail)

# regression results
mod1.df <- response.f(mod1.glm, dat.df)
taxon.df <- response.f(taxon.glm, dat.df)
realm.df <- response.f(realm.glm, dat.df)

write.csv(mod1.df, "./Data/04_Model_results_PDF_n10.csv", row.names = FALSE)
write.csv(taxon.df, "./Data/04_Model_results_taxon_n10.csv", row.names = FALSE)
write.csv(realm.df, "./Data/04_Model_results_realm_n10.csv", row.names = FALSE)