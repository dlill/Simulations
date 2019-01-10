library(conveniencefunctions)

dir_final <- "~/Promotion/Writing/Papers/2017 02 MRA optimization procedure Paper/Code/"
dir_models <- "~/Promotion/Projects/MRA/Simulations/Models/"
setwd(dir_models)

file_index <- c(Cascade.Rmd = "Cascade/Cascade.Rmd",
  `Cascade with cxz complex2.Rmd` = "Cascade/Cascade with cxz different perturbations.Rmd",
  `Cascade with zy feedback.Rmd` = "Cascade/Cascade with gainlf different perturbations.Rmd",
  `Cascade with zx feedback.Rmd` = "Cascade/Cascade with gainuf different perturbations.Rmd",
  Hub.Rmd = "Hub/Hub.Rmd",
  `Hub in cascade.Rmd` = "HubInCascade/HubInCascade.Rmd",
  Phosphatase.Rmd = "Phosphatase/Phosphatase.Rmd",
  Prabakaran.Rmd = "Prabakaran/Prabakaran.Rmd",
  `Prabakaran_model_com_spec.nb` = "Prabakaran/Mathematica/Prabakaran_model_com_spec.nb")
iwalk(file_index, ~file.copy(.x, paste0(dir_final, "/Scripts/",.y)))
