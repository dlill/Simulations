try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
library(tidyverse)
library(conveniencefunctions)
library(MRAr)
matrix_2_wide_df <- function(x,nm) {
  as.data.frame(setNames(as.list(x),nm))
}

results <- readRDS("results.rds")


.outputFolder <- "Models/Cascade/Noise/plots"
dir.create(.outputFolder)
setwd(.outputFolder)

# ..  -----
print_distribution_plots <- function(results, a, s, perturbation_strength) {
  
  plotopt <- results %>%
    filter(alpha == a) %>%
    filter(srel == s) %>%
    filter(map_lgl(pars_perturbed, ~near(.x[1], perturbation_strength))) %>%
    .$r_optupstream %>%
    {.[!map_lgl(., is.null)]} %>%
    map(matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>%
    do.call(rbind,.) %>% gather(rij, value, everything()) %>%
    ggplot(aes(x = rij, y = value)) +
    # geom_violin(draw_quantiles = 0.5) +
    geom_boxplot() +
    facet_wrap(~rij, scales = "free")+
    # coord_cartesian(ylim = c(-2,2))+
    # coord_cartesian(ylim = c(-5,5))+
    ggtitle(paste0("r(a=a_opt), srel = ",s, ", parameter pertubation = ", round(exp(perturbation_strength),2), " * p")) +
    geom_blank()
  # ggsave(paste0("r_opt_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, ".png"), plotopt)
  # ggsave(paste0("srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "_r_opt", ".png"), plotopt)
  # ggsave(paste0("r_opt_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "coord2", ".png"), plotopt +
  #          coord_cartesian(ylim = c(-2,2)))
  print(plotopt)
  
  plot0 <- results %>%
    filter(alpha == a) %>%
    filter(srel == s) %>%
    filter(map_lgl(pars_perturbed, ~near(.x[1], perturbation_strength))) %>%
    .$r_0 %>%
    {.[!map_lgl(., is.null)]} %>%
    map(matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>%
    do.call(rbind,.) %>% gather(rij, value, everything()) %>%
    ggplot(aes(x = rij, y = value)) +
    # geom_violin(draw_quantiles = 0.5) +
    geom_boxplot() +
    facet_wrap(~rij, scales = "free")+
    # coord_cartesian(ylim = c(-2,2))+
    # coord_cartesian(ylim = c(-5,5))+
    ggtitle(paste0("r(a=0), srel = ",s, ", parameter pertubation = ", round(exp(perturbation_strength),2), " * p")) +
    geom_blank()
  # ggsave(paste0("r_0_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, ".png"), plot0)
  # ggsave(paste0("srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "_r_0", ".png"), plot0)
  # ggsave(paste0("r_0_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "coord2", ".png"), plot0 +
  #          coord_cartesian(ylim = c(-2,2)))
  print(plot0)
  
  plot1 <- results %>%
    filter(alpha == a) %>%
    filter(srel == s) %>%
    filter(map_lgl(pars_perturbed, ~near(.x[1], perturbation_strength))) %>%
    .$r_1upstream %>%
    {.[!map_lgl(., is.null)]} %>%
    map(matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>%
    do.call(rbind,.) %>% gather(rij, value, everything()) %>%
    ggplot(aes(x = rij, y = value)) +
    # geom_violin(draw_quantiles = 0.5) +
    geom_boxplot() +
    facet_wrap(~rij, scales = "free")+
    # coord_cartesian(ylim = c(-2,2))+
    # coord_cartesian(ylim = c(-5,5))+
    ggtitle(paste0("r(a=",a,"), srel = ",s, ", parameter pertubation = ", round(exp(perturbation_strength),2), " * p")) +
    geom_blank()
  # ggsave(paste0("r_1_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, ".png"), plot1)
  # ggsave(paste0("srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "_r_1", ".png"), plot1)
  # ggsave(paste0("r_1_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "coord2", ".png"), plot1 +
  #          coord_cartesian(ylim = c(-2,2)))
  # 
  # 
  # a <- asp$a[1]
  # s <- asp$s[1]
  # perturbation_strength <- asp$perturbation_strength[1]
  print(plot1)
  
  plotkept <-
    results %>% 
    filter(alpha == a) %>% 
    filter(srel == s) %>% 
    filter(map_lgl(pars_perturbed, ~near(.x[1], perturbation_strength))) %>% 
    .$r_keptupstream %>% 
    {.[!map_lgl(., is.null)]} %>% 
    map(matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>% 
    do.call(rbind,.) %>% gather(rij, value, everything()) %>% 
    ggplot(aes(x = value, group = rij)) + 
    geom_bar(aes(y = ..prop..)) + 
    facet_wrap(~rij)+
    ggtitle(paste0("Matrix element is minimized, srel = ",s, ", parameter pertubation = ", round(exp(perturbation_strength),2), " * p")) +
    geom_blank() 
  # ggsave(paste0("r_kept_srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, ".png"), plotkept)
  # ggsave(paste0("srel_",s, "delta_p_", round(exp(perturbation_strength),2), "_alpha_", a, "_r_kept", ".png"), plotkept)
  
  print(plotkept)
  print(plotopt + coord_cartesian(ylim = c(-2,2)))
  print(plot0 + coord_cartesian(ylim = c(-2,2)))
  print(plot1 + coord_cartesian(ylim = c(-2,2)))
  
  return(invisible(NULL))
}


asp <- expand.grid(a = 0, s = unique(results$srel), perturbation_strength = map_dbl(unique(results$pars_perturbed),1))

walk(seq_len(nrow(asp)), function(.x) {
  pdf(paste0("boxplots_srel_",asp[.x,]$s, "delta_p_", round(exp(asp[.x,]$perturbation_strength),2), "_alpha_", asp[.x,]$a, ".pdf"))
  with(asp[.x,], print_distribution_plots(results, a, s, perturbation_strength))
  dev.off()
})
# ..  -----

