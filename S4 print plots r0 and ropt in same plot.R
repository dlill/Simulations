print_distribution_plots2 <- function(results_unnested_matrices, s, perturbation_strength) {

p1 <-
  results_unnested_matrices %>% 
  filter(srel == s) %>%
  filter(map_lgl(pars_perturbed, ~near(.x[1], perturbation_strength))) %>%
  filter(par_opt_setting == "upstream"|is.na(par_opt_setting)) %>%
  filter(!mat %in% c("r_kept","r(a=1)")) %>%
  mutate(mat = factor(mat))%>% 
  ggplot(aes(x = mat, y = value)) +
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot() +
  facet_wrap(~element, scales = "free")+
  xlab("Local response Matrix") + 
  ylab("Local response coefficient") +
  theme_dMod()+
  ggtitle(paste0("Noise level = ",s*100, "%")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_blank()
  
  ggsave(paste0("boxcombined_srel-", s, "_deltap-", round(exp(perturbation_strength),2), ".png"), p1)
  ggsave(paste0("boxcombined_srel-", s, "_deltap-", round(exp(perturbation_strength),2), "zoom.png"), p1 +
           coord_cartesian(ylim = c(-2,2)))
}

asp <- expand.grid(a = 0, s = map_dbl(error_pars, 2), perturbation_strength = map_dbl(perturbations,1))

# dir.create("~/Promotion/Projects/MRA/Simulations/Models/Cascade/Noise/plots")
setwd("~/Promotion/Projects/MRA/Simulations/Models/Cascade/Noise/plots")

walk(seq_len(nrow(asp)), function(.x) {
  with(asp[.x,], print_distribution_plots2(results_unnested_matrices, s, perturbation_strength))
})





