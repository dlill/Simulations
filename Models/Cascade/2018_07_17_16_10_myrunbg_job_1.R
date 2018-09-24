try(library(parallel))
try(library(xtable))
try(library(scales))
try(library(MRAr))
try(library(miniUI))
try(library(shiny))
try(library(conveniencefunctions))
try(library(forcats))
try(library(stringr))
try(library(dplyr))
try(library(purrr))
try(library(readr))
try(library(tidyr))
try(library(tibble))
try(library(ggplot2))
try(library(tidyverse))
try(library(numDeriv))
try(library(dMod))
try(library(cOde))
try(library(deSolve))
try(library(rootSolve))
try(library(stats))
try(library(graphics))
try(library(grDevices))
try(library(utils))
try(library(datasets))
try(library(methods))
try(library(base))
setwd('~/2018_07_17_16_10_myrunbg_job_1_folder')
rm(list = ls())

load('2018_07_17_16_10_myrunbg_job.RData')
files <- list.files(pattern = '.so')
for (f in files) dyn.load(f)
.node <- 1
.runbgOutput <- try({
    cores <- detectFreeCores()
    alphas <- seq(-0.5, 5, 0.5)
    map(alphas, function(alpha) {
        mclapply(pp_list, function(pp) {
            p_pert <- p_pert_fun(pp, pars = pars, modelname = NULL)
            myframe <- dMod.frame("1", g, x, p_log, NULL, NULL, p_pert = list(p_pert), xs = list(xs), pars = list(pars))
            g <- myframe$g[[1]]
            p_pert <- myframe$p_pert[[1]]
            p <- myframe$p[[1]]
            xs <- myframe$xs[[1]]
            pars <- myframe$pars[[1]]
            perturbation_prediction <- (xs * p * p_pert)(c(0, Inf), pars, deriv = F)
            try({
                r0 <- R_fun(pars_opt = pars_opt, perturbation_prediction = perturbation_prediction, obs_fun = g, p_fun = (p * p_pert), pars = pars) %>% local_response_matrix_eq10()
            })
            if (inherits(r0, "try-error")) 
                return(NULL)
            algo <- function(which_pars_opt, prefix = NULL) {
                rkept <- r_kept_fun(pars_opt = pars_opt[which_pars_opt], perturbation_prediction = perturbation_prediction, obs_fun = g, p_fun = (p * p_pert), pars = pars, alpha = -log(.Machine$double.eps) + alpha)
                myfits <- mstrust(obj_alpha, center = structure(rep(0, length(which_pars_opt)), names = which_pars_opt), studyname = "Fits", cores = 1, fits = 10, sd = 1, mypars = pars, perturbation_prediction = perturbation_prediction, r_kept = rkept, p_fun = (p * p_pert), obs_fun = g)
                try(myfits %>% as.parframe()) %>% {
                  e <- .
                  if (inherits(e, "try-error")) 
                    return(NULL)
                }
                ropt <- r_alpha_fun(pars_opt = myfits %>% as.parframe() %>% as.parvec() %>% unclass(), pars = pars, perturbation_prediction = perturbation_prediction, p_fun = (p * p_pert))
                out <- tibble(list(rkept), list(myfits %>% as.parframe), list(ropt %>% round(2)))
                names(out) <- paste0(prefix, c("r_kept", "parframes", "ropt"))
                return(out)
            }
            to_upstream_module <- algo(c("a_tox_Cxy", "a_toy_Cyz", "a_toz_Czx"), "upstream")
            to_both_modules <- algo(c("a_Cxy", "a_Cyz", "a_Czx"), "both")
            bind_cols(tibble(r0 = list(r0 %>% round(2)), pars_perturbed = list(pp), alpha = alpha), to_upstream_module, to_both_modules)
        }, mc.cores = cores) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
})
save(.runbgOutput, file = '2018_07_17_16_10_myrunbg_job_1_result.RData')