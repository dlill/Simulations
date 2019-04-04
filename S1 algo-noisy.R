algo_noisy <- function(pars, # parameters to simulate the model
         which_pars_perturbed, # names of parameters which are perturbed
         which_alpha_pars, # names of alpha parameters
         alpha_pars, # named alpha parameters ?
         alpha, # the alpha such that ralpha is compared against r0
         
         g,  # the functions
         xs,
         p_log,
         
         srel = 0.1, # parameters for error parameters sig = srel * y + sabs
         sabs = 0.1, 
         N = 1000,   # N simulated experiments
         nreplicates = 1, # nreplicates in each experiment
         
         perturbation_effect = 0.9 # fold change of the inner which_pars_perturbed
         ) {
  if (any(!(which_alpha_pars %in% names(alpha_pars))))
    stop("some alpha_pars in alpha_pars_settings which don't exist in alpha_pars")
  
  myperturbation <- rep(log(perturbation_effect), length(which_pars_perturbed)) %>% `names<-`(which_pars_perturbed)
  p_pert <- p_pert_fun(myperturbation, pars = pars, modelname = NULL)
  perturbation_prediction <- (xs*p_log*p_pert)(c(0,Inf), pars, deriv = F)
  
  out <- map(seq_len(N), function(i) {
    set.seed(i)
    perturbation_prediction <- map(perturbation_prediction, function(.x) .x + rnorm(length(.x), 0, sabs + srel * .x))
    out.i <- tibble(r_0 = list(NULL),
                  r_1 = list(NULL),
                  r_kept = list(NULL),
                  parframes = list(NULL),
                  r_opt = list(NULL),
                  perturbation_predictioni = list(perturbation_prediction))
    
    
    r_0 <- R_fun(pars_opt = alpha_pars[which_alpha_pars] * 0 + log(.Machine$double.eps),
                 perturbation_prediction = perturbation_prediction,
                 obs_fun = g,
                 p_fun = (p_log*p_pert),
                 pars = pars) %>% local_response_matrix_eq10()
    out.i$r_0 <- list(r_0 %>% round(2))
    
    r_1 <- R_fun(pars_opt = alpha_pars[which_alpha_pars] * 0 + alpha,
                 perturbation_prediction = perturbation_prediction,
                 obs_fun = g,
                 p_fun = (p_log*p_pert),
                 pars = pars) %>% local_response_matrix_eq10()
    out.i$r_1 <- list(r_1 %>% round(2))
    
    r_kept <- r_kept_fun2(r_0, r_1)
    out.i$r_kept <- list(r_kept)
    
    myfits <- mstrust(obj_alpha,
                      center =  structure(rep(0, length(which_alpha_pars)), names = which_alpha_pars),
                      studyname = "Fits",
                      cores = 1,
                      fits = 3,
                      sd = 1,
                      mypars = pars,
                      perturbation_prediction = perturbation_prediction,
                      r_kept = r_kept,
                      p_fun = (p_log * p_pert),
                      obs_fun = g)
    myparframe <- try(myfits %>% as.parframe())
    out.i$parframes <- list(myparframe)
    if (inherits(myparframe, "try-error")) {
      out.i$parframes <- list(myfits)
      return(out.i)
    }
    
    r_opt <- r_alpha_fun(pars_opt = myfits %>% as.parframe() %>% as.parvec() %>% unclass() ,
                         pars = pars,
                         perturbation_prediction = perturbation_prediction,
                         p_fun = (p_log*p_pert),
                         obs_fun = g)
    
    out.i$r_opt <- list(r_opt %>% round(2))
    return(out.i)
  })
  out <- do.call(rbind, out)
  # out <- bind_cols(out, perturbation_prediction0 = map(seq_len(N),~perturbation_prediction))
  return(out)
}


# perturbation_prediction0 <- perturbation_prediction
# map(perturbation_prediction[-1], ~.x/perturbation_prediction[[1]] - 1) %>% do.call(rbind,.) %>% abs %>% apply(2,max)
# map(perturbation_prediction[-1], ~.x/perturbation_prediction[[1]] - 1) %>% do.call(rbind,.) %>% abs %>% apply(2,min)
# changes <- map(perturbation_prediction[-1], ~.x/perturbation_prediction[[1]] - 1) %>% do.call(rbind,.) %>% abs
# changes[,c(3,6,9)] + cbind(changes[,c(5,8)],0)

test <- algo_noisy(pars, c("x1", "y1", "z1"), c("a_tox_Cxy", "a_toy_Cyz"), alpha_pars, 0, g, xs, p_log, 0.01,0.0, 100, perturbation_effect = 0.2)

map(test$r_opt, matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>% 
  do.call(rbind,.) %>% gather(rij, value, everything()) %>% 
  ggplot(aes(x = rij, y = value)) + 
  geom_violin(draw_quantiles = 0.5) + 
  facet_wrap(~rij, scales = "free")+
  # coord_cartesian(ylim = c(-2,2))+
  ggtitle("ropt") + 
  geom_blank() 

map(test$r_0, matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>% 
  do.call(rbind,.) %>% gather(rij, value, everything()) %>% 
  ggplot(aes(x = rij, y = value)) + 
  geom_violin(draw_quantiles = 0.5) + 
  facet_wrap(~rij, scales = "free")+
  # coord_cartesian(ylim = c(-2,2))+
  ggtitle("r0") + 
  geom_blank() 

map(test$r_1, matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>% 
  do.call(rbind,.) %>% gather(rij, value, everything()) %>% 
  ggplot(aes(x = rij, y = value)) + 
  geom_violin(draw_quantiles = 0.5) + 
  facet_wrap(~rij, scales = "free")+
  coord_cartesian(ylim = c(-2,2))+
  ggtitle("r1") + 
  geom_blank() 

map(test$r_kept, matrix_2_wide_df, nm = (paste0("r", outer(1:3,1:3, paste0)))) %>% 
  do.call(rbind,.) %>% gather(rij, value, everything()) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(aes(y = ..prop..), stat = "count") + 
  facet_wrap(~rij, scales = "free")+
  # coord_cartesian(ylim = c(-2,2))+
  geom_blank() 

