# -------------------------------------------------------------------------#
# 0 Header ----
# -------------------------------------------------------------------------#
#
# S001-Simulations.R
#
# [PURPOSE]
# 
# Run some simulations to scan the influence of the complex
# "At how much complex is network reconstruction affected?"
#
#
# [AUTHOR]
# Daniel Lill
#
# [Date]
# Wed Feb 23 18:55:16 2022
#
# [Git-hash]
# 6229212d12b394b9e9250be76374556896c688f0
#
rm(list = ls(all.names = TRUE))
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
library(MRAr)
library(scales)
library(xtable)


.outputFolder <- paste0("../04-Output/", "S001-Simulations")
for(folder in c(.outputFolder)) 
if(!dir.exists(folder)) dir.create(folder)

# -------------------------------------------------------------------------#
# 1 Set up model ----
# -------------------------------------------------------------------------#
# .. Prediction function ----- #
modelpath  <- paste0("./")
modelname  <- "MassActionCascade2Modules"
mymodel    <- as.eqnlist(read.csv(paste0(modelpath, modelname, ".csv")))
myodemodel <- odemodel(mymodel, modelname = modelname, deriv = T)

x  <- Xs(odemodel = myodemodel)
xs <- Xs_steady(myodemodel)

# .. Observation function ----- #
modules <- modules0 <- c("X","Y")
obs <- obs0 <- c("m2+a*Cm2e1+a_1*Cm2e1",
                 "e2+a*Cm2e1+a_2*Cm2e1")
g <- g0 <- Y(as.eqnvec(structure(obs0, names = modules0)), mymodel, compile = TRUE, modelname = "obsfn0", attach.input = FALSE)

# .. Parameters ----- #
pars_raw <- read.csv(paste0(modelpath, "pars.csv"), header = FALSE)  # dynamical parameters
pars_raw <- structure(pars_raw$V2, names = as.character(pars_raw$V1))
ic_raw <- read.csv(paste0(modelpath, "IC.csv"), header = F)          # initial conditions
ic_raw <- structure(ic_raw$V2, names = as.character(ic_raw$V1))
# sort ics like they appear in x, so it can be used in stode()
vars_x <- attr(x, "parameters")[attr(x, "parameters") %in% names(ic_raw)]
ic_raw <- ic_raw[vars_x]

pars_inner_opt_0 <- setdiff(attr(g, "parameters"),names(c(pars_raw,ic_raw)))
pars_inner_opt_0 <- structure(rep(0, length(pars_inner_opt_0)), names = pars_inner_opt_0)
pars_inner_0 <- pars_inner <- c(ic_raw, pars_raw, pars_inner_opt_0)

# .. Parameter transformation function ----- #
# Logtrafo of all pars
logtrafo <- structure(paste0("exp( log", names(pars_inner_0), ")"), names = names(pars_inner_0))
# initial values
logtrafo[c("m2", "Cm2e1", "e2")] <- 0

p_log <- P(trafo = logtrafo, compile = TRUE, modelname = "p_log")

pars_opt_0 <- pars_opt <- structure(rep(log(0+.Machine$double.eps), length(pars_inner_opt_0)), names = paste0("log", names(pars_inner_opt_0)))
pars_0     <- pars     <- structure(log(pars_inner_0+.Machine$double.eps), names = paste0("log", names(pars_inner_0))) 
pars_0[names(pars_opt_0)] <- log(0+.Machine$double.eps)

# .. Choose perturbed parameters, set up p_pert ----- #
# Which pars shall be perturbed?
pars_perturbed_0 <- pars_perturbed <- c(logm1 = log(0.9), loge1 = log(0.9))
p_pert <- p_pert_fun(pars_perturbed = pars_perturbed_0, pars = pars_0)

# .. Run first prediction ----- #
# Simulate
perturbation_prediction <- (xs*p_log*p_pert)(times = c(0,Inf), pars = pars, deriv = T)
# Calculate r
MRAr::r_alpha_fun(pars_opt = pars_opt_0*c(1,0,1) + c(0,0,0),
                  perturbation_prediction = perturbation_prediction,
                  obs_fun = g, 
                  p_fun = (p_log*p_pert),
                  pars = pars)


# # -------------------------------------------------------------------------#
# # Run predictions to scan complexes ----
# # -------------------------------------------------------------------------#
# # .. Define parameters ----- #
# scan_pars <- expand.grid(
#   logk3on  = seq(-5,5,1), 
#   logk3off = seq(-5,5,1), 
#   logk4    = seq(-5,5,1)
#   )
# pars_list <- lapply(seq_len(nrow(scan_pars)), function(i) {
#   p <- pars
#   pi <- unlist(scan_pars[i,])
#   p[names(pi)] <- pi
#   p
# })
# 
# # .. Run ODE predictions ----- #
# ncores <- 8
# prediction_list <- parallel::mclapply(X = seq_len(nrow(scan_pars)), mc.cores = ncores, FUN = function(i) {
#   pred_i <- (xs*p_log*p_pert)(times = c(0,Inf), pars = pars_list[[i]], deriv = T)
# })
# 
# # .. Transform predictions ----- #
# ncores <- 8
# prediction_list_long <- parallel::mclapply(X = seq_len(nrow(scan_pars)), mc.cores = ncores, FUN = function(i) {
#   pred_i <- prediction_list[[i]]
#   pred_i <- as.data.frame(pred_i)
#   
#   r <- r_alpha_fun(pars_opt = pars_opt_0*c(1,0,1) + c(0,0,0),
#                    perturbation_prediction = prediction_list[[i]],
#                    obs_fun = g, 
#                    p_fun = (p_log*p_pert),
#                    pars = pars_list[[i]])
#   
#   pred_i <- cbind(pred_i, r21 = r[2], r12 = r[3], parameterSetId = i, scan_pars[i,])
#   pred_i <- data.table(pred_i)
#   pred_i
# })
# predictions <- data.table::rbindlist(prediction_list_long)
# 
# # .. Calculate influence of complex ----- #
# d <- copy(predictions)
# d <- d[condition == "Ctr"]
# d <- d[name == "Cm2e1"]
# d[,`:=`(value = log10(value))]
# d <- d[between(r12, 0, 2)]
# dplot <- melt(d, measure.vars = c("r12", "r21"), variable.name = "rId", variable.factor = FALSE, value.name = "rVal")
# 
# cfggplot(dplot, aes(value, rVal)) +
#   facet_wrap(~rId) + 
#   geom_point(aes(color = logk4), size = 0.1) + 
#   scale_color_viridis_c() + 
#   labs(x = "log10 species")
# 



# -------------------------------------------------------------------------#
# Run predictions to scan complex for coupling k3on and k3off  ----
# -------------------------------------------------------------------------#
scan_pars <- data.frame(
  logk3on = pars["logk3on"]   + seq(-10,5,0.5),
  logk3off = pars["logk3off"] + seq(-10,5,0.5) - 0 ,
  logk4 = pars["logk4"] - 2
  )
i <- (seq_len(nrow(scan_pars)))[[1]]
pars_list <- lapply(seq_len(nrow(scan_pars)), function(i) {
  p <- pars
  pi <- unlist(scan_pars[i,,drop=FALSE])
  p[names(pi)] <- pi
  p
})

ncores <- 8
prediction_list <- parallel::mclapply(X = seq_len(nrow(scan_pars)), mc.cores = ncores, FUN = function(i) {
  pred_i <- (xs*p_log*p_pert)(times = c(0,Inf), pars = pars_list[[i]], deriv = T)
})

ncores <- 8
prediction_list_long <- parallel::mclapply(X = seq_len(nrow(scan_pars)), mc.cores = ncores, FUN = function(i) {
  pred_i <- prediction_list[[i]]
  pred_i <- as.data.frame(pred_i)
  
  r <- r_alpha_fun(pars_opt = pars_opt_0*c(1,0,1) + c(0,0,0),
                   perturbation_prediction = prediction_list[[i]],
                   obs_fun = g, 
                   p_fun = (p_log*p_pert),
                   pars = pars_list[[i]])
  
  pred_i <- cbind(pred_i, r21 = r[2], r12 = r[3], parameterSetId = i, scan_pars[i,,drop=FALSE])
  pred_i <- data.table(pred_i)
  pred_i
})
predictions <- data.table::rbindlist(prediction_list_long)

# .. Calculate influence of complex ----- #
d <- copy(predictions)
d <- d[condition == "Ctr"]
d <- d[name == "Cm2e1"]
d <- d[between(r12, 0, 2)]
dplot <- melt(d, measure.vars = c("r12", "r21"), variable.name = "rId", variable.factor = FALSE, value.name = "rVal")

# Depending on the amount of complex, the connection strength varies
# But this is due to saturation, as shown in the next plot
cfggplot(dplot, aes(value, rVal)) +
  facet_wrap(~rId) +
  geom_point(aes(color = logk3on), size = 0.1) +
  scale_color_viridis_c() +
  scale_x_log10() +
  labs(x = "log10 Complex")

# .. plot all species vs k4 ----- #
d <- copy(predictions)
d <- d[condition == "Ctr"]
# d[,`:=`(value = log10(value))]

# Here it shows that everything is saturated for high logk3on and logk3off values
cfggplot(d, aes(logk3on, value)) +
  facet_wrap(~name, scales = "free") + 
  scale_y_log10(expand = c(0.1,0.1)) + 
  geom_point() +
  geom_blank()

# .. Calculate kms of scan_pars -----
sp <- data.table(scan_pars)
sp[,`:=`(km = (exp(logk3off) + exp(logk4))/ exp(logk3on))]
sp

# -------------------------------------------------------------------------#
# Run predictions to scan complexes by scanning k4 only  ----
# -------------------------------------------------------------------------#
scan_pars <- expand.grid(
  logk4    = seq(-5,5,0.1)
  )
i <- (seq_len(nrow(scan_pars)))[[1]]
pars_list <- lapply(seq_len(nrow(scan_pars)), function(i) {
  p <- pars
  pi <- unlist(scan_pars[i,,drop=FALSE])
  p[names(pi)] <- pi
  p
})

ncores <- 8
prediction_list <- parallel::mclapply(X = seq_len(nrow(scan_pars)), mc.cores = ncores, FUN = function(i) {
  pred_i <- (xs*p_log*p_pert)(times = c(0,Inf), pars = pars_list[[i]], deriv = T)
})

ncores <- 8
prediction_list_long <- parallel::mclapply(X = seq_len(nrow(scan_pars)), mc.cores = ncores, FUN = function(i) {
  pred_i <- prediction_list[[i]]
  pred_i <- as.data.frame(pred_i)
  
  r <- r_alpha_fun(pars_opt = pars_opt_0*c(1,0,1) + c(0,0,0),
                   perturbation_prediction = prediction_list[[i]],
                   obs_fun = g, 
                   p_fun = (p_log*p_pert),
                   pars = pars_list[[i]])
  
  pred_i <- cbind(pred_i, r21 = r[2], r12 = r[3], parameterSetId = i, scan_pars[i,,drop=FALSE])
  pred_i <- data.table(pred_i)
  pred_i
})
predictions <- data.table::rbindlist(prediction_list_long)

# .. Calculate influence of complex ----- #
d <- copy(predictions)
d <- d[condition == "Ctr"]
d <- d[name == "Cm2e1"]
d[,`:=`(value = log10(value))]
d <- d[between(r12, 0, 2)]
dplot <- melt(d, measure.vars = c("r12", "r21"), variable.name = "rId", variable.factor = FALSE, value.name = "rVal")

# Apparently the less complex, the less influence of module 1 on module 2
# But this is misleading, see next subsection
cfggplot(dplot, aes(value, rVal)) +
  facet_wrap(~rId) + 
  geom_point(aes(color = logk4), size = 0.1) + 
  scale_color_viridis_c() + 
  labs(x = "log10 species")

# .. plot all species vs k4 ----- #
d <- copy(predictions)
d <- d[condition == "Ctr"]
d[,`:=`(value = log10(value))]

# Here it shows that the whole e2 is saturated for high logk4 values
cfggplot(d, aes(logk4, value)) +
  facet_wrap(~name, scales = "free") + 
  geom_point() + 
  geom_blank()

# 



# -------------------------------------------------------------------------#
# OLD ----
# -------------------------------------------------------------------------#
# .. Evaluate the Jacobian ----- #
# The value of m2.k3on should be the same as (- m2*e1)
myfunc <- myodemodel$extended %>% attr("eq") %>% funC0()
jacargs1 <- c(perturbation_prediction[[1]][-1],exp(pars)[-(1:5)]) %>% set_names(names(pars_inner_0))
jacargs2 <- perturbation_prediction[[1]] %>% attr("der") %>% drop() %>% .[-1]
do.call(myfunc, as.list(c(jacargs1, jacargs2)))[,"m2.k3on"] == -perturbation_prediction[[1]][3]*perturbation_prediction[[1]][5]

# .. Generate the perturbation data and run the algorithm ----- #
perturbation_prediction <- perturbation_prediction_0 <-  (xs*p_log*p_pert)(times = c(0,Inf), pars = pars, deriv = F)
# perturbation_prediction
r_kept_0 <- r_kept_fun(pars_opt = pars_opt_0, obs_fun = g, p_fun = (p_log*p_pert))
# r_kept_0 <- c(F,F,T,F)

# perturbation_prediction_0 <-  (xs*p_pert)(times = c(0,Inf), pars = pars_0, deriv = F)

r_0 <-   R_fun(pars_opt = pars_opt[c("loga_1")],
               perturbation_prediction = perturbation_prediction,
               obs_fun = g,
               p_fun = (p_log*p_pert),
               pars = pars) %>% local_response_matrix_eq10()
r_0%>% round(3)

myfits <- mstrust(obj_alpha,
                  center =  (pars_opt_0-pars_opt_0)["loga_1"],
                  studyname = "Fits",
                  cores = 3,
                  fits = 3,
                  sd = 3, mypars = pars)
# myfits %>% as.parframe() %>% plotValues()

(best_fit <- myfits %>% as.parframe() %>% as.parvec())
r_alpha_fun(best_fit, pars = pars) %>% round(3)




# .. Compute the local response matrix in dependence of a ----- #
# Use par_opt_setting to define, to which module(s) the complex should be added.
par_opt_settings <- list("loga", "loga_1", "loga_2", c("loga_1", "loga_2")) %>% set_names(1:4)
par_opt_settings <- par_opt_settings[1:2]
dr_list <- lapply(par_opt_settings, function(i) {
  r_opt_dose_response(which_par = "logk1", 
                      dosages = pars_0["logk1"]+seq(-2,2,by = 2),
                      pars_opt = pars_opt_0[i],
                      alpha_scan_range = sort(unique(seq(-10,10,by = .1))),
                      pars = pars
  )
})




# .. Plot the local response matrix in dependence of a ----- #
# - Different opacities correspond to different communicating species
myplot <- dr_list %>%   
  # extract2(1) %>%
  combine_dr_list() %>% 
  filter(par_opt_setting %in% 1:2) %>%
  filter(near(dose,log(0.1))) %>%
  filter(! (r_element %in% paste0("r", c(11,22,33)))) %>% 
  # filter(alpha >-5, alpha <3) %>% 
  filter(matrix == "r_alpha") %>%
  apply_expression({
    x$matrix[x$matrix == "r_alpha"] <- "r"
    x$matrix[x$matrix == "r_opt"] <- "r(a = a_opt)"
    x$matrix <- x$matrix %>% factor(levels = c("r", "r(a = a_opt)", "r_kept"))
    x
  }) %>% 
  mutate(Element = r_element) %>% 
  mutate(points = value) %>% 
  apply_expression({
    nears <- seq(-10,10,by=0.1)[apply(sapply(c(-8.5,-0.7,2.1,9), function(i) near(seq(-10,10,by=0.1),i)),1,any)]
    x[!(x$alpha %in% nears),"points"] <- NA
    x[(x$alpha < (-6))|(x$alpha > (7)),"value"] <- NA
    filter(x, alpha >= (-9))
  }) %>% 
  apply_expression({
    x$par_opt_setting[x$par_opt_setting == 1] <- "X* and Y*"
    x$par_opt_setting[x$par_opt_setting == 2] <- "X*"
    x$par_opt_setting <- x$par_opt_setting %>% factor(levels = c("X* and Y*","X*"))
    x
  }) %>%
  
  ggplot(mapping = aes(x = alpha, y = value, color = Element, linetype = Element, alpha = par_opt_setting)) +
  geom_line() +
  theme_dMod() +
  scale_color_dMod() +
  geom_abline(slope = 0,intercept = 0, linetype = 3) +
  # geom_rect(xmin = -0.7, xmax = 3.2, ymin = -10, ymax = 10, alpha = 0.01, size = 0)+
  geom_point(aes(y = points), size = 3)+
  # facet_grid(.~par_opt_setting, scales = "free")+
  # facet_grid(matrix~par_opt_setting, scales = "free")+
  scale_alpha_discrete(name = "Complex \nadded to", range = c(0.3,1)) +
  scale_x_continuous(breaks = c(-5,0,5))+
  xlab("log(a)")
# plotly::ggplotly(myplot)
myplot
# ggsave(plot = myplot,
#        filename = paste0("small width algorithm pos12",".png"),
#        width = 10,
#        height = 8,
#        units = "cm",
#        device = "png")



# Exit ----
future::plan("sequential")


