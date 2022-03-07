# -------------------------------------------------------------------------#
# Header ----
# -------------------------------------------------------------------------#
# Mo Mär 07 2022
# git hash da55dfa4bb71e0b6cd5a77da0021fbf03026425b
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
library(tidyverse)
library(conveniencefunctions)
library(MRAr)
matrix_2_wide_df <- function(x,nm) {
  as.data.frame(setNames(as.list(x),nm))
}

# .. read results -----
results <- readRDS("results.rds")


.outputFolder <- "Models/Cascade/Noise/plots"
dir.create(.outputFolder)
setwd(.outputFolder)

# .. unnest matrices  -----
results_unnested_matrices <- results %>%
  filter(alpha == 0) %>% 
  filter(!map_lgl(r_0, is.null) && 
           !map_lgl(r_1upstream, is.null) && 
           !map_lgl(r_keptupstream, is.null) && 
           !map_lgl(r_optupstream, is.null) 
         # &&
         # !map_lgl(r_1both, is.null) &&
         # !map_lgl(r_keptboth, is.null) &&
         # !map_lgl(r_optboth, is.null)
  ) %>% 
  mutate(r_optupstream = map(r_optupstream, matrix_2_wide_df, nm = (paste0("upstream_ropt", outer(1:3,1:3, paste0)))),
         # r_optboth = map(r_optboth, . %>% matrix_2_wide_df(paste0("both_ropt", outer(1:3,1:3, paste0)))),
         
         r_keptupstream = map(r_keptupstream, matrix_2_wide_df, nm = (paste0("upstream_rkept", outer(1:3,1:3, paste0)))),
         # r_keptboth = map(r_keptboth, . %>% matrix_2_wide_df(paste0("both_rkept", outer(1:3,1:3, paste0)))),
         
         r_1upstream = map(r_1upstream, matrix_2_wide_df, nm = (paste0("upstream_r1", outer(1:3,1:3, paste0)))),
         # r_1both = map(r_1both, . %>% matrix_2_wide_df(paste0("both_r1", outer(1:3,1:3, paste0)))),
         
         r_0 = map(r_0, . %>% matrix_2_wide_df(paste0("r_0", outer(1:3,1:3, paste0)))),
         
         alpha = as.numeric(alpha)) %>% 
  # unnest(c(r_optupstream,# r_optboth,
  #        r_keptupstream, #r_keptboth,
  #        r_1upstream, #r_1both, 
  #        r_0))
  {.}

# .. Further post processing  -----
d0 <- data.table(results_unnested_matrices)
d0[,`:=`(rowid = 1:.N)]
d <- 
  rbindlist(list(
    r_0            = setnames(unnest_listcol_dt(d0[,list(r_0, rowid)], "r_0")                      , c(paste0("r", outer(1:3,1:3, paste0)), "rowid")),
    r_1upstream    = setnames(unnest_listcol_dt(d0[,list(r_1upstream, rowid)], "r_1upstream")      , c(paste0("r", outer(1:3,1:3, paste0)), "rowid")),
    r_keptupstream = setnames(unnest_listcol_dt(d0[,list(r_keptupstream, rowid)], "r_keptupstream"), c(paste0("r", outer(1:3,1:3, paste0)), "rowid")),
    r_optupstream  = setnames(unnest_listcol_dt(d0[,list(r_optupstream, rowid)], "r_optupstream")  , c(paste0("r", outer(1:3,1:3, paste0)), "rowid"))
  ),
  idcol = "matrixType")
d <- merge(d0[,list(realization, srel, sabs, rowid)], d)
d <- melt(d, id.vars = c("rowid", "realization", "srel", "sabs", "matrixType"), variable.name = "relement", variable.factor = FALSE, value.name = "rvalue")

# -------------------------------------------------------------------------#
# Plot ----
# -------------------------------------------------------------------------#
# .. Plot distributions  -----
srels <- setNames(c(0.001,0.01,0.05,0.1), paste0(LETTERS[1:4], " srel = ",c(0.001,0.01,0.05,0.1) * 100, " %"))
plots <- lapply(seq_along(srels), function(x) {
  
  y = srels[x]
  
  dplot <- copy(d)
  dplot <- dplot[srel == y]
  dplot <- dplot[matrixType %in% c("r_0", "r_optupstream")]
  # dplot[,`:=`(matrixType = ifelse(matrixType == "r_0", "r(0)", "r(a_opt)"))]
  dplot[,`:=`(matrixType = ifelse(matrixType == "r_0", "a=0", "a_opt"))]
  dplot[,`:=`(matrixType = factor(matrixType, c("a=0", "a_opt")))]
  # dplot[,`:=`(matrixType2 = factor(matrixType, c("r(a_opt)", "r(0)")))]
  dplot[,`:=`(elementType = ifelse(relement %in% c("r12", "r13", "r23"), "sequestration", "regulatory"))]
  
  ylim <- NULL
  if (y >= 0.05) ylim = c(-2,2)
  
  cfggplot(dplot, aes(matrixType, rvalue, color = elementType)) + 
    facet_wrap(~relement, scales = "free") + 
    geom_hline(yintercept = 0, lty = "11", size = 0.2) + 
    geom_boxplot(outlier.size = 0.2) +
    scale_y_continuous(n.breaks = 5, limits = egg::symmetric_range) + 
    scale_color_cf() + 
    coord_cartesian(ylim = ylim) + 
    guides(color = FALSE) + 
    labs(x = NULL, y = NULL, title = names(srels)[x]) + 
    theme_cf(7, FLAGbold = FALSE) + 
    theme(title = element_text(face = "bold"))
})

pl <- gridExtra::arrangeGrob(grobs = plots)
# grid::grid.draw(pl)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
ggsave(file.path(.outputFolder, "006-Distributions.pdf"), pl, width = 15.5, height = 16, scale = 1, units = "cm")

# .. plot rkept  -----
srels <- setNames(c(0.001,0.01,0.05,0.1), paste0(LETTERS[1:4], " srel = ",c(0.001,0.01,0.05,0.1) * 100, " %"))
x <- (seq_along(srels))[[1]]
plots <- lapply(seq_along(srels), function(x) {
  
  y = srels[x]
  
  dplot <- copy(d)
  dplot <- dplot[srel == y]
  dplot <- dplot[matrixType %in% c("r_keptupstream")]
  # dplot[,`:=`(matrixType = ifelse(matrixType == "r_0", "r(0)", "r(a_opt)"))]
  dplot[,`:=`(matrixType = ifelse(matrixType == "r_0", "a=0", "a_opt"))]
  dplot[,`:=`(matrixType = factor(matrixType, c("a=0", "a_opt")))]
  # dplot[,`:=`(matrixType2 = factor(matrixType, c("r(a_opt)", "r(0)")))]
  dplot[,`:=`(classification = c(ifelse(rvalue == 1,"sequestration", "regulatory")))]
  dplot <- dplot[,`list`(class = c("seq.", "regul."), 
                         frequency = c(mean(rvalue), 1-mean(rvalue))), by = c("relement")]
  dplot[,`:=`(elementType = ifelse(relement %in% c("r12", "r13", "r23"), "sequestration", "regulatory"))]
  dplot[relement %in% c("r11", "r22", "r33") ,`:=`(frequency  = NA)]
  
  cfggplot(dplot, aes(class, frequency, fill = elementType)) + 
    facet_wrap(~relement, scales = "free") + 
    # geom_hline(yintercept = 0, lty = "11", size = 0.2) + 
    geom_col(color = NA) +
    scale_y_continuous(n.breaks = 3, limits = c(0,1.1)) + 
    scale_color_cf(aesthetics = c("fill", "color")) + 
    guides(color = FALSE, fill = FALSE) + 
    labs(x = NULL, y = NULL, title = names(srels)[x]) + 
    theme_cf(7, FLAGbold = FALSE) + 
    theme(title = element_text(face = "bold"))
})

pl <- gridExtra::arrangeGrob(grobs = plots)
# grid::grid.draw(pl)
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
ggsave(file.path(.outputFolder, "007-Classification.pdf"), pl, width = 15.5, height = 16, scale = 1, units = "cm")



