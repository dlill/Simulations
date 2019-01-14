print_matrices_at_point_locations <- function(point_locations, plot_data, path) {
  map(point_locations, function(i) {
    plot_data  %>%
      filter(near(alpha, i, tol = 1e-6)) %>%
      .$value %>% c(-1, ., -1) %>%  round(2) %>% paste0(c("", "\\textcolor{red}{", "\\textcolor{blue}{", ""),.,c("", "}", "}", "")) %>% matrix(ncol = 2) %>% 
      xtable() %>% print(tabular.enivronment = "pmatrix", floating = FALSE, hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, sanitize.text.function = function(x){x}) %>%
      str_replace_all(c("tabular"= "pmatrix", "\\{rr\\}" = "", "\\\\begin" = "$\\\\begin", "\\{ll\\}" = "")) %>% paste0(., "$") 
  }) %>% 
    paste0(collapse = "\n\n") %>% 
    paste0("\\documentclass[10pt,a4paper]{article}
\\usepackage[latin1]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{graphicx}
\\usepackage{color}
\\begin{document}",., "\n\n\\end{document}") %>% 
    write_lines(paste0(path, ".tex"))
  
  system2("pandoc", c(paste0("-o ", path, ".pdf"), paste0(path, ".tex")))
}

print_matrices_at_point_locations(point_locations, plot_data, path = "matrices_of_plot2upper_Prabakaran")
