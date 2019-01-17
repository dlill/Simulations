library(magick)
setwd("~/Promotion/Projects/MRA/Simulations/")
setwd("Models/Prabakaran/")

print_matrices_at_point_locations <- function(point_locations, plot_data, path) {
  map(point_locations, function(i) {
    plot_data  %>%
      filter(near(alpha, i, tol = 1e-6)) %>%
      .$value %>% c(-1, ., -1) %>%  round(2) %>% 
      paste0("\\mathbf{",., "}") %>% 
      paste0(c("", "\\textcolor{red}{", "\\textcolor{blue}{", ""),.,c("", "}", "}", "")) %>% matrix(ncol = 2) %>% 
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
  
  system('pdflatex -synctex=1 -interaction=nonstopmode "matrices_of_plot2upper_Prabakaran".tex')
  # system2("pandoc", c(paste0("-o ", path, ".pdf"), paste0(path, ".tex")))
}

# ---------------------------------------------------------- #
# plots2upper ----
# ---------------------------------------------------------- #

# ----------------------------------------------- #
# .. Plot ----
# ----------------------------------------------- #
plot_data <-
  dr_list %>%   
  combine_dr_list() %>% 
  filter(!is.na(value)) %>%
  remove_diags %>% 
  
  filter(matrix == "r_alpha") %>%
  
  mutate(alpha = exp(alpha)) %>%
  filter(alpha<5) %>% #.$alpha
  
  fancify_labels %>%
  filter(par_opt_setting %in% "X*") %>%
  
  assign_point_locations %>% 
  {.}

geom_rect_data <- plot_data %>% 
  filter(alpha %in% point_locations) %>% 
  group_by(alpha) %>% 
  summarise(maxval = max(value), minval = min(value)) %>% 
  ungroup

text_data <- geom_rect_data %>% 
  mutate(label=rank(alpha)) %>% 
  # group_by(alpha) %>% 
  mutate(value = c(minval-0.25, maxval+0.25)[c(4,5,3)]) #%>% 
# ungroup

plot_data <- plot_data %>% left_join(geom_rect_data)


myplot <- plot_data %>%   
  ggplot(mapping = aes(x = alpha, y = value, color = Element
                       # uncomment for one par_opt_setting only
                       # , alpha = par_opt_setting
  )) +
  
  
  geom_rect(aes(xmin = alpha - 0.17, xmax = alpha + 0.17, ymin = minval - 0.13, ymax = maxval + 0.13), alpha = 1, fill = "#E9E9E9", color = "black", size = 0.3) + 
  geom_text(aes(label = label), data = text_data, color = "black", size = 5, fontface = "bold") +
  
  geom_line() +
  theme_dMod() +
  scale_colour_manual(values = c("blue", "red", "green")) +
  geom_abline(slope = 0,intercept = 0, linetype = 3) +
  geom_point(aes(y = points), size = 3)+
  scale_alpha_discrete(name = "Complex \nadded to", range = c(0.3,1)) +
  xlab("       Weight parameter a") +
  ylab("Connection coefficient") +
  scale_x_continuous(breaks = sort(c(0:5, point_locations[2])), labels = c("0", expression(a[opt]), 1:5))

# ----------------------------------------------- #
# .. matrices ----
# ----------------------------------------------- #

print_matrices_at_point_locations(point_locations, plot_data, path = "matrices_of_plot2upper_Prabakaran")

# ----------------------------------------------- #
# .. imagemagick ----
# ----------------------------------------------- #
mat1 <- image_read_pdf("matrices_of_plot2upper_Prabakaran.pdf") %>%
  image_crop("284x100+590+524")
mat2 <- image_read_pdf("matrices_of_plot2upper_Prabakaran.pdf") %>%
  image_crop("234x100+590+626")
mat3 <- image_read_pdf("matrices_of_plot2upper_Prabakaran.pdf") %>%
  image_crop("284x100+590+728")
# ----
perc1 <- 80
perc2 <- 48

X <- (c(x1 = 100, x2 = 10, x3 = 166, x4 = 185, x5 = 725, x6 = 700) + 10) * perc1/100
DX = c(dx1 = 45, dx2.1 = 155, dx2 = 125, dx4 = -10, dx5 = 0, dx6 = -10) * perc1/100
Y <- c(y1 = 685, y2 = 752, y3 = 564, y4 = 490, y5 = 727) * perc1/100  - 12

fig <- image_graph(1000 * perc1 / 100, 800 * perc1/100, res = 96 /3 * 5)
myplot
fig <- image_draw(fig)

lines(seq(X[1], X[2], length.out = 100), seq(Y[1],Y[2], length.out = 100), lty = 2)
lines(seq(X[1]+DX[1], X[2] + DX[2], length.out = 100), seq(Y[1],Y[2], length.out = 100), lty = 2)

lines(seq(X[3], X[4], length.out = 100), seq(Y[3],Y[5], length.out = 100), lty = 2)
lines(seq(X[3] + DX[1], X[4] + DX[3], length.out = 100), seq(Y[3],Y[2], length.out = 100), lty = 2)

lines(seq(X[5], X[6], length.out = 100), seq(Y[4],Y[2], length.out = 100), lty = 2)
lines(seq(X[5] + DX[1], X[6] + DX[3], length.out = 100), seq(Y[4],Y[2], length.out = 100), lty = 2)
dev.off()
fig %>% 
  image_composite(image_scale(mat1, geometry_size_percent(perc2)), offset = paste0("+",X[2]+DX[4],"+",Y[2])) %>% 
  image_composite(image_scale(mat2, geometry_size_percent(perc2)), offset = paste0("+",X[4]+DX[5], "+",Y[2])) %>% 
  image_composite(image_scale(mat3, geometry_size_percent(perc2)), offset = paste0("+",X[6] + DX[6],"+",Y[2])) %>% 
image_write("algotoupper.png")

# ---------------------------------------------------------- #
# plots2upper ----
# ---------------------------------------------------------- #

# ----------------------------------------------- #
# .. Plot ----
# ----------------------------------------------- #

plot_data <-
  dr_list %>%   
  combine_dr_list() %>% 
  filter(!is.na(value)) %>%
  remove_diags %>% 
  
  filter(matrix == "r_alpha") %>%
  
  mutate(alpha = exp(alpha)) %>%
  filter(alpha<5) %>% #.$alpha
  
  fancify_labels %>%
  filter(par_opt_setting %in% "X* and Y*") %>%
  
  assign_point_locations %>% 
  {.}

geom_rect_data <- plot_data %>% 
  filter(alpha %in% point_locations) %>% 
  group_by(alpha) %>% 
  summarise(maxval = max(value), minval = min(value))


text_data <- geom_rect_data %>% 
  mutate(label=rank(alpha)) %>% 
  # group_by(alpha) %>% 
  mutate(value = c(minval-0.25, maxval+0.25)[c(4,2,3)]) #%>% 
# ungroup

plot_data <- plot_data %>% left_join(geom_rect_data)

# print_matrices_at_point_locations(point_locations, plot_data)

myplot <- plot_data %>%   
  ggplot(mapping = aes(x = alpha, y = value, color = Element
                       # uncomment for one par_opt_setting only
                       # , alpha = par_opt_setting
  )) +
  
  geom_rect(aes(xmin = alpha - 0.17, xmax = alpha + 0.17, ymin = minval - 0.13, ymax = maxval + 0.13), alpha = 1, fill = "#E9E9E9", color = "black", size = 0.3) + 
  geom_text(aes(label = label), data = text_data, color = "black", size = 5, fontface = "bold") +
  
  geom_line() +
  theme_dMod() +
  scale_colour_manual(values = c("blue", "red", "green")) +
  geom_abline(slope = 0,intercept = 0, linetype = 3) +
  geom_point(aes(y = points), size = 3)+
  scale_alpha_discrete(name = "Complex \nadded to", range = c(0.3,1)) +
  xlab("      Weight parameter a") +
  ylab("Connection coefficient") +
  scale_x_continuous(breaks = sort(c(0:5, point_locations[2])), labels = c("0", expression(a[opt]), 1:5))


myplot

# ----------------------------------------------- #
# .. matrices ----
# ----------------------------------------------- #


print_matrices_at_point_locations(point_locations, plot_data, path = "matrices_of_plot2both")

# ----------------------------------------------- #
# .. imagemagick ----
# ----------------------------------------------- #
mat1 <- image_read_pdf("matrices_of_plot2both.pdf") %>%
  image_crop("284x100+590+524")
mat2 <- image_read_pdf("matrices_of_plot2both.pdf") %>%
  image_crop("234x100+590+626")
mat3 <- image_read_pdf("matrices_of_plot2both.pdf") %>%
  image_crop("284x100+590+728")
# ----
perc1 <- 80
perc2 <- 48

X <- (c(x1 = 100, x2 = 10, x3 = 166, x4 = 185, x5 = 725, x6 = 700) + 10) * perc1/100
DX = c(dx1 = 45, dx2.1 = 155, dx2 = 125, dx4 = -10, dx5 = 0, dx6 = -5) * perc1/100
Y <- c(y1 = 685, y2 = 752, y3 = 530, y4 = 400, y5 = 727) * perc1/100  - 12

fig <- image_graph(1000 * perc1 / 100, 800 * perc1/100, res = 96 /3 * 5)
myplot
fig <- image_draw(fig)

lines(seq(X[1], X[2], length.out = 100), seq(Y[1],Y[2], length.out = 100), lty = 2)
lines(seq(X[1]+DX[1], X[2] + DX[2], length.out = 100), seq(Y[1],Y[2], length.out = 100), lty = 2)

lines(seq(X[3], X[4], length.out = 100), seq(Y[3],Y[5], length.out = 100), lty = 2)
lines(seq(X[3] + DX[1], X[4] + DX[3], length.out = 100), seq(Y[3],Y[2], length.out = 100), lty = 2)

lines(seq(X[5], X[6], length.out = 100), seq(Y[4],Y[2], length.out = 100), lty = 2)
lines(seq(X[5] + DX[1], X[6] + DX[3], length.out = 100), seq(Y[4],Y[2], length.out = 100), lty = 2)
dev.off()
fig %>% 
  image_composite(image_scale(mat1, geometry_size_percent(perc2)), offset = paste0("+",X[2]+DX[4],"+",Y[2])) %>% 
  image_composite(image_scale(mat2, geometry_size_percent(perc2)), offset = paste0("+",X[4]+DX[5], "+",Y[2])) %>% 
  image_composite(image_scale(mat3, geometry_size_percent(perc2)), offset = paste0("+",X[6] + DX[6],"+",Y[2])) %>% 
  image_write("algotoboth.png")

# ---------------------------------------------------------- #
# old ----
# ---------------------------------------------------------- #
perc1 <- 80
perc2 <- 40

X <- (c(x1 = 100, x2 = 20, x3 = 166, x4 = 160, x5 = 725, x6 = 700) + 10) * perc1/100
DX = c(dx1 = 45, dx2 = 125, dx3 = -15) * perc1/100
Y <- c(y1 = 685, y2 = 752, y3 = 524, y4 = 390) * perc1/100  - 12

fig <- image_graph(1000 * perc1 / 100, 800 * perc1/100, res = 96 /3 * 5)
myplot
fig <- image_draw(fig)

lines(seq(X[1], X[2], length.out = 100), seq(Y[1],Y[2], length.out = 100), lty = 2)
lines(seq(X[1]+DX[1], X[2] + DX[2], length.out = 100), seq(Y[1],Y[2], length.out = 100), lty = 2)

lines(seq(X[3], X[4], length.out = 100), seq(Y[3],Y[2], length.out = 100), lty = 2)
lines(seq(X[3] + DX[1], X[4] + DX[3], length.out = 100), seq(Y[3],Y[2], length.out = 100), lty = 2)

lines(seq(X[5], X[6], length.out = 100), seq(Y[4],Y[2], length.out = 100), lty = 2)
lines(seq(X[5] + DX[1], X[6] + DX[3], length.out = 100), seq(Y[4],Y[2], length.out = 100), lty = 2)
dev.off()
image_composite(fig, image_scale(mat1, geometry_size_percent(perc2)), offset = paste0("+",X[2]+DX[4],"+",Y[2])) %>% 
  image_composite(image_scale(mat2, geometry_size_percent(perc2)), offset = paste0("+",X[4], "+",Y[2])) %>% 
  image_composite(image_scale(mat3, geometry_size_percent(perc2)), offset = paste0("+",X[6],"+",Y[2])) %>% 
image_write("algotoboth.png")

# ---------------------------------------------------------- #
# Documentation ----
# ---------------------------------------------------------- #

writeLines("algotoboth.png was created by Script S5 in the main folder.", "algotoboth.txt")
writeLines("algotoupper.png was created by Script S5 in the main folder.", "algotoupper.txt")

file.copy("algotoupper.png", "~/Promotion/Writing/Papers/2017 02 MRA optimization procedure Paper/Revision1/Figure2-1.png")
file.copy("algotoboth.png", "~/Promotion/Writing/Papers/2017 02 MRA optimization procedure Paper/Revision1/Figure2-2.png")
file.copy("algotoupper.txt", "~/Promotion/Writing/Papers/2017 02 MRA optimization procedure Paper/Revision1/Figure2-1.txt")
file.copy("algotoboth.txt", "~/Promotion/Writing/Papers/2017 02 MRA optimization procedure Paper/Revision1/Figure2-2.txt")

unlink_dMod(c("aux", "log", "gz"))
