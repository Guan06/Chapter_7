#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
source("../01.common_scripts/plot_setting.R")

density_file <- "../00.data/Global_density.txt"
density <- read.table(density_file, header = T, sep = "\t")

c_method <- c("Pearson" = "#56B4E9",
              "Spearman" = "#009E73",
              "SparCC" = "#E69F00")

## for panel A --> bacteria
den_bac <- density[density$Kingdom == "Bacteria", ]

pa_1 <- ggplot(den_bac, aes(Network, Density, fill = Network)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    main_theme +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(y = "Network density", x = "") +
    scale_x_discrete(limits = c("Pearson", "Spearman", "SparCC")) +
    scale_fill_manual(values = c_method) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_line(color = NA),
          legend.position = "top")

d1 <- read.table("Bacteria_Pearson_node.txt", header = T, sep = "\t")
d2 <- read.table("Bacteria_Spearman_node.txt", header = T, sep = "\t")
d3 <- read.table("Bacteria_SparCC_node.txt", header = T, sep = "\t")

d1 <- d1[, c(1, 2)]
d1$Network <- "Pearson"

d2 <- d2[, c(1, 2)]
d2$Network <- "Spearman"

d3 <- d3[, c(1, 2)]
d3$Network <- "SparCC"

degree <- rbind(d1, d2, d3)
pa_2 <- ggplot(degree, aes(Degree, fill = Network, color = Network)) +
    geom_density(alpha = 0.7) +
    main_theme +
    labs(y = "", x = "Degree distribution") +
    scale_fill_manual(values = c_method) +
    scale_color_manual(values = c_method) +
    theme(legend.position = "none")

p_legend <- get_legend(pa_1)
pa_1 <- pa_1 + theme(legend.position = "none")

pa_12 <- plot_grid(pa_1, pa_2, nrow = 1, rel_widths = c(0.4, 1),
                   labels = c('a', ''))

## for panel B --> fungi
den_fun <- density[density$Kingdom == "Fungi", ]
pb_1 <- ggplot(den_fun, aes(Network, Density, fill = Network)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    main_theme +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(y = "Network density", x = "") +
    scale_x_discrete(limits = c("Pearson", "Spearman", "SparCC")) +
    scale_fill_manual(values = c_method) +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_line(color = NA),
          legend.position = "none")

d1 <- read.table("Fungi_Pearson_node.txt", header = T, sep = "\t")
d2 <- read.table("Fungi_Spearman_node.txt", header = T, sep = "\t")
d3 <- read.table("Fungi_SparCC_node.txt", header = T, sep = "\t")

d1 <- d1[, c(1, 2)]
d1$Network <- "Pearson"

d2 <- d2[, c(1, 2)]
d2$Network <- "Spearman"

d3 <- d3[, c(1, 2)]
d3$Network <- "SparCC"

degree <- rbind(d1, d2, d3)
pb_2 <- ggplot(degree, aes(Degree, fill = Network, color = Network)) +
    geom_density(alpha = 0.7) +
    main_theme +
    labs(y = "", x = "Degree distribution") +
    scale_fill_manual(values = c_method) +
    scale_color_manual(values = c_method) +
    theme(legend.position = "none")

pb_12 <- plot_grid(pb_1, pb_2, nrow = 1,
                   labels = c('b', ''),
                   rel_widths = c(0.4, 1))

p <- plot_grid(p_legend, pa_12, pb_12, ncol = 1,
               rel_heights = c(0.1, 1, 1))
ggsave("Figure_7_5.pdf", p, height = 3, width = 6)
