#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(parallelDist)
source("../01.common_scripts/plot_setting.R")
source("../01.common_scripts/pcoa_function.R")
source("../01.common_scripts/order.R")

## for figure panel A -> finally changed to panel C
bc_file <- "../00.data/Bac_ASV_bc.rds"
design_file <- "../00.data/Bac_design_3809.txt"

bc <- readRDS(bc_file)
design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)

p_a <- pcoa(dmr, design, 12, "Compartment", "Host.Species")

## for figure panel B -> finally changed to panel D
bc_file <- "../00.data/Fun_ASV_bc.rds"
design_file <- "../00.data/Fun_design_2232.txt"

bc <- readRDS(bc_file)
design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_b <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_b <- p_b + theme(legend.position = "none")

## for panel C -> finally changed to panel A
alpha_file <- "../00.data/Bac_alpha_mean.txt"
design_file <- "../00.data/Bac_design_3809.txt"

alpha <- read.table(alpha_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana", "chlamydomonas_reinhardtii")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)

inter <- intersect(alpha$Sample_ID, design$Original.SampleID)
alpha <- alpha[alpha$Sample_ID %in% inter, ]
design <- design[design$Original.SampleID %in% inter, ]

design <- design[match(alpha$Sample_ID, design$Original.SampleID), ]
design$Com_ST <- paste0(design$Compartment, "_", design$Soil)

alpha$Sample_ID <- design$Sample_ID
alpha_design <- merge(alpha, design)

alpha_design$Shannon <- as.numeric(alpha_design$Shannon)

p_c <- ggplot(alpha_design, aes(x = Group, y = Shannon, color = Compartment,
                       fill = Compartment, shape = Host.Species)) +
    geom_boxplot(outlier.size = 0.3, outlier.shape = NULL) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_x_discrete(limits = plot_order, labels = bac_labels) +
    main_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

## for figure panel D -> finally changed to panel B
alpha_file <- "../00.data/Fun_alpha_mean.txt"
design_file <- "../00.data/Fun_design_2232.txt"
alpha <- read.table(alpha_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)

design$Com_ST <- paste0(design$Compartment, "_", design$Soil)
alpha_design <- merge(alpha, design)

alpha_design$Shannon <- as.numeric(alpha_design$Shannon)

p_d <- ggplot(alpha_design, aes(x = Group, y = Shannon, color = Compartment,
                       fill = Compartment, shape = Host.Species)) +
    geom_boxplot(outlier.size = 0.3, outlier.shape = NULL) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_x_discrete(limits = fun_plot_order, labels = fun_labels) +
    main_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

## organize the panels with cowplot
library(cowplot)
p_a_legend <- get_legend(p_a)
p_a <- p_a + theme(legend.position = "none")
lower <- plot_grid(p_a, p_b, p_a_legend, nrow = 1,
                   labels = c('c', 'd', ''),
                   rel_widths = c(1, 1, 0.3),
                   rel_heights = c(1, 1, 1))

upper <- plot_grid(p_c, p_d, nrow = 1, labels = c('a', "b"),
            rel_widths = c(1.25, 1))

p_aligned <- plot_grid(upper, lower, nrow = 2,
                       rel_heights = c(1.25, 2),
                        align = "v", axis = "l")
ggsave("Figure_7_1.pdf", p_aligned, width = 10)
