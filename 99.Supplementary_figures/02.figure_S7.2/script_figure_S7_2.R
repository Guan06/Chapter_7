#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(ggridges)

source("../../01.common_scripts/plot_setting.R")
source("../../01.common_scripts/order.R")

## for bacteria
asv_file <- "../../00.data/Bac_ASV_raref.rds"
rep_file <- "../../00.data/Bac_ASV_rep.rds"
design_file <- "../../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")

asv <- asv[rownames(asv) %in% rownames(rep), ]

RA0 <- data.frame(Sample_ID = colnames(asv),
                 RA_rep = colSums(asv))

design <- design[, colnames(design) %in% c("Sample_ID", "Host.Species",
                                           "Compartment", "Soil")]

RA <- merge(RA0, design)
ra_file <- "Figure_S7_2a.txt"
write.table(RA, ra_file, quote = F, sep = "\t", row.names = F)

RA_p1 <- RA[RA$Compartment %in%
            c("soil", "rhizosphere", "rhizoplane", "root"), ]

p1 <- ggplot(RA_p1, aes(x = RA_rep, y = Compartment, fill = Compartment),
             color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species),
                size = 1.4, color = "gray48", alpha = 0.6,
                position = position_jitter(width = 0.1)) +
    scale_fill_manual(values = c_Com) +
    scale_color_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = c("soil", "rhizosphere",
                                "rhizoplane", "root")) +
    xlim(c(0,1)) + main_theme +
    theme(legend.position = "none") +
    coord_flip()

host_lst <- c("zea_mays_l", "lotus_japonicus",
              "arabidopsis_thaliana", "chlamydomonas_reinhardtii")

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)
RA2 <- merge(RA0, design)

p3 <- ggplot(RA2, aes(x = RA_rep, y = Group, fill = Compartment,
                     shape = Host.Species), color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species, color = Compartment),
                size = 1.4, alpha = 0.6,
                position = position_jitter(width = 0.1)) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = plot_order, labels = bac_labels) +
    main_theme + xlim(c(0,1)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust=1, size = 8)) +
    coord_flip()

## for fungi
asv_file <- "../../00.data/Fun_ASV_raref.rds"
rep_file <- "../../00.data/Fun_ASV_rep.rds"
design_file <- "../../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")

asv <- asv[rownames(asv) %in% rownames(rep), ]

RA0 <- data.frame(Sample_ID = colnames(asv),
                 RA_rep = colSums(asv))

design <- design[, colnames(design) %in% c("Sample_ID", "Host.Species",
                                           "Compartment", "Soil")]

RA <- merge(RA0, design)
ra_file <- "Figure_S7_2b.txt"
write.table(RA, ra_file, quote = F, sep = "\t", row.names = F)

RA_p2 <- RA[RA$Compartment %in%
            c("soil", "rhizosphere", "rhizoplane", "root"), ]

p2 <- ggplot(RA_p2, aes(x = RA_rep, y = Compartment, fill = Compartment),
             color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species),
                size = 1.4, color = "gray48", alpha = 0.7,
                position = position_jitter(width = 0.1)) +
    scale_fill_manual(values = c_Com) +
    scale_color_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = c("soil", "rhizosphere",
                                "rhizoplane", "root")) +
    xlim(c(0,1)) + main_theme +
    theme(legend.position = "none") +
    coord_flip()

design <- design[design$Host.Species %in% host_lst, ]
design$Group <- paste0(design$Host.Species, "_", design$Soil, "_",
                       design$Compartment)
RA4 <- merge(RA0, design)

p4 <- ggplot(RA4, aes(x = RA_rep, y = Group, fill = Compartment,
                     shape = Host.Species), color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
    geom_jitter(aes(shape = Host.Species, color = Compartment),
                size = 1.4, alpha = 0.7,
                position = position_jitter(width = 0.1)) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = fun_plot_order, labels = fun_labels) +
    main_theme + xlim(c(0,1)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust=1, size = 8)) +
    coord_flip()

p12 <- plot_grid(p1, p2, nrow = 1, labels = 'auto')

p <- plot_grid(p12, p3, p4, ncol = 1, labels = c('', 'c', 'd'), rel_heights= c(1.4, 1,1))

ggsave("Figure_S7_2.pdf", height = 10, width = 10)

