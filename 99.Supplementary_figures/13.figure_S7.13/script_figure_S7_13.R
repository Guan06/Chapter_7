#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(ggplot2)
library(dplyr)
library(ggridges)
source("../../01.common_scripts/plot_setting.R")

tax_file <- "../../00.data/Bac_ASV_tax.txt"
design_file <- "../../00.data/design_854.txt"
fam_file <- "../../13.figure_7.15/aRA_occupancy_families.txt"
rep_file <- "../../00.data/Bac_ASV_rep.rds"
asv_file <- "../../00.data/Bac_ASV_raref.rds"

tax <- read.table(tax_file, header = T, sep = "\t")
design <- read.table(design_file, header = T, sep = "\t")
rep <- readRDS(rep_file)
asv <- readRDS(asv_file)

### get the list of families
fam <- read.table(fam_file, header = T, sep = "\t")
df_sort <- fam %>% group_by(Taxonomy) %>% summarise(avg_RA = mean(RA))
df_sort <- df_sort[order(df_sort$avg_RA, decreasing = T), ]
## average RA > 0.5 %
df_sort <- df_sort[df_sort$avg_RA > 0.005, ]
df_sort <- df_sort[df_sort$Taxonomy != "Unassigned", ]

lst <- as.character(df_sort$Taxonomy)

tax <- tax[tax$Family %in% lst, ]
asv_lst <- intersect(rownames(tax), rownames(rep))
asv <- asv[rownames(asv) %in% asv_lst, ]
asv <- asv[rowSums(asv) > 0, colSums(asv) > 0]

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)

ra <- data.frame(Sample_ID = colnames(asv),
                 RA = colSums(asv))
ra <- merge(ra, des)

this_order <- c("soil",
                "lotus_japonicus_rhizosphere",
                "lotus_japonicus_root",
                "chlamydomonas_reinhardtii_phycosphere",
                "arabidopsis_thaliana_rhizosphere",
                "arabidopsis_thaliana_rhizoplane",
                "arabidopsis_thaliana_root")

this_label <- c("CAS", "Lj_rhizos", "Lj_root",
                "Chlamy_phyco",
                "At_rhizos",
                "At_rhizop", "At_root")

p <- ggplot(ra, aes(x = RA, y = Group, fill = Compartment,
                     shape = Host.Species), color = "gray48") +
    geom_density_ridges(scale = 0.5, quantile_lines = TRUE, quantiles = 2,
                        color = "gray48") +
#    geom_point(aes(shape = Host.Species, color = Compartment),
#                size = 1, alpha = 0.7) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_y_discrete(limits = this_order, labels = this_label) +
    main_theme +
    labs(y = "") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust=1, size = 8)) +
    coord_flip()

ggsave("Figure_S7_13.pdf", height = 3)
