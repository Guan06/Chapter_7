#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(dplyr)
library(mina)

source("../01.common_scripts/plot_setting.R")
source("./function_plot_CAS.R")

### for ASV
asv_file <- "../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../00.data/design_854.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

## reformat design
des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)
###

## for taxonomy at the Family level
tax_file <- "../00.data/Bac_ASV_tax.txt"
tax <- read.table(tax_file, header = T, sep = "\t")
tax <- tax[rownames(tax) %in% rownames(asv), ]

tax[is.na(tax$Family), ]$Family <- "Unassigned"
tax$Family <- as.factor(tax$Family)

tab_tax <- apply(asv, 2, function(x) rowsum(x, as.factor(tax$Family)))
rownames(tab_tax) <- levels(tax$Family)

## for network cluster table
cls_ap_file <- "../00.data/CAS_spearman_ap_tab.txt"
cls_ap <- read.table(cls_ap_file, header = T, sep = "\t")
cls_ap_tab <- get_net_cls_tab(asv, cls_ap)

### plot for asv, family and network cluster
p1 <- plot_CAS(asv)
p2 <- plot_CAS(tab_tax)
p3 <- plot_CAS(cls_ap_tab)

p <- plot_grid(p1, p2, p3, nrow = 3, ncol = 1,
               align = "v", axis = "l", labels = 'auto')

ggsave("Figure_7_12.pdf", p, height = 3, width = 4)
