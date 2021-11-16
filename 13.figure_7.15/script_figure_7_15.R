#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
source("../01.common_scripts/plot_setting.R")

asv_file <- "../00.data/Bac_ASV_raref.rds"
rep_file <- "../00.data/Bac_ASV_rep.rds"
tax_file <- "../00.data/Bac_ASV_tax.txt"
design_file <- "../00.data/design_854.txt"
cnt_file <- "./CAS_node_connectivity_mean_tax.txt"

asv <- readRDS(asv_file)
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")
tax <- read.table(tax_file, header = T, sep = "\t")
cnt <- read.table(cnt_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)

g1 <- c("soil", "arabidopsis_thaliana_rhizosphere",
        "arabidopsis_thaliana_rhizoplane", "arabidopsis_thaliana_root")
g2 <- c("soil", "chlamydomonas_reinhardtii_phycosphere",
        "lotus_japonicus_rhizosphere", "lotus_japonicus_root")

TAX_DF <- c()
g_lst <- unique(des$Group)

for (g in g_lst) {
    this_des <- des[des$Group == g, ]
    this_asv <- asv[rownames(asv) %in% rownames(rep), colnames(asv) %in% this_des$Sample_ID]
    this_asv <- this_asv[rowSums(this_asv) > 0, colSums(this_asv) > 0]

    this_fill <- unique(this_des$Compartment)
    this_tax <- tax[rownames(tax) %in% rownames(this_asv), ]
    this_tax$Family[is.na(this_tax$Family)] <- "Unassigned"

    TAX <- as.factor(this_tax$Family)
    tab_TAX <- apply(this_asv, 2, function(x) rowsum(x, TAX))
    rownames(tab_TAX) <- levels(TAX)

	tax1_df <- data.frame(Taxonomy = rownames(tab_TAX),
	                    Occupancy = rowSums(tab_TAX > 0) / ncol(tab_TAX),
	                    Sample_number = ncol(tab_TAX),
	                    RA = rowSums(tab_TAX) / ncol(tab_TAX),
	                    Group = g,
                        Compartment = this_fill,
	                    Rank = "Family",
	                    Tax_number = nrow(tab_TAX))
    TAX_DF <- rbind(TAX_DF, tax1_df)
}

write.table(TAX_DF, "aRA_occupancy_families.txt", quote = F, sep = "\t", row.names = F)

## sort by the average RA of each family in each condition
df_sort <- TAX_DF %>% group_by(Taxonomy) %>% summarise(avg_RA = mean(RA))
df_sort <- df_sort[order(df_sort$avg_RA, decreasing = T), ]
## average RA > 0.5 %
df_sort <- df_sort[df_sort$avg_RA > 0.005, ]
df_sort <- df_sort[df_sort$Taxonomy != "Unassigned", ]

x_order <- as.character(df_sort$Taxonomy)

dat <- TAX_DF[TAX_DF$Taxonomy %in% x_order, ]
dat2 <- dat[dat$Group %in% g1, ]
dat1 <- dat[dat$Group %in% g2, ]

dat1$Compartment <- factor(dat1$Compartment,
                              levels = c("soil", "phycosphere",
                                         "rhizosphere", "root"),
                              ordered = "TRUE")

dat2$Compartment <- factor(dat2$Compartment,
                           levels = c("soil", "rhizosphere",
                                      "rhizoplane", "root"),
                           ordered = "TRUE")

p1_a <- ggplot(dat1, aes(fill = Compartment)) +
    geom_bar(aes(x = Taxonomy, y = Occupancy),
             stat = 'identity', position = "dodge", color = NA) +
    main_theme +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    labs(x = "Family", y = "Occupancy") +
    theme(legend.position = "none",
            axis.ticks = element_blank(),
            axis.line.x = element_blank(),
            axis.text = element_text(size = 6),
            axis.text.x = element_text(angle = 90,
                                       vjust = 0.5,
                                       hjust=1),
            axis.title = element_text(size = 6),
            plot.title = element_text(size = 6))

p1_b <- ggplot(dat1, aes(x = Taxonomy, y = RA,
                         color = Compartment)) +
    geom_line(aes(group = Compartment), size = 0.8) +
    geom_point() +
    scale_y_continuous(limits = c(0, 0.6), position = "right") +
    main_theme +
     scale_color_manual(values = c_Com) +
    scale_x_discrete(limits = x_order) +
    theme_minimal_hgrid(color = "gray85", line_size = 0.5) +
    labs(y = "RA") +
    theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.text.y = element_text(size = 6),
              axis.title.y = element_text(size = 6),
              axis.title.x = element_blank())

aligned_plots <- align_plots(p1_a, p1_b, align="hv", axis="tblr")
p1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

p2_a <- ggplot(dat2, aes(fill = Compartment)) +
    geom_bar(aes(x = Taxonomy, y = Occupancy),
             stat = 'identity', position = "dodge", color = NA) +
    main_theme +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_Com) +
    scale_fill_manual(values = c_Com) +
    labs(x = "Family", y = "Occupancy") +
    theme(legend.position = "none",
            axis.ticks = element_blank(),
            axis.line.x = element_blank(),
            axis.text = element_text(size = 6),
            axis.text.x = element_text(angle = 90,
                                       vjust = 0.5,
                                       hjust=1),
            axis.title = element_text(size = 6),
            plot.title = element_text(size = 6))

p2_b <- ggplot(dat2, aes(x = Taxonomy, y = RA,
                         color = Compartment)) +
    geom_line(aes(group = Compartment), size = 0.8) +
    geom_point() +
    scale_y_continuous(limits = c(0, 0.6), position = "right") +
    main_theme +
     scale_color_manual(values = c_Com) +
    scale_x_discrete(limits = x_order) +
    theme_minimal_hgrid(color = "gray85", line_size = 0.5) +
    labs(y = "RA") +
    theme(legend.position = "none",
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.text.y = element_text(size = 6),
              axis.title.y = element_text(size = 6),
              axis.title.x = element_blank())

aligned_plots <- align_plots(p2_a, p2_b, align="hv", axis="tblr")
p2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

### plot the panel c
tax <- tax[tax$Family %in% x_order, ]
cnt <- cnt[cnt$Node %in% rownames(tax), ]

tax2 <- data.frame(Node = rownames(tax), Family = tax$Family)
cnt <- merge(cnt, tax2)

dat_tax <- cnt %>% group_by(Family) %>%
        summarise(At_rhizos = mean(arabidopsis_thaliana_rhizosphere),
                  CAS = mean(CAS),
                  Lj_rhizos = mean(lotus_japonicus_rhizosphere),
                  Lj_root = mean(lotus_japonicus_root),
                  At_rhizop = mean(arabidopsis_thaliana_rhizoplane),
                  At_root = mean(arabidopsis_thaliana_root),
                  Cr_phyco = mean(chlamydomonas_reinhardtii_phycosphere))

dat_tax <- as.data.frame(dat_tax)

dat <- dat_tax %>%
    pivot_longer(At_rhizos:Cr_phyco,
                 names_to = "Group",
                 values_to = "Connectivity")

dat <- as.data.frame(dat)

dat_des <- data.frame(Group = c("At_rhizos", "CAS", "Lj_rhizos", "Lj_root",
                                "At_rhizop", "At_root", "Cr_phyco"),
                      Compartment = c("rhizosphere", "soil", "rhizosphere",
                                      "root", "rhizoplane", "root",
                                      "phycosphere"),
                      Host.Species = c("arabidopsis_thaliana", "host_soil",
                                       rep("lotus_japonicus", 2),
                                       rep("arabidopsis_thaliana", 2),
                                       "chlamydomonas_reinhardtii"))

dat <- merge(dat, dat_des)

p3 <- ggplot(dat, aes(x = Family, y = Connectivity, color = Compartment)) +
    geom_point(aes(shape = Host.Species)) +
    geom_line(aes(group = Group), size = 0.8) +
    scale_x_discrete(limits = x_order) +
    scale_color_manual(values = c_Com) +
    scale_shape_manual(values = s_Host) +
    scale_fill_manual(values = c_Com) +
    labs(x = "Family", y = "Connectivity") +
    main_theme +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
            axis.text.x = element_text(angle = -60,
                                       vjust = 0.5,
                                       hjust=0),
            axis.title = element_text(size = 6),
            plot.title = element_text(size = 6))

p <- plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "l", labels = 'auto')
ggsave("Figure_7_15.pdf", p, width = 8, height = 10)
