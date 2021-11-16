#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(ggplot2)
library(cowplot)
source("../01.common_scripts/plot_setting.R")

tax_file <- "../00.data/Bac_ASV_tax.txt"
rep_file <- "../00.data/Bac_ASV_rep.rds"

tax <- read.table(tax_file, header = T, sep = "\t")
rep <- readRDS(rep_file)

tax <- tax[rownames(tax) %in% rownames(rep), ]
tax$Node <- rownames(tax)
tax[is.na(tax)] <- "Unassigned"
tax <- tax[, c("Node", "Phylum")]

plot_cmp <- function(x) {
    this_cmp_file <- paste0("./node_cmp_stat/", x, ".rds")
    this_cmp <- readRDS(this_cmp_file)

    this_cmp$Diff <- as.numeric(this_cmp$Diff)
    this_cmp$Avg_BS <- as.numeric(this_cmp$Avg_BS)
    this_cmp$Avg_PM <- as.numeric(this_cmp$Avg_PM)

    this_cmp$Diff_abs <- abs(this_cmp$Diff)
    this_cmp <- this_cmp[order(this_cmp$Diff_abs, decreasing = T), ]

    this_cmp <- this_cmp[1:300, ]

    dat <- merge(this_cmp, tax)
    dat$FC <- ifelse(dat$Avg_BS > dat$Avg_PM,
                (dat$Avg_BS / dat$Avg_PM - 1),
                (1 - dat$Avg_PM/dat$Avg_BS))
    
    p <- ggplot(dat, aes(x = FC, Avg_BS)) +
     geom_point(aes(size = abs(Diff), color = Phylum, fill = Phylum),
                shape = 21, alpha = 0.8) +
        main_theme +
        scale_color_manual(values = c_phyla) +
        scale_fill_manual(values = c_phyla) +
        xlim(-5, 5) +
        scale_size(range = c(0, 5)) +
        theme(legend.position = "bottom", axis.line.y = element_blank()) +
        labs(y = x) +
        geom_vline(xintercept = 0, color = "gray")
}

p1 <- plot_cmp("arabidopsis_thaliana_root_vs_CAS_contrast")
#p2 <- plot_cmp("arabidopsis_thaliana_root_vs_CAS_ja")
p3 <- plot_cmp("arabidopsis_thaliana_root_vs_CAS_ja0")

p4 <- plot_cmp("CAS_vs_lotus_japonicus_root_contrast")
#p5 <- plot_cmp("CAS_vs_lotus_japonicus_root_ja")
p6 <- plot_cmp("CAS_vs_lotus_japonicus_root_ja0")

p1_legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

#p2_legend <- get_legend(p2)
#p2 <- p2 + theme(legend.position = "none")

p3_legend <- get_legend(p3)
p3 <- p3 + theme(legend.position = "none")

p4_legend <- get_legend(p4)
p4 <- p4 + theme(legend.position = "none")

#p5_legend <- get_legend(p5)
#p5 <- p5 + theme(legend.position = "none")

p6_legend <- get_legend(p6)
p6 <- p6 + theme(legend.position = "none")

#p <- plot_grid(p1, p4, p2, p5, p3, p6, nrow = 3, labels = "auto")
p <- plot_grid(p1, p4, p3, p6, nrow = 2, labels = "auto")
ggsave("Figure_7_14.pdf", height = 6, width = 10)

ggsave("Figure_7_14_legend_a.pdf", p1_legend, height = 1.5, width = 10)
ggsave("Figure_7_14_legend_b.pdf", p4_legend, height = 1.5, width = 10)
ggsave("Figure_7_14_legend_c.pdf", p3_legend, height = 1.5, width = 10)
ggsave("Figure_7_14_legend_d.pdf", p6_legend, height = 1.5, width = 10)


