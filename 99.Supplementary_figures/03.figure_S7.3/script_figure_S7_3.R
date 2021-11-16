#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(parallelDist)
source("../../01.common_scripts/plot_setting.R")
source("../../01.common_scripts/pcoa_function.R")

## for panel A
asv_file <- "../../00.data/Bac_ASV_rep.rds"
design_file <- "../../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- as.matrix(bc)

design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_a <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
bac_rep_bc <- bc

## for panel B
asv_file <- "../../00.data/Fun_ASV_rep.rds"
design_file <- "../../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
bc <- parDist(t(asv), method = "bray", threads = 40)
bc <- as.matrix(bc)

design <- read.table(design_file, header = T, sep = "\t")
dmr <- cmdscale(bc, k = 4, eig = T)
p_b <- pcoa(dmr, design, 12, "Compartment", "Host.Species")
p_b <- p_b + theme(legend.position = "none")
fun_rep_bc <- bc

## for panel C
bac_bc_file <- "../../00.data/Bac_ASV_bc.rds"
fun_bc_file <- "../../00.data/Fun_ASV_bc.rds"

bac_bc <- readRDS(bac_bc_file)
fun_bc <- readRDS(fun_bc_file)
c_bac <- bac_rep_bc / bac_bc
c_fun <- fun_rep_bc / fun_bc

flatten_mat <- function(x) {
    ut <- upper.tri(x)
    data.frame(
        sample_1 = rownames(x)[row(x)[ut]],
        sample_2 = rownames(x)[col(x)[ut]],
        C = x[ut]
    )
}

c_bac_df <- flatten_mat(c_bac)
c_fun_df <- flatten_mat(c_fun)

c_bac_df$Kingdom <- "Bacteria"
c_fun_df$Kingdom <- "Fungi"
c_df <- rbind(c_bac_df, c_fun_df)

p_c <- ggplot(c_df, aes(Kingdom, C)) +
    geom_violin(color = "gray64", fill = "gray48") +
    ylim(0, 1.2) +
    main_theme

## format panels
library(cowplot)
p_legend <- get_legend(p_a)
p_a <- p_a + theme(legend.position = "none")

p_ab <- plot_grid(p_a, p_b, p_legend, nrow = 1,
                   labels = c('a', 'b', ''),
                   rel_widths = c(1, 1, 0.3),
                   rel_heights = c(1, 1, 1))

p <- plot_grid(p_ab, p_c, nrow = 1, labels = c('', 'c'),
               rel_widths = c(1, 0.2))

ggsave("Figure_S7_3.pdf", p, width = 14, height = 5)
