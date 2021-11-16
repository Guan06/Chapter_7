#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
library(cowplot)

source("../../01.common_scripts/plot_setting.R")

bacteria_stat <- read.table("./dis_stat_spectra_bacteria_33.txt", header = T, sep = "\t")
fungi_stat <- read.table("./dis_stat_spectra_fungi_33.txt", header = T, sep = "\t")

bac <- new("mina")
fun <- new("mina")

dis_stat(bac) <- bacteria_stat
dis_stat(fun) <- fungi_stat

p1 <- net_dis_plot(bac)
p2 <- net_dis_plot(fun)

p <- plot_grid(p1, p2, nrow = 1)
ggsave("Figure_S7_8.pdf", p)
