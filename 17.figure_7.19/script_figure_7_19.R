#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

###############################################################################
## run script get_accu_pm.R and get_accu_dis.R to get the accumulative
## permutation of groups ordered by the single group permutation distance
## change before plotting.
###############################################################################

library(cowplot)
library(mina)
source("function_plot_accu.R")
source("../01.common_scripts/plot_setting.R")

######### for family ##########################################################
dis_mean_file <- "../15.figure_7.17/Family_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)

x_order <- c("All_BS", family, "All_PM")
dis_folder <- "./Family/spectra_dis/"

p1 <- plot_accu_indi("arabidopsis_thaliana_root_vs_CAS", dis_folder)
p3 <- plot_accu_indi("CAS_vs_lotus_japonicus_root", dis_folder)

######## for network cluster ##################################################
dis_mean_file <- "../15.figure_7.17/Cluster_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Cluster)

x_order <- c("All_BS", family, "All_PM")

dis_folder <- "./Network_cluster/spectra_dis/"

p2 <- plot_accu_indi("arabidopsis_thaliana_root_vs_CAS", dis_folder)
p4 <- plot_accu_indi("CAS_vs_lotus_japonicus_root", dis_folder)

p <- plot_grid(p1, p2, p3, p4, ncol = 1, labels = "auto",
               rel_heights = c(1.2, 1, 1.2, 1))
ggsave("Figure_7_19.pdf", p, width = 10, height = 12)
