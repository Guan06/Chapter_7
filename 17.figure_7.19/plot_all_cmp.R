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
source("cmp.lst")
######### for family ##########################################################
dis_mean_file <- "../15.figure_7.17/Family_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Family)

x_order <- c("All_BS", family, "All_PM")

dis_folder <- "./Family/spectra_dis/"
out_dir <- "./Family/plots/"
if (!dir.exists(out_dir)) dir.create(out_dir)

lapply(cmp_lst1, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst1"))
lapply(cmp_lst2, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst2"))
lapply(cmp_lst3, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst3"))
lapply(cmp_lst4, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst4"))

######## for network cluster ##################################################
dis_mean_file <- "../15.figure_7.17/Cluster_single_pm_dis_mean.txt"

dis_mean <- read.table(dis_mean_file, header = T, sep = "\t")
dis_mean <- dis_mean[order(dis_mean$Decrease_Mean, decreasing = T), ]
family <- as.character(dis_mean$Cluster)

x_order <- c("All_BS", family, "All_PM")

dis_folder <- "./Network_cluster/spectra_dis/"
out_dir <- "./Network_cluster/plots/"
if (!dir.exists(out_dir)) dir.create(out_dir)

lapply(cmp_lst1, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst1"))
lapply(cmp_lst2, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst2"))
lapply(cmp_lst3, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst3"))
lapply(cmp_lst4, FUN=function(x) plot_accu(x, dis_folder, out_dir, "lst4"))
