#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

#library(mina)

library(devtools)
load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/mina/")

source("function_grp_sum_mean.R")
bs_pm_folder <- "../12.figure_7.14/CAS_vs_Root/"

### for network cluster features
grp_file <- "../00.data/CAS_spearman_ap_tab.txt"
out_dir <- "./network_cluster_cmp/"
out_dir_mean <- "./network_cluster_cmp_mean/"

grp <- read.table(grp_file, header = T, sep = "\t")

grp$Group <- grp$Color <- grp$Cluster
grp$Group <- as.factor(grp$Group)

if (!dir.exists(out_dir)) dir.create(out_dir)
net_grp_cmp(bs_pm_folder, "ja0", dir = out_dir, grp)

grp <- unique(grp[, c("Group", "Cluster_size", "Color")])

if (!dir.exists(out_dir_mean)) dir.create(out_dir_mean)
get_cmp_sum_mean_grp(out_dir, pattern = "ja0", out_dir_mean)

dat3 <- add_FC(paste0(out_dir_mean, "/ja0_arabidopsis_thaliana_root_vs_CAS.txt"),
               group = grp)
dat4 <- add_FC(paste0(out_dir_mean, "/ja0_CAS_vs_lotus_japonicus_root.txt"),
               group = grp)

write.table(dat3, "NetCls_CAS_vs_At_root.txt",
            quote = F, sep = "\t", row.names = F)
write.table(dat4, "NetCls_CAS_vs_Lj_root.txt",
            quote = F, sep = "\t", row.names = F)

### for family features
grp_file <- "CAS_ASV_1774_tax.txt"
out_dir <- "./family_cmp/"
out_dir_mean <- "./family_cmp_mean/"

grp <- read.table(grp_file, header = T, sep = "\t")
grp$ID <- rownames(grp)
grp$Group <- grp$Family
grp$Group[is.na(grp$Group)] <- "Unassigned"
grp$Group <- as.factor(grp$Group)
grp$Color <- grp$Phylum

if (!dir.exists(out_dir)) dir.create(out_dir)
net_grp_cmp(bs_pm_folder, "ja0", dir = out_dir, grp)

grp <- unique(grp[, c("Group", "Color")])

if (!dir.exists(out_dir_mean)) dir.create(out_dir_mean)
get_cmp_sum_mean_grp(out_dir, pattern = "ja0", out_dir_mean)

dat1 <- add_FC(paste0(out_dir_mean, "/ja0_arabidopsis_thaliana_root_vs_CAS.txt"),
               group = grp)
dat2 <- add_FC(paste0(out_dir_mean, "/ja0_CAS_vs_lotus_japonicus_root.txt"),
               group = grp)

write.table(dat1, "Family_CAS_vs_At_root.txt",
            quote = F, sep = "\t", row.names = F)
write.table(dat2, "Family_CAS_vs_Lj_root.txt",
            quote = F, sep = "\t", row.names = F)
