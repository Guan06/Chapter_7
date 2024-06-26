#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
library(devtools)
#load_all("/biodata/dep_psl/grp_rgo/guan/Rpackage_Mina/mina/")

library(Hmisc)
source("../../01.common_scripts/plot_setting.R")
source("../../01.common_scripts/pcoa_function.R")

asv_file <- "../../00.data/Bac_ASV_rep.rds"
design_file <- "../../00.data/Bac_design_3809.txt"

## factors for R2 calculation
bio_factors <- c("Compartment", "Soil", "Soil.Batch",
                  "Host.Species", "Host.Genotype", "Conditions")
tech_factors <- c("Study", "Run")
###

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

net <- rcorr(t(asv), type = "spearman")
mat_cor <- net$r
mat_p <- net$P
mat_cor[mat_p > 0.05] <- 0

#cls_ap <- net_cls(mat_cor, method = "ap", cutoff = 0)

ap <- apcluster(mat_cor, p = 0)
cls_ap <- re_format_AP(ap)

cls_ap_tab <- get_net_cls_tab(asv, cls_ap)
dis_ap <- com_dis(cls_ap_tab, method = "bray")

print("R2 of AP cluster(bio_factors) for bacterial community: ")
print(get_r2(dis_ap, design, group = bio_factors))
print("R2 of AP cluster(bio and tech facotrs) for bacterial community: ")
print(get_r2(dis_ap, design, group = c(bio_factors, tech_factors)))

ap_file <- "Figure_S7_4a_tab_spearman_ap.txt"
write.table(cls_ap, ap_file, quote = F, row.names = F, sep = "\t")

## diversity analysis based on network clusters
dmr_ap <- dmr(dis_ap)

p1 <- pcoa(dmr_ap, design, 12, "Compartment", "Host.Species")

## for Fungi
asv_file <- "../../00.data/Fun_ASV_rep.rds"
design_file <- "../../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

net <- rcorr(t(asv), type = "spearman")
mat_cor <- net$r
mat_p <- net$P
mat_cor[mat_p > 0.05] <- 0

#cls_ap <- net_cls(mat_cor, method = "ap", cutoff = 0)
ap <- apcluster(mat_cor, p = 0)
cls_ap <- re_format_AP(ap)

cls_ap_tab <- get_net_cls_tab(asv, cls_ap)
dis_ap <- com_dis(cls_ap_tab, method = "bray")

print("R2 of AP cluster(bio_factors) for fungal community: ")
print(get_r2(dis_ap, design, group = bio_factors))
print("R2 of AP cluster(bio and tech facotrs) for fungal community: ")
print(get_r2(dis_ap, design, group = c(bio_factors, tech_factors)))

ap_file <- "Figure_S7_4_tab_spearman_ap.txt"
write.table(cls_ap, ap_file, quote = F, row.names = F, sep = "\t")

## diversity analysis based on network clusters
dmr_ap <- dmr(dis_ap)

p2 <- pcoa(dmr_ap, design, 12, "Compartment", "Host.Species")


## put panels together
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

p <- plot_grid(p1, p2, nrow = 1)
ggsave("Figure_S7_4.pdf", p, width = 9)
