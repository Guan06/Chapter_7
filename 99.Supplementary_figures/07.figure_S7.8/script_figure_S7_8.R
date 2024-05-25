#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
source("../../01.common_scripts/plot_setting.R")
source("../../01.common_scripts/pcoa_function.R")

design_file <- "../../00.data/design_854.txt"
net_dis_file <- "../../00.data/CAS_spearman_ap_dis.rds"

net_dis <- readRDS(net_dis_file)
design <- read.table(design_file, header = T, sep = "\t")

## r2
bio_factors <- c("Compartment", "Soil.Batch",
                  "Host.Species", "Host.Genotype", "Conditions")
tech_factors <- c("Study", "Run")

print(get_r2(net_dis, design, group = bio_factors))
print(get_r2(net_dis, design, group = c(bio_factors, tech_factors)))

dmr <- cmdscale(net_dis, k = 4, eig = T)
p <- pcoa(dmr, design, 12, "Compartment", "Host.Species", size = 2.4)
p <- p + theme(legend.position = "none")

ggsave("Figure_S7_8.pdf", p, width = 5, height = 5)
