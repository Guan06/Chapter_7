#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(cowplot)
library(dplyr)

source("../01.common_scripts/plot_setting.R")

###############################################################################
## run script get_grp_size.R to get the size of families and network clusters
## under each condition before plotting.
###############################################################################

#### function to plot the scatters from data frame
plot_feature <- function(x, cmp, color){
    dat <- merge(x, decrease)
    
    t <- cor.test(dat[[cmp]], dat$Decrease, method = "spearman", exact = FALSE)
    edge <- t$estimate
    sig <- t$p.value

    edge <- round(edge, 2)
    sig <- round(sig, 4)

    p <- ggplot(dat, aes_string(cmp, "Decrease")) +
          scale_x_log10() +
        geom_smooth(method = "lm", se = FALSE,
                    color = "#017F97") +
        geom_point(color = color, alpha = 0.8) + main_theme +
        theme(text = element_text(colour = "black", size = 8),
              axis.text = element_text(size = 8)) +
        labs(title = paste0("r = ", edge, " P = ", sig))
}
#######################

######### for family ##########################################################
dis_file <- "../15.figure_7.17/Figure_7_17a.txt"
size_file <- "size_families.txt"
ra_occu_file <- "../13.figure_7.15/aRA_occupancy_families.txt"
ja0_at_file <- "../14.figure_7.16/Family_CAS_vs_At_root.txt"
ja0_lj_file <- "../14.figure_7.16/Family_CAS_vs_Lj_root.txt"

dis <- read.table(dis_file, header = T, sep = "\t")
colnames(dis)[1] <- "Group"
dis <- dis[dis$Group != "Unassigned", ]

size <- read.table(size_file, header = T, sep = "\t")
ra_occu <- read.table(ra_occu_file, header = T, sep = "\t")
ra_occu <- ra_occu[, c("Taxonomy", "Occupancy", "RA", "Group")]
colnames(ra_occu) <- c("Group", "Occupancy", "RA", "Compartment")

at <- read.table(ja0_at_file, header = T, sep = "\t")
lj <- read.table(ja0_lj_file, header = T, sep = "\t")

### Plot for panel a: CAS_vs_At_root, Family
at_color <- "#1B9E77"

dis_at <- dis[dis$Comparison == "arabidopsis_thaliana_root_vs_CAS", ]
decrease <- dis_at[, c("Group", "Decrease")]

size_CAS <- size[size$Compartment == "soil", ]
size_at <- size[size$Compartment == "arabidopsis_thaliana_root", ]

p0 <- plot_feature(size_CAS, "Number_of_ASVs", color = at_color) +
    labs(x = "Number of ASVs in CAS")
p1 <- plot_feature(size_at, "Number_of_ASVs", color = at_color) +
    labs(x = "Number of ASVs in At root")

ra <- ra_occu[, c("Group", "RA", "Compartment")]
ra_CAS <- ra[ra$Compartment == "soil", ]
ra_At <- ra[ra$Compartment == "arabidopsis_thaliana_root", ]

ra_CAS_norm <- ra_CAS %>% inner_join(size_CAS, by = c("Group", "Compartment"))
ra_CAS_norm$aRA_norm <- ra_CAS_norm$RA / ra_CAS_norm$Number_of_ASVs

ra_At_norm <- ra_At %>% inner_join(size_at, by = c("Group", "Compartment"))
ra_At_norm$aRA_norm <- ra_At_norm$RA / ra_At_norm$Number_of_ASVs

p2 <- plot_feature(ra_CAS, "RA", color = at_color) +
    labs(x = "aRA in CAS")

p3 <- plot_feature(ra_CAS_norm, "aRA_norm", color = at_color) +
    labs(x = "Normalized aRA in CAS")

p4 <- plot_feature(ra_At, "RA", color = at_color) +
    labs(x = "aRA in At root")

p5 <- plot_feature(ra_At_norm, "aRA_norm", color = at_color) +
    labs(x = "Normalized aRA in At root")

occu <- ra_occu[, c("Group", "Occupancy", "Compartment")]
occu_CAS <- occu[occu$Compartment == "soil", ]
occu_At <- occu[occu$Compartment == "arabidopsis_thaliana_root", ]

p6 <- plot_feature(occu_CAS, "Occupancy", color = at_color) +
    labs(x = "Occupancy in CAS")
p7 <- plot_feature(occu_At, "Occupancy", color = at_color) +
    labs(x = "Occupancy in At root")

at <- at[, c("Group", "Avg_BS")]
p8 <- plot_feature(at, "Avg_BS", color = at_color) +
    labs(x = "Jaccard distance between CAS and At root")

p_a <- plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8,
                 nrow = 1, align = "h", axis = "b",
                 labels = c('a', '', '', '', '', '', '', ''))

#ggsave("Figure_S7_16.pdf", p_a, height = 2, width = 16)

### Plot for panel c: CAS_vs_Lj_root,  Family
lj_color <- "#D95F02"

dis_lj <- dis[dis$Comparison == "CAS_vs_lotus_japonicus_root", ]
decrease <- dis_lj[, c("Group", "Decrease")]

size_CAS <- size[size$Compartment == "soil", ]
size_lj <- size[size$Compartment == "lotus_japonicus_root", ]

p0 <- plot_feature(size_CAS, "Number_of_ASVs", color = lj_color) +
    labs(x = "Number of ASVs in CAS")
p1 <- plot_feature(size_lj, "Number_of_ASVs", color = lj_color) +
    labs(x = "Number of ASVs in Lj root")

ra <- ra_occu[, c("Group", "RA", "Compartment")]
ra_CAS <- ra[ra$Compartment == "soil", ]
ra_Lj <- ra[ra$Compartment == "lotus_japonicus_root", ]

ra_CAS_norm <- ra_CAS %>% inner_join(size_CAS, by = c("Group", "Compartment"))
ra_CAS_norm$aRA_norm <- ra_CAS_norm$RA / ra_CAS_norm$Number_of_ASVs

ra_Lj_norm <- ra_Lj %>% inner_join(size_lj, by = c("Group", "Compartment"))
ra_Lj_norm$aRA_norm <- ra_Lj_norm$RA / ra_Lj_norm$Number_of_ASVs

p2 <- plot_feature(ra_CAS, "RA", color = lj_color) +
    labs(x = "aRA in CAS")

p3 <- plot_feature(ra_CAS_norm, "aRA_norm", color = lj_color) +
    labs(x = "Normalized aRA in CAS")

p4 <- plot_feature(ra_Lj, "RA", color = lj_color) +
    labs(x = "aRA in Lj root")

p5 <- plot_feature(ra_Lj_norm, "aRA_norm", color = lj_color) +
    labs(x = "Normalized aRA in Lj root")

occu <- ra_occu[, c("Group", "Occupancy", "Compartment")]
occu_CAS <- occu[occu$Compartment == "soil", ]
occu_Lj <- occu[occu$Compartment == "lotus_japonicus_root", ]

p6 <- plot_feature(occu_CAS, "Occupancy", color = lj_color) +
    labs(x = "Occupancy in CAS")
p7 <- plot_feature(occu_Lj, "Occupancy", color = lj_color) +
    labs(x = "Occupancy in Lj root")

lj <- lj[, c("Group", "Avg_BS")]
p8 <- plot_feature(lj, "Avg_BS", color = lj_color) +
    labs(x = "Jaccard distance between CAS and Lj root")

p_c <- plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8,
                 nrow = 1, align = "h", axis = "b",
                 labels = c('c', '', '', '', '', '', '', ''))


######### for network clusters ################################################
dis_file <- "../15.figure_7.17/Figure_7_17b.txt"
size_file <- "size_network_clusters.txt"
ra_occu_file <- "../14.figure_7.16/aRA_occupancy_network_clusters.txt"
ja0_at_file <- "../14.figure_7.16/NetCls_CAS_vs_At_root.txt"
ja0_lj_file <- "../14.figure_7.16/NetCls_CAS_vs_Lj_root.txt"

dis <- read.table(dis_file, header = T, sep = "\t")
colnames(dis)[1] <- "Group"
size <- read.table(size_file, header = T, sep = "\t")
ra_occu <- read.table(ra_occu_file, header = T, sep = "\t")
ra_occu <- ra_occu[, c("Group", "Occupancy", "RA", "Group.1")]
colnames(ra_occu) <- c("Group", "Occupancy", "RA", "Compartment")
ra_occu$Group <- paste0("Cluster_", ra_occu$Group)

at <- read.table(ja0_at_file, header = T, sep = "\t")
lj <- read.table(ja0_lj_file, header = T, sep = "\t")

### plot for panel b: CAS_vs_At_root, Network clusters
at_color <- "#1B9E77"

dis_at <- dis[dis$Comparison == "arabidopsis_thaliana_root_vs_CAS", ]
decrease <- dis_at[, c("Group", "Decrease")]

size_CAS <- size[size$Compartment == "soil", ]
size_CAS$Group <- paste0("Cluster_", size_CAS$Group)

size_at <- size[size$Compartment == "arabidopsis_thaliana_root", ]
size_at$Group <- paste0("Cluster_", size_at$Group)

p0 <- plot_feature(size_CAS, "Number_of_ASVs", color = at_color) +
    labs(x = "Number of ASVs in CAS")
p1 <- plot_feature(size_at, "Number_of_ASVs", color = at_color) +
    labs(x = "Number of ASVs in At root")

ra <- ra_occu[, c("Group", "RA", "Compartment")]
ra_CAS <- ra[ra$Compartment == "soil", ]

ra_At <- ra[ra$Compartment == "arabidopsis_thaliana_root", ]

ra_CAS_norm <- ra_CAS %>% inner_join(size_CAS, by = c("Group", "Compartment"))
ra_CAS_norm$aRA_norm <- ra_CAS_norm$RA / ra_CAS_norm$Number_of_ASVs

ra_At_norm <- ra_At %>% inner_join(size_at, by = c("Group", "Compartment"))
ra_At_norm$aRA_norm <- ra_At_norm$RA / ra_At_norm$Number_of_ASVs

p2 <- plot_feature(ra_CAS, "RA", color = at_color) +
    labs(x = "aRA in CAS")

p3 <- plot_feature(ra_CAS_norm, "aRA_norm", color = at_color) +
    labs(x = "Normalized aRA in CAS")

p4 <- plot_feature(ra_At, "RA", color = at_color) +
    labs(x = "aRA in At root")

p5 <- plot_feature(ra_At_norm, "aRA_norm", color = at_color) +
    labs(x = "Normalized aRA in At root")

occu <- ra_occu[, c("Group", "Occupancy", "Compartment")]
occu_CAS <- occu[occu$Compartment == "soil", ]
occu_At <- occu[occu$Compartment == "arabidopsis_thaliana_root", ]

p6 <- plot_feature(occu_CAS, "Occupancy", color = at_color) +
    labs(x = "Occupancy in CAS")
p7 <- plot_feature(occu_At, "Occupancy", color = at_color) +
    labs(x = "Occupancy in At root")

at <- at[, c("Group", "Avg_BS")]
at$Group <- paste0("Cluster_", at$Group)

p8 <- plot_feature(at, "Avg_BS", color = at_color) +
    labs(x = "Jaccard distance between CAS and At root")

p_b <- plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8,
                 nrow = 1, align = "h", axis = "b",
                 labels = c('b', '', '', '', '', '', '', ''))

### Plot for panel d: CAS_vs_Lj_root,  Network cluster
lj_color <- "#D95F02"

dis_lj <- dis[dis$Comparison == "CAS_vs_lotus_japonicus_root", ]
decrease <- dis_lj[, c("Group", "Decrease")]

size_CAS <- size[size$Compartment == "soil", ]
size_CAS$Group <- paste0("Cluster_", size_CAS$Group)
size_lj <- size[size$Compartment == "lotus_japonicus_root", ]
size_lj$Group <- paste0("Cluster_", size_lj$Group)

p0 <- plot_feature(size_CAS, "Number_of_ASVs", color = lj_color) +
    labs(x = "Number of ASVs in CAS")
p1 <- plot_feature(size_lj, "Number_of_ASVs", color = lj_color) +
    labs(x = "Number of ASVs in Lj root")

ra <- ra_occu[, c("Group", "RA", "Compartment")]
ra_CAS <- ra[ra$Compartment == "soil", ]
ra_Lj <- ra[ra$Compartment == "lotus_japonicus_root", ]

ra_CAS_norm <- ra_CAS %>% inner_join(size_CAS, by = c("Group", "Compartment"))
ra_CAS_norm$aRA_norm <- ra_CAS_norm$RA / ra_CAS_norm$Number_of_ASVs

ra_Lj_norm <- ra_Lj %>% inner_join(size_lj, by = c("Group", "Compartment"))
ra_Lj_norm$aRA_norm <- ra_Lj_norm$RA / ra_Lj_norm$Number_of_ASVs

p2 <- plot_feature(ra_CAS, "RA", color = lj_color) +
    labs(x = "aRA in CAS")

p3 <- plot_feature(ra_CAS_norm, "aRA_norm", color = lj_color) +
    labs(x = "Normalized aRA in CAS")

p4 <- plot_feature(ra_Lj, "RA", color = lj_color) +
    labs(x = "aRA in Lj root")

p5 <- plot_feature(ra_Lj_norm, "aRA_norm", color = lj_color) +
    labs(x = "Normalized aRA in Lj root")

occu <- ra_occu[, c("Group", "Occupancy", "Compartment")]
occu_CAS <- occu[occu$Compartment == "soil", ]
occu_Lj <- occu[occu$Compartment == "lotus_japonicus_root", ]

p6 <- plot_feature(occu_CAS, "Occupancy", color = lj_color) +
    labs(x = "Occupancy in CAS")
p7 <- plot_feature(occu_Lj, "Occupancy", color = lj_color) +
    labs(x = "Occupancy in Lj root")

lj <- lj[, c("Group", "Avg_BS")]
lj$Group <- paste0("Cluster_", lj$Group)

p8 <- plot_feature(lj, "Avg_BS", color = lj_color) +
    labs(x = "Jaccard distance between CAS and Lj root")

p_d <- plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8,
                 nrow = 1, align = "h", axis = "b",
                 labels = c('d', '', '', '', '', '', '', ''))

p <- plot_grid(p_a, p_b, p_c, p_d, nrow = 4)
ggsave("Figure_S7_16.pdf", p, width = 18, height = 8)
