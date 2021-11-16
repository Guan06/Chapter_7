#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
source("../../01.common_scripts/plot_setting.R")

bs_pm_folder <- "../../11.figure_7.13/CAS_bs_pm/"
tax_file <- "../../00.data/Bac_ASV_tax.txt"
design_file <- "../../00.data/design_854.txt"
fam_file <- "../../13.figure_7.15/aRA_occupancy_families.txt"

### get the list of x_order
fam <- read.table(fam_file, header = T, sep = "\t")
df_sort <- fam %>% group_by(Taxonomy) %>% summarise(avg_RA = mean(RA))
df_sort <- df_sort[order(df_sort$avg_RA, decreasing = T), ]
## average RA > 0.5 %
df_sort <- df_sort[df_sort$avg_RA > 0.005, ]
df_sort <- df_sort[df_sort$Taxonomy != "Unassigned", ]

x_order <- as.character(df_sort$Taxonomy)
###

tax <- read.table(tax_file, header = T, sep = "\t")

pm1_files <- sort(list.files(bs_pm_folder, pattern = "_pm1.rds",
                                 full.names = TRUE))
pm2_files <- sort(list.files(bs_pm_folder, pattern = "_pm2.rds",
                                 full.names = TRUE))

## function for getting node connectivity, i.e., sum of |edge strength|
get_node_cnt <- function(x) {
    pm1 <- readRDS(pm1_files[i])
    pm2 <- readRDS(pm2_files[i])
    group_mn <- strsplit(basename(pm1_files[i]), "_pm1.rds")[[1]][1]
    group_m <- strsplit(group_mn, "_vs_")[[1]][1]
    group_n <- strsplit(group_mn, "_vs_")[[1]][2]
    
    y_pm <- list()
    y_pm[1] <- pm1
    y_pm[2] <- pm2

    this_m <- y_pm[[1]]
    this_n <- y_pm[[2]]

    pm_len <- length(this_m)
    cnt <- c()
    for (j in 1 : pm_len) {
            adj_m <- unlist(this_m[[j]])
            adj_n <- unlist(this_n[[j]])
            adj_m[is.na(adj_m)] <- 0
            adj_n[is.na(adj_n)] <- 0

            this_m_cnt <- rowSums(abs(adj_m))
            this_n_cnt <- rowSums(abs(adj_n))

            this_cnt <- data.frame(Node = c(rownames(adj_m), rownames(adj_n)),
                                   Network = c(rep(group_m, nrow(adj_m)),
                                           rep(group_n, nrow(adj_n))),
                                   Connectivity = c(this_m_cnt, this_n_cnt))
            cnt <- rbind(cnt, this_cnt)
    }
    return(cnt)
}

len <- length(pm1_files)
CNT <- c()
for (i in 1:len) {
    CNT <- rbind(CNT, get_node_cnt(i))
}

## get the mean CNT
CNT$Group <- paste0(CNT$Network, "_", CNT$Node)
design <- unique(CNT[, c("Group", "Network", "Node")])
CNT2 <- as.data.frame(CNT %>%
                      group_by(Group) %>%
                      summarise(Mean_CNT = mean(Connectivity)))
CNT2 <- merge(CNT2, design)

## read in taxonomy file
tax <- read.table(tax_file, header = T, sep = "\t")
tax <- data.frame(Node = rownames(tax), Family = tax$Family)
tax <- tax[tax$Family %in% x_order, ]

## merge and reformat data frame for plotting
dat <- merge(CNT2, tax)
dat <- dat[, colnames(dat) != "Group"]
dat_wide <- pivot_wider(dat, names_from = Network, values_from = Mean_CNT)
cnt <- as.data.frame(dat_wide)

design <- read.table(design_file, header = T, sep = "\t")
des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "CAS"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)

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

p <- ggplot(dat, aes(x = Family, y = Connectivity, color = Compartment)) +
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

ggsave("Figure_S7_14.pdf", p, width = 8, height = 3.5)
