#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

bs_pm_folder <- "../11.figure_7.13/CAS_bs_pm/"
tax_file <- "../00.data/Bac_ASV_tax.txt"

bs1_files <- sort(list.files(bs_pm_folder, pattern = "_bs1.rds",
                                 full.names = TRUE))
bs2_files <- sort(list.files(bs_pm_folder, pattern = "_bs2.rds",
                                 full.names = TRUE))

## function for getting node connectivity, i.e., sum of |edge strength|
get_node_cnt <- function(x) {
    bs1 <- readRDS(bs1_files[i])
    bs2 <- readRDS(bs2_files[i])
    group_mn <- strsplit(basename(bs1_files[i]), "_bs1.rds")[[1]][1]
    group_m <- strsplit(group_mn, "_vs_")[[1]][1]
    group_n <- strsplit(group_mn, "_vs_")[[1]][2]
    
    y_bs <- list()
    y_bs[1] <- bs1
    y_bs[2] <- bs2

    this_m <- y_bs[[1]]
    this_n <- y_bs[[2]]

    bs_len <- length(this_m)
    cnt <- c()
    for (j in 1 : bs_len) {
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

len <- length(bs1_files)
CNT <- c()
for (i in 1:len) {
    CNT <- rbind(CNT, get_node_cnt(i))
}

out_file <- "CAS_node_connectivity.txt"
write.table(CNT, out_file, quote = F, sep = "\t", row.names = F)

## get the mean CNT
library(dplyr)
library(tidyr)

CNT$Group <- paste0(CNT$Network, "_", CNT$Node)
design <- unique(CNT[, c("Group", "Network", "Node")])
CNT2 <- as.data.frame(CNT %>%
                      group_by(Group) %>%
                      summarise(Mean_CNT = mean(Connectivity)))
CNT2 <- merge(CNT2, design)

## read in taxonomy file
tax <- read.table(tax_file, header = T, sep = "\t")
tax <- tax[rownames(tax) %in% CNT2$Node, c("Phylum", "Class", "Order")]
tax$Node <- rownames(tax)
tax[is.na(tax$Phylum), ]$Phylum <- "Unassigned"

## merge and reformat data frame for plotting
dat <- merge(CNT2, tax)
dat <- dat[, colnames(dat) != "Group"]
dat_wide <- pivot_wider(dat, names_from = Network, values_from = Mean_CNT)
dat_wide <- as.data.frame(dat_wide)

out_file2 <- "CAS_node_connectivity_mean_tax.txt"
write.table(dat_wide, out_file2, quote = F, sep = "\t", row.names = F)
