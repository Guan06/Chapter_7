#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dplyr)
asv_file <- "../00.data/Bac_ASV_raref.rds"
rep_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/design_854.txt"
group_file <- "../00.data/CAS_spearman_ap_tab.txt"

asv <- readRDS(asv_file)
rep <- readRDS(rep_file)
design <- read.table(design_file, header = T, sep = "\t")
des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- "soil"

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)

des <- rbind(des_soil, des_nsoil)
g_lst <- unique(des$Group)

group <- read.table(group_file, header = T, sep = "\t")
group$Group <- group$Cluster <- as.character(group$Cluster)

df_ra <- c()

for (g in g_lst) {
    this_des <- des[des$Group == g, ]
    this_asv <- asv[rownames(asv) %in% rownames(rep),
                    colnames(asv) %in% this_des$Sample_ID]
    this_asv <- this_asv[rowSums(this_asv) > 0, colSums(this_asv) > 0]

    this_fill <- unique(this_des$Compartment)
    this_grp <- group[group$ID %in% rownames(this_asv), ]

    this_asv <- this_asv[match(this_grp$ID, rownames(this_asv)), ]

    GRP <- as.factor(this_grp$Group)
    tab_GRP <- rowsum(this_asv, GRP)
#    tab_GRP <- apply(this_asv, 2, function(x) rowsum(x, GRP))
#    rownames(tab_GRP) <- levels(GRP)

	this <- data.frame(Group = rownames(tab_GRP),
	                    Occupancy = rowSums(tab_GRP > 0) / ncol(tab_GRP),
	                    Sample_number = ncol(tab_GRP),
	                    RA = rowSums(tab_GRP) / ncol(tab_GRP),
	                    Group = g,
                        Compartment = this_fill,
	                    Rank = "NetCls",
	                    Tax_number = nrow(tab_GRP))
    df_ra <- rbind(df_ra, this)
}
write.table(df_ra, "aRA_occupancy_network_clusters.txt", quote = F, sep = "\t",
            row.names = F)
