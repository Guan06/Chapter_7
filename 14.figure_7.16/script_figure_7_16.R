#!/biodata/dep_psl/grp_rgo/guan/bin//Rscript

library(dplyr)
source("../01.common_scripts/plot_setting.R")

### Before plotting ###########################################################
# Run the scripts get_grp_feature.R and get_network_cluster_aRA.R in the same
# folder before plotting.
###############################################################################

### function for plotting
plot_feature <- function(x, ra) {
    dat <- read.table(x, header = T, sep = "\t")
    dat <- merge(dat, ra)
    dat$Color <- as.character(dat$Color)

    lim <- (round(max(abs(dat$FC)) / 0.1) + 1) * 0.1
#    p <- ggplot(dat, aes(x = FC, abs(Accu_Diff))) +
    p <- ggplot(dat, aes(x = FC,  avg_RA)) +
        geom_point(aes(size = abs(Diff), color = Color), shape = 16, alpha = 0.7) +
        main_theme +
        xlim(-lim, lim) +
        scale_y_log10() +
        scale_size(limits = c(0, 0.5), range = c(0, 3)) +
        theme(legend.position = "top", axis.line.y = element_blank()) +
        labs(x = "Fold Change", y = "aRA of the group") +
        geom_vline(xintercept = 0, color = "gray") +
        geom_hline(yintercept = 0.005, color = "salmon", linetype = "dashed") +
        geom_vline(xintercept = 1, color = "gold", linetype = "dashed") +
        geom_vline(xintercept = -1, color = "gold", linetype = "dashed")

    return(p)
}
###

## for panel a, b, family comparison
fam_file <- "../13.figure_7.15/aRA_occupancy_families.txt"
fam <- read.table(fam_file, header = T, sep = "\t")

fam_ra <- fam %>% group_by(Taxonomy) %>% summarise(avg_RA = mean(RA))
fam_ra <- as.data.frame(fam_ra[fam_ra$Taxonomy != "Unassigned", ])
colnames(fam_ra)[1] <- "Group"

p1 <- plot_feature("./Family_CAS_vs_At_root.txt", ra = fam_ra)
p1 <- p1 + scale_color_manual(values = c_phyla, guide = "none")

p2 <- plot_feature("./Family_CAS_vs_Lj_root.txt", ra = fam_ra)
p2 <- p2 + scale_color_manual(values = c_phyla, guide = "none")

## for panel c, d, network cluster comparison
nc_file <- "aRA_occupancy_network_clusters.txt"
nc <- read.table(nc_file, header = T, sep = "\t")
nc_ra <- nc %>% group_by(Group) %>% summarise(avg_RA = mean(RA))
nc_ra <- as.data.frame(nc_ra)

## get avg_RA for each network cluster
cc <- length(nc_ra$Group)

p3 <- plot_feature("NetCls_CAS_vs_At_root.txt", ra = nc_ra)
p3 <- p3 + scale_color_manual(values = getDark2(cc), guide = "none")

## remove the outlier
dat <- read.table("NetCls_CAS_vs_At_root.txt", header = T, sep = "\t")
dat <- merge(dat, nc_ra)
dat$Color <- as.character(dat$Color)
lim <- 2
p3_2 <- ggplot(dat, aes(x = FC,  avg_RA)) +
        geom_point(aes(size = abs(Diff), color = Color), shape = 16, alpha = 0.7) +
        main_theme +
        xlim(-lim, lim) +
        scale_y_log10() +
        scale_size(limits = c(0, 0.5), range = c(0, 3)) +
        theme(legend.position = "top", axis.line.y = element_blank()) +
        labs(x = "Fold Change", y = "aRA of the network cluster") +
        geom_vline(xintercept = 0, color = "gray") +
        geom_hline(yintercept = 0.005, color = "salmon", linetype = "dashed") +
        geom_vline(xintercept = 1, color = "gold", linetype = "dashed") +
        geom_vline(xintercept = -1, color = "gold", linetype = "dashed") +
        scale_color_manual(values = getDark2(cc), guide = "none")

p4 <- plot_feature("NetCls_CAS_vs_Lj_root.txt", ra= nc_ra)
p4 <- p4 + scale_color_manual(values = getDark2(cc), guide = "none")

p <- plot_grid(p1, p2, p3_2, p4, nrow = 2, labels = 'auto')
ggsave("Figure_7_16.pdf", p, width = 10)
ggsave("Figure_S7_15.pdf", p3, width = 5, height = 3.5)
