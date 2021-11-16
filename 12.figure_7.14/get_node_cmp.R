#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)

net_node_cmp("./CAS_vs_Root/", cmp = "contrast", dir = "./node_cmp/")
net_node_cmp("./CAS_vs_Root/", cmp = "ja", dir = "./node_cmp/")
net_node_cmp("./CAS_vs_Root/", cmp = "ja0", dir = "./node_cmp/")

if (!dir.exists("./node_cmp_stat/")) dir.create("./node_cmp_stat/")

get_cmp_stat <- function(x, cmp = "contrast") {

    bs_files <- list.files(x, pattern = paste0("dis_bs_", cmp, "_"),
                           full.names = TRUE)

    for (bs in bs_files){
        this_bs <- readRDS(bs)
        vs <- unlist(strsplit(bs, paste0(cmp, "_")))[2]

        pm <- paste0(x, "dis_pm_", cmp, "_", vs)
        this_pm <- readRDS(pm)

        this_bs[is.na(this_bs)] <- 0
        this_pm[is.na(this_pm)] <- 0

        nc <- ncol(this_bs)

        ## compare the contrast/distance of each node in bs and pm
        bs_pm <- c()
        for (i in 1:nc) {
#            if (sum(this_bs[, i]) == 0) next
#            if (sum(this_pm[, i]) == 0) next
            p <- wilcox.test(this_bs[, i], this_pm[, i],
                alternative = "two.sided")$p.value
            if (is.na(p)) p <- 1
            avg_bs <- mean(this_bs[, i])
            avg_pm <- mean(this_pm[, i])
            diff <- avg_bs - avg_pm
            bs_pm <- rbind(bs_pm, c(colnames(this_bs)[i], avg_bs, avg_pm, diff, p))
        }
        colnames(bs_pm) <- c("Node", "Avg_BS", "Avg_PM", "Diff", "P_value")
        bs_pm <- as.data.frame(bs_pm)
        bs_pm$FDR <- p.adjust(bs_pm$P_value, method = "fdr")

        bs_pm_sig <- bs_pm[bs_pm$FDR < 0.05, ]
        suffix <- unlist(strsplit(vs, ".rds"))[1]
        saveRDS(bs_pm_sig, paste0("./node_cmp_stat/", suffix, "_", cmp, ".rds"))
    }
}

cmp1 <- get_cmp_stat("./node_cmp/", cmp = "contrast")
cmp2 <- get_cmp_stat("./node_cmp/", cmp = "ja")
cmp3 <- get_cmp_stat("./node_cmp/", cmp = "ja0")
