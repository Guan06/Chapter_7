get_cmp_sum_mean_grp <- function(x, pattern = "ja0", out_dir){

    bs_files <- list.files(x, pattern = paste0("dis_bs_", pattern),
                       full.names = TRUE)
    for (bs in bs_files) {
	    this_bs <- readRDS(bs)
	    cmp <- unlist(strsplit(bs, "dis_bs_"))[2]
	    pm <- paste0(x, "dis_pm_", cmp)
	    this_pm <- readRDS(pm)

	    this_bs[is.na(this_bs)] <- 0
	    this_pm[is.na(this_pm)] <- 0

	    nc <- ncol(this_bs)
	    ## compare the features of each group in bs and pm
	    bs_pm <- c()
	    for (i in 1:nc) {
	        p <- wilcox.test(this_bs[, i], this_pm[, i],
	            alternative = "two.sided")$p.value
	        if (is.na(p)) p <- 1
	        avg_bs <- mean(this_bs[, i])
	        avg_pm <- mean(this_pm[, i])
	        avg_diff <- avg_bs - avg_pm

            accu_bs <- sum(this_bs[, i])
            accu_pm <- sum(this_pm[, i])
            accu_diff <- accu_bs - accu_pm

	        bs_pm <- rbind(bs_pm, c(colnames(this_bs)[i],
                                    accu_bs, accu_pm, accu_diff,
                                    avg_bs, avg_pm, avg_diff, p))
	    }
	    colnames(bs_pm) <- c("Group", "Accu_BS", "Accu_PM", "Accu_Diff",
                             "Avg_BS", "Avg_PM", "Diff", "P_value")
	    bs_pm <- as.data.frame(bs_pm)
	    bs_pm$FDR <- p.adjust(bs_pm$P_value, method = "fdr")

	    bs_pm_sig <- as.data.frame(bs_pm[bs_pm$FDR < 0.05, ])
        bs_pm_sig <- bs_pm_sig[order(bs_pm_sig$Accu_Diff, decreasing = T), ]
	    out_file <- unlist(strsplit(cmp, ".rds"))[1]
	    out_file <- paste0(out_dir, "/", out_file, ".txt")

	    write.table(bs_pm_sig, out_file, quote = F, sep = "\t", row.names = F)
    }
}

add_FC <- function(x, group) {
    dat <- read.table(x, header = T, sep = "\t")
    dat <- merge(dat, group)

    dat$FC <- ifelse(dat$Avg_BS > dat$Avg_PM,
                 (dat$Avg_BS / dat$Avg_PM - 1),
                 (1 - dat$Avg_PM / dat$Avg_BS))
    dat <- dat[order(abs(dat$Diff), decreasing = T), ]
    return(dat)
}
