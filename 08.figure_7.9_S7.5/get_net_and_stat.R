#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)
library(readr)

## for bacteria
asv_file <- "../00.data/Bac_ASV_rep.rds"
design_file <- "../00.data/Bac_design_3809.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

out_dir <- "./spectra_bacteria_bs_pm_99"
if (!dir.exists(out_dir)) dir.create(out_dir)

dat <- new("mina", norm = asv, des = design)
dat <- bs_pm(dat, group = "Compartment", g_size = 80, s_size = 50,
        rm = FALSE, sig = TRUE, bs = 99, pm = 99)
write_rds(dat, "bacteria_net_99.rds.xz", "xz")

### get Spectra distance
dat <- net_dis(dat, method = "spectra", dir = out_dir)
dat_stat <- dis_stat(dat)
write.table(dat_stat, "dis_stat_spectra_bacteria.txt", quote = F,
            sep = "\t", row.names = F)
### get Jaccard distance
dat <- net_dis(dat, method = "Jaccard")
dat_stat <- dis_stat(dat)
write.table(dat_stat, "dis_stat_jaccard_bacteria.txt", quote = F,
            sep = "\t", row.names = F)

## for fungi
asv_file <- "../00.data/Fun_ASV_rep.rds"
design_file <- "../00.data/Fun_design_2232.txt"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

out_dir <- "./spectra_fungi_bs_pm_99"
if (!dir.exists(out_dir)) dir.create(out_dir)

dat <- new("mina", norm = asv, des = design)
dat <- bs_pm(dat, group = "Compartment", g_size = 80, s_size = 50,
        rm = FALSE, sig = TRUE, bs = 99, pm = 99)
write_rds(dat, "fungi_net_99.rds.xz", "xz")

### get Spectra distance
dat <- net_dis(dat, method = "spectra", dir = out_dir)
dat_stat <- dis_stat(dat)
write.table(dat_stat, "dis_stat_spectra_fungi.txt", quote = F,
            sep = "\t", row.names = F)
### get Jaccard distance
dat <- net_dis(dat, method = "Jaccard")
dat_stat <- dis_stat(dat)
write.table(dat_stat, "dis_stat_jaccard_fungi.txt", quote = F,
            sep = "\t", row.names = F)
