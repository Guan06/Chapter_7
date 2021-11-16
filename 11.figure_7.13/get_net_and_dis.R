#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(mina)

asv_file <- "../00.data/CAS_rep_ASV_reNorm.rds"
design_file <- "../00.data/design_854.txt"
bs_pm_dir <- "./CAS_bs_pm/"
out_dir <- "./CAS_spectra_egv/"

asv <- readRDS(asv_file)
design <- read.table(design_file, header = T, sep = "\t")

des_soil <- design[design$Compartment == "soil", ]
des_soil$Group <- des_soil$Soil

des_nsoil <- design[design$Compartment != "soil", ]
des_nsoil$Group <- paste0(des_nsoil$Host.Species, "_", des_nsoil$Compartment)
des <- rbind(des_soil, des_nsoil)

dat <- new("mina", norm = asv, des = des)

if (!dir.exists(bs_pm_dir)) dir.create(bs_pm_dir)
if (!dir.exists(out_dir)) dir.create(out_dir)

## for writing the correlation matrices of each subsampling and network
## inference in both origianl and permuation datasets

#dat <- bs_pm(dat, group = "Group", g_size = 40, s_size = 20,
#            rm = FALSE, sig = TRUE, bs = 33, pm = 33,
#            individual = TRUE, out_dir = bs_pm_dir)

## for the stat calculation
dat <- bs_pm(dat, group = "Group", g_size = 40, s_size = 20,
            rm = FALSE, sig = TRUE, bs = 33, pm = 33,
            individual = FALSE)

dat <- net_dis(dat, method = "spectra", egv = TRUE, dir = out_dir)
dis1 <- dis_stat(dat)
write.table(dis1, "net_dis_spectra.txt", quote = F, sep = "\t", row.names = F)

dat <- net_dis(dat, method = "Jaccard")
dis2 <- dis_stat(dat)
write.table(dis2, "net_dis_jaccard.txt", quote = F, sep = "\t", row.names = F)
