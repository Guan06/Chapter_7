#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

source("../01.common_scripts/order.R")
source("../01.common_scripts/plot_setting.R")

bac_file <- "Bacteria_net_features.txt"
fun_file <- "Fungi_net_features.txt"
bac <- read.table(bac_file, header = T, sep = "\t")
fun <- read.table(fun_file, header = T, sep = "\t")

info_lst <- c("Group", "Soil", "Host", "Compartment")

bac$Neg_density <- -bac$Neg_density
bac_pos <- bac[, c(info_lst, "Pos_density")]
bac_neg <- bac[, c(info_lst, "Neg_density")]

colnames(bac_pos)[5] <- colnames(bac_neg)[5] <- "Density"
bac_both <- rbind(bac_pos, bac_neg)
bac_both$Kingdom <- "Bacteria"

fun$Neg_density <- -fun$Neg_density
fun_pos <- fun[, c(info_lst, "Pos_density")]
fun_neg <- fun[, c(info_lst, "Neg_density")]

colnames(fun_pos)[5] <- colnames(fun_neg)[5] <- "Density"
fun_both <- rbind(fun_pos, fun_neg)
fun_both$Kingdom <- "Fungi"

all <- rbind(bac_both, fun_both)

p0 <-ggplot(all, aes(Group, Density, fill = Compartment, color = Kingdom)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    main_theme + theme(legend.position = "top") +
    scale_fill_manual(values = c_Com) +
    scale_color_manual(values = c("NA", "black")) +
    theme(axis.title.x = element_blank(), axis.line.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 0, size = 0.5)

p2 <- p0 + scale_x_discrete(limits = order_3a, labels = order_3b)

ggsave("Figure_7_6.pdf", p2, width = 10, height = 5)
