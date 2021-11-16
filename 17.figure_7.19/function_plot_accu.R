plot_accu <- function(x, dis_folder, out_dir, prefix) {
    this_cmp_file <- paste0(dis_folder, "accu_dis_", x, ".txt")
    this_cmp <- read.table(this_cmp_file, header = T)

    p1 <- ggplot(this_cmp) +
            geom_boxplot(aes(Family, Distance), color = "gray64", alpha = 0.4,
                         outlier.shape = NA) +
            geom_point(aes(Family, Mean), color = "navyblue") +
            geom_path(aes(Family, Mean), color = "navyblue", group = 1) +
            scale_x_discrete(limits = x_order) +
            main_theme +
            theme(plot.background = element_blank(),
                  plot.title = element_text(size = 8),
                axis.text.x = element_text(colour = "black", angle = 90,
                size = 5, hjust = 1))

    this_mean_file <- paste0(dis_folder, "accu_dis_", x, "_mean.txt")
    this_mean <- read.table(this_mean_file, header = T)
    this_mean <- this_mean[match(x_order, this_mean$Family), ]

    ### plot dual for sig
    p2 <- ggplot(this_mean, aes(x = Family)) +
        geom_point(aes(y = P_adjust_W), color = "#00ADDA", alpha = 0.7) +
        geom_path(aes(y = P_adjust_W), color = "#00ADDA", group = 1, alpha = 0.7) +
        geom_point(aes(y = P_value_F), color = "gold3") +
        geom_path(aes(y = P_value_F), color = "gold3", group = 1) +
        geom_point(aes(y = P_adjust_T), color = "deeppink", alpha = 0.7) +
        geom_path(aes(y = P_adjust_T), color = "deeppink", group = 1, alpha = 0.7) +
        scale_x_discrete(limits = x_order) +
        scale_y_continuous(position = "right") +
        labs(y = "P value") +
        geom_hline(yintercept = 0.05, color = "salmon", linetype = "dashed") +
        geom_hline(yintercept = 0.5, color = "gray", alpha = 0.7) +
        geom_hline(yintercept = 0.4, color = "gray", alpha = 0.7,
                   linetype = "dashed") +
        geom_hline(yintercept = 0.6, color = "gray", alpha = 0.7,
                   linetype = "dashed") +
        main_theme +
        theme(plot.background = element_blank(),
              legend.position = "none",
              plot.title = element_text(size = 8),
                axis.text.x = element_text(colour = "black", angle = 90,
                size = 5, hjust = 1))

    aligned_plots <- align_plots(p1, p2, align = "hv", axis = "tblr")
    p <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])

    plot_file <- paste0(out_dir, prefix, "_", x, ".pdf")
    ggsave(plot_file, p, height = 4, width = 10)
}


### function for script_figure_7_19.R

plot_accu_indi <- function(x, dis_folder) {
    this_cmp_file <- paste0(dis_folder, "accu_dis_", x, ".txt")
    this_cmp <- read.table(this_cmp_file, header = T)

    p1 <- ggplot(this_cmp) +
            geom_boxplot(aes(Family, Distance), color = "gray64", alpha = 0.4,
                         outlier.shape = NA) +
            geom_point(aes(Family, Mean), color = "navyblue") +
            geom_path(aes(Family, Mean), color = "navyblue", group = 1) +
            scale_x_discrete(limits = x_order) +
            labs(x = "") +
            main_theme +
            theme(plot.background = element_blank(),
                  plot.title = element_text(size = 8),
                axis.text.x = element_text(colour = "black", angle = 90,
                size = 5, hjust = 1))

    this_mean_file <- paste0(dis_folder, "accu_dis_", x, "_mean.txt")
    this_mean <- read.table(this_mean_file, header = T)
    this_mean <- this_mean[match(x_order, this_mean$Family), ]

    ### plot dual for sig
    p2 <- ggplot(this_mean, aes(x = Family)) +
        geom_point(aes(y = P_adjust_W), color = "#00ADDA", alpha = 0.8) +
        geom_path(aes(y = P_adjust_W), color = "#00ADDA", group = 1, alpha = 0.8) +
        geom_point(aes(y = P_value_F), color = "gold3", alpha = 0.8) +
        geom_path(aes(y = P_value_F), color = "gold3", group = 1, alpha = 0.8) +
#        geom_point(aes(y = P_adjust_T), color = "deeppink", alpha = 0.7) +
#        geom_path(aes(y = P_adjust_T), color = "deeppink", group = 1, alpha = 0.7) +
        scale_x_discrete(limits = x_order) +
        scale_y_continuous(position = "right") +
        labs(y = "P value", x = "") +
        geom_hline(yintercept = 0.05, color = "salmon", linetype = "dashed") +
        geom_hline(yintercept = 0.5, color = "gray", alpha = 0.7) +
        geom_hline(yintercept = 0.4, color = "gray", alpha = 0.7,
                   linetype = "dashed") +
        geom_hline(yintercept = 0.6, color = "gray", alpha = 0.7,
                   linetype = "dashed") +
        main_theme +
        theme(plot.background = element_blank(),
              legend.position = "none",
              plot.title = element_text(size = 8),
                axis.text.x = element_text(colour = "black", angle = 90,
                size = 5, hjust = 1))

    aligned_plots <- align_plots(p1, p2, align = "hv", axis = "tblr")
    p <- ggdraw(aligned_plots[[2]]) + draw_plot(aligned_plots[[1]])
}
