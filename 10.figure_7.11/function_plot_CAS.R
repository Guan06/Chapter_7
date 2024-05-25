plot_CAS <- function(x) {
    glst <- unique(des$Group)

    occu_ra_df <- c()
    for (g in glst){
        this_des <- des[des$Group == g, ]
        this_x <- x[, colnames(x) %in% this_des$Sample_ID]
        this_x <- this_x[rowSums(this_x) > 0, ]

        this_com <- unique(this_des$Compartment)
        this_host <- ifelse(g == "soil", "host_soil",
                            unique(this_des$Host.Species))

        this <- data.frame(Taxon = rownames(this_x),
                           Occupancy = rowSums(this_x > 0) / ncol(this_x),
                           RA = rowSums(this_x) / ncol(this_x),
                           Group = g,
                           Compartment = this_com,
                           Host = this_host)

        this <- this[order(this$Occupancy), ]
        occu_ra_df <- rbind(occu_ra_df, this)
    }
    p <- plot_dual(occu_ra_df)
}

plot_dual <- function(x) {
    lt_host <-c("arabidopsis_thaliana" = "solid",
                 "chlamydomonas_reinhardtii" = "solid",
                 "lotus_japonicus" = "twodash", "host_soil" = "solid")

     p1 <- ggplot(x, aes(x = Occupancy,
                        color = Compartment, linetype = Host)) +
        geom_density() +
        main_theme +
        scale_color_manual(values = c_Com) +
        scale_linetype_manual(values = lt_host) +
        scale_fill_manual(values = c_Com) +
        labs(x = "Occupancy", y = "Density of taxa") +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.line.x = element_blank(),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 6),
              plot.title = element_text(size = 6))

    lt_host <- c("arabidopsis_thaliana" = "dashed",
                 "chlamydomonas_reinhardtii" = "dashed",
                 "lotus_japonicus" = "dotted", "host_soil" = "dashed")

    p2 <- ggplot(x %>% group_by(Group) %>% mutate(cv = cumsum(RA)),
            aes(x = Occupancy, y = cv,
                colour = factor(Compartment), linetype = Host)) +
        geom_line() +
        scale_y_continuous(limits = c(0,1), position = "right") +
        main_theme +
        scale_color_manual(values = c_Com) +
        scale_linetype_manual(values = lt_host) +
        theme_minimal_grid(color = "gray88", line_size = 0.3) +
        labs(y = "aRA of taxa") +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.line.x = element_line(color = "black"),
              axis.line.y = element_line(color = "black"),
              axis.text.y = element_text(size = 6),
              axis.title.y = element_text(size = 6),
              axis.title.x = element_blank())

    aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
    ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
}
