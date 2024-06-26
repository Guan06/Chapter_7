plot_order <- c("zea_mays_l_DEMO_soil",
                "zea_mays_l_DEMO_rhizosphere",
                "zea_mays_l_DEMO_root",
                "zea_mays_l_DOK_soil",
                "zea_mays_l_DOK_rhizosphere",
                "zea_mays_l_DOK_root",
                "lotus_japonicus_CAS_soil",
                "lotus_japonicus_CAS_rhizosphere",
                "lotus_japonicus_CAS_root",
                "chlamydomonas_reinhardtii_CAS_phycosphere",
                "arabidopsis_thaliana_CAS_soil",
                "arabidopsis_thaliana_CAS_rhizosphere",
                "arabidopsis_thaliana_CAS_rhizoplane",
                "arabidopsis_thaliana_CAS_root",
                "arabidopsis_thaliana_HAS_soil",
                "arabidopsis_thaliana_HAS_rhizosphere",
                "arabidopsis_thaliana_HAS_root",
                "arabidopsis_thaliana_Germany_soil",
                "arabidopsis_thaliana_Germany_rhizosphere",
                "arabidopsis_thaliana_Germany_rhizoplane",
                "arabidopsis_thaliana_Germany_root",
                "arabidopsis_thaliana_France_soil",
                "arabidopsis_thaliana_France_rhizosphere",
                "arabidopsis_thaliana_France_rhizoplane",
                "arabidopsis_thaliana_France_root",
                "arabidopsis_thaliana_Spain_soil",
                "arabidopsis_thaliana_Spain_rhizosphere",
                "arabidopsis_thaliana_Spain_rhizoplane",
                "arabidopsis_thaliana_Spain_root",
                "arabidopsis_thaliana_Sweden_soil",
                "arabidopsis_thaliana_Sweden_rhizosphere",
                "arabidopsis_thaliana_Sweden_rhizoplane",
                "arabidopsis_thaliana_Sweden_root",
                "arabidopsis_thaliana_Italy_soil",
                "arabidopsis_thaliana_Italy_rhizosphere",
                "arabidopsis_thaliana_Italy_rhizoplane",
                "arabidopsis_thaliana_Italy_root")
bac_labels <- c("Zm_DEMO_soil", "Zm_DEMO_rhizos", "Zm_DEMO_root",
                "Zm_DOK_soil", "Zm_DOK_rhizos", "Zm_DOK_root",
                "Lj_CAS_soil", "Lj_CAS_rhizos", "Lj_CAS_root",
                "Chlamy_phyco",
                "At_CAS_soil", "At_CAS_rhizos", "At_CAS_rhizop", "At_CAS_root",
                "At_HAS_soil", "At_HAS_rhizos", "At_HAS_root",
                "At_DE_soil", "At_DE_rhizos", "At_DE_rhizop", "At_DE_root",
                "At_FR_soil", "At_FR_rhizos", "At_FR_rhizop", "At_FR_root",
                "At_ES_soil", "At_ES_rhizos", "At_ES_rhizop", "At_ES_root",
                "At_SE_soil", "At_SE_rhizos", "At_SE_rhizop", "At_SE_root",
                "At_IT_soil", "At_IT_rhizos", "At_IT_rhizop", "At_IT_root")

fun_plot_order <- c("zea_mays_l_DEMO_soil",
                "zea_mays_l_DEMO_rhizosphere",
                "zea_mays_l_DEMO_root",
                "zea_mays_l_DOK_soil",
                "zea_mays_l_DOK_rhizosphere",
                "zea_mays_l_DOK_root",
                "lotus_japonicus_CAS_soil",
                "lotus_japonicus_CAS_rhizosphere",
                "lotus_japonicus_CAS_root",
                "arabidopsis_thaliana_Germany_soil",
                "arabidopsis_thaliana_Germany_rhizosphere",
                "arabidopsis_thaliana_Germany_rhizoplane",
                "arabidopsis_thaliana_Germany_root",
                "arabidopsis_thaliana_France_soil",
                "arabidopsis_thaliana_France_rhizosphere",
                "arabidopsis_thaliana_France_rhizoplane",
                "arabidopsis_thaliana_France_root",
                "arabidopsis_thaliana_Spain_soil",
                "arabidopsis_thaliana_Spain_rhizosphere",
                "arabidopsis_thaliana_Spain_rhizoplane",
                "arabidopsis_thaliana_Spain_root",
                "arabidopsis_thaliana_Sweden_soil",
                "arabidopsis_thaliana_Sweden_rhizosphere",
                "arabidopsis_thaliana_Sweden_rhizoplane",
                "arabidopsis_thaliana_Sweden_root",
                "arabidopsis_thaliana_Italy_soil",
                "arabidopsis_thaliana_Italy_rhizosphere",
                "arabidopsis_thaliana_Italy_rhizoplane",
                "arabidopsis_thaliana_Italy_root")

fun_labels <- c("Zm_DEMO_soil", "Zm_DEMO_rhizos", "Zm_DEMO_root",
                "Zm_DOK_soil", "Zm_DOK_rhizos", "Zm_DOK_root",
                "Lj_CAS_soil", "Lj_CAS_rhizos", "Lj_CAS_root",
                "At_DE_soil", "At_DE_rhizos", "At_DE_rhizop", "At_DE_root",
                "At_FR_soil", "At_FR_rhizos", "At_FR_rhizop", "At_FR_root",
                "At_ES_soil", "At_ES_rhizos", "At_ES_rhizop", "At_ES_root",
                "At_SE_soil", "At_SE_rhizos", "At_SE_rhizop", "At_SE_root",
                "At_IT_soil", "At_IT_rhizos", "At_IT_rhizop", "At_IT_root")

order_3a <- c("DEMO",
                "zea_mays_l_DEMO_rhizosphere",
                "zea_mays_l_DEMO_root",
                "DOK",
                "zea_mays_l_DOK_rhizosphere",
                "zea_mays_l_DOK_root",
                "CAS",
                "lotus_japonicus_CAS_rhizosphere",
                "lotus_japonicus_CAS_root",
                "chlamydomonas_reinhardtii_CAS_phycosphere",
                "arabidopsis_thaliana_CAS_rhizosphere",
                "arabidopsis_thaliana_CAS_rhizoplane",
                "arabidopsis_thaliana_CAS_root",
                "HAS",
                "arabidopsis_thaliana_HAS_rhizosphere",
                "arabidopsis_thaliana_HAS_root",
                "Germany",
                "arabidopsis_thaliana_Germany_rhizosphere",
                "arabidopsis_thaliana_Germany_rhizoplane",
                "arabidopsis_thaliana_Germany_root",
                "France",
                "Spain",
                "arabidopsis_thaliana_Spain_rhizosphere",
                "arabidopsis_thaliana_Spain_root",
                "Sweden",
                "arabidopsis_thaliana_Sweden_rhizosphere",
                "arabidopsis_thaliana_Sweden_root",
                "Italy",
                "arabidopsis_thaliana_Italy_root")

order_3b <- c("DEMO", "Zm_DEMO_rhizo", "Zm_DEMO_root",
               "DOK", "Zm_DOK_rhizo", "Zm_DOK_root",
                 "CAS", "Lj_CAS_rhizo", "Lj_CAS_root", "Chlamy_CAS_phyco",
                 "At_CAS_rhizo", "At_CAS_rhizop", "At_CAS_root",
                 "HAS", "At_HAS_rhizo", "At_HAS_root",
                 "Germany", "At_DE_rhizo", "At_DE_rhizop", "At_DE_root",
                 "France", "Spain", "At_ES_rhizo", "At_ES_root",
                 "Sweden", "At_SE_rhizo", "At_SE_root",
                 "Italy", "At_IT_root")
