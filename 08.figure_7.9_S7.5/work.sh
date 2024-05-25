cp /biodata/dep_psl/grp_rgo/guan/20200616_figures/19.MANUSCRIPT/23.network_comparison_whole_dataset_compartment/stat_spectra_bacteria.txt dis_stat_spectra_bacteria_33.txt
cp /biodata/dep_psl/grp_rgo/guan/20200616_figures/19.MANUSCRIPT/23.network_comparison_whole_dataset_compartment/stat_spectra_fungi_v2.txt dis_stat_spectra_fungi_33.txt

Rscript  script_figure_7_9_and_S7_5.R  ./spectra_bacteria_bs_pm_33/ Figure_7_9_bacteria.pdf
Rscript script_figure_7_9_and_S7_5.R ./spectra_fungi_bs_pm_33/ Figure_7_9_fungi.pdf
Rscript script_figure_7_9_and_S7_5.R  ./spectra_bacteria_bs_pm_99/ Figure_S7_5_bacteria.pdf
Rscript script_figure_7_9_and_S7_5.R  ./spectra_fungi_bs_pm_99/ Figure_S7_5_fungi.pdf
