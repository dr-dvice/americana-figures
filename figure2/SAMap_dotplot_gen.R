#!/usr/bin/env Rscript
# Figure 2B — SAMap alignment dotplots, S. americana vs 3 fly atlases

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
})

options(Seurat.object.assay.version = "v3")

samap_dir  <- "/home/davidbellini/files/workspace/americana_atlas3/results/samap/comparisons"
output_dir <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure2/fig2b_candidates_aug17"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 110
panel_height_mm <- 60
mm_to_in        <- 1 / 25.4

save_eps <- function(p, path) {
  postscript(path,
             width  = panel_width_mm  * mm_to_in,
             height = panel_height_mm * mm_to_in,
             family = "Helvetica", paper = "special",
             onefile = FALSE, horizontal = FALSE)
  print(p); dev.off()
}

save_png <- function(p, path) {
  ggsave(path, p,
         width = panel_width_mm, height = panel_height_mm,
         units = "mm", dpi = 300)
}

# Load Seurat object for cluster color mapping (colors keyed by cluster number)
so <- readRDS("/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_aug17.rds")
cluster_ids    <- sort(as.numeric(as.character(unique(so@meta.data$SCT_snn_res.1))))
cluster_colors <- hue_pal()(n = length(cluster_ids))
names(cluster_colors) <- as.character(cluster_ids)

# fly2020 = Kurmangaliyev 2020, fly2018 = Davie 2018, fly2021 = Ozel 2021
datasets <- list(
  list(csv = "sa3_fly2020_mapping_table_ntop100.csv", out = "fig2B_dotplot_kurmangaliyev2020_ntop100"),
  list(csv = "sa3_fly2020_mapping_table_all.csv",     out = "fig2B_dotplot_kurmangaliyev2020_all"),
  list(csv = "sa3_fly2018_mapping_table_ntop100.csv", out = "fig2B_dotplot_davie2018_ntop100"),
  list(csv = "sa3_fly2018_mapping_table_all.csv",     out = "fig2B_dotplot_davie2018_all"),
  list(csv = "sa3_fly2021_mapping_table_ntop100.csv", out = "fig2B_dotplot_ozel2021_ntop100"),
  list(csv = "sa3_fly2021_mapping_table_all.csv",     out = "fig2B_dotplot_ozel2021_all")
)

for (ds in datasets) {
  in_file <- file.path(samap_dir, ds$csv)
  cat("Processing:", in_file, "\n")
  table_format <- read_csv(in_file, show_col_types = FALSE)

  table_long <- table_format %>%
    pivot_longer(cols = -`...1`, names_to = "Fly_Cluster", values_to = "Alignment_Score") %>%
    rename(SA_Cluster = `...1`) %>%
    filter(str_detect(SA_Cluster, "sa"), str_detect(Fly_Cluster, "fl")) %>%
    mutate(SA_Cluster  = gsub("sa_", "", SA_Cluster),
           Fly_Cluster = gsub("fl_", "", Fly_Cluster))

  table_long_filter <- table_long %>%
    filter(Alignment_Score > 0.3) %>%
    mutate(SA_Cluster = as.numeric(SA_Cluster))
  custom_order <- sort(unique(table_long_filter$SA_Cluster))

  table_color <- table_long_filter %>%
    mutate(SA_Cluster = as.character(SA_Cluster)) %>%
    group_by(SA_Cluster) %>%
    arrange(desc(Alignment_Score)) %>%
    mutate(rank = row_number(), score_bin = floor(Alignment_Score * 5) / 5) %>%
    group_by(SA_Cluster, score_bin) %>%
    mutate(bin_rank = row_number()) %>%
    ungroup() %>%
    group_by(SA_Cluster) %>%
    filter(rank <= 3 | bin_rank == 1) %>%
    ungroup() %>%
    mutate(Dot_Color = ifelse(rank <= 2, cluster_colors[SA_Cluster], "#BEBEBE"))

  p <- ggplot(table_color,
              aes(x = factor(SA_Cluster, levels = custom_order),
                  y = Alignment_Score)) +
    geom_text_repel(aes(label = Fly_Cluster),
                    size            = 3 / .pt,
                    color           = "black",
                    max.overlaps    = Inf,
                    box.padding     = 0.2,
                    min.segment.length = 0,
                    family          = "Helvetica") +
    geom_point(aes(fill = Dot_Color), shape = 21, size = 1.5,
               color = "black", stroke = 0.3, show.legend = FALSE) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = expression(italic("S. americana") ~ "cluster"),
         y = "Alignment score") +
    theme_minimal(base_family = "Helvetica") +
    theme(
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 5, family = "Helvetica"),
      axis.text.y        = element_text(size = 5, family = "Helvetica"),
      axis.title.x       = element_text(size = 6, family = "Helvetica"),
      axis.title.y       = element_text(size = 6, family = "Helvetica"),
      legend.position    = "right",
      legend.text        = element_text(size = 5, family = "Helvetica"),
      legend.key.size    = unit(3, "mm"),
      legend.spacing.y   = unit(0.5, "mm"),
      plot.margin        = margin(1, 1, 1, 1, "mm")
    )

  save_eps(p, file.path(output_dir, paste0(ds$out, ".eps")))
  save_png(p, file.path(output_dir, paste0(ds$out, ".png")))
  cat("Done:", ds$out, "\n")
}
