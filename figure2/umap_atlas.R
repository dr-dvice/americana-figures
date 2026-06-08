#!/usr/bin/env Rscript
# Figure 2C — UMAP of S. americana optic lobe atlas (sa3, apr14 labels)
# Nature Ecology and Evolution formatting

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

options(Seurat.object.assay.version = "v3")

rds_path <- "/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.rds"
out_dir  <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 110
panel_height_mm <- 80
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

so <- readRDS(rds_path)
Idents(so) <- "cell_type"

# Extract coords + labels directly — avoids DimPlot pt.size scaling issues
umap_coords <- as.data.frame(so@reductions$umap@cell.embeddings)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$cell_type <- as.character(Idents(so))
umap_coords$cell_type[umap_coords$cell_type == "Transmedullary_neurons1"] <- "Tm_mix1"
umap_coords$cell_type[umap_coords$cell_type == "Transmedullary_neurons2"] <- "Tm_mix2"

cat("Clusters:", paste(sort(unique(umap_coords$cell_type)), collapse = ", "), "\n")

# Cluster centroids for labels
centroids <- umap_coords %>%
  group_by(cell_type) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

# Match Seurat's default color order: clusters sorted by count descending
cluster_levels <- names(sort(table(umap_coords$cell_type), decreasing = TRUE))
n_clusters     <- length(cluster_levels)
pal            <- setNames(scales::hue_pal()(n_clusters), cluster_levels)

p <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.3, stroke = 0, shape = 16, show.legend = FALSE) +
  geom_label_repel(data       = centroids,
                   aes(label  = cell_type, fill = cell_type),
                   color          = "black",
                   size           = 5 / .pt,
                   fontface       = "bold",
                   family         = "Helvetica",
                   label.padding  = unit(0.08, "lines"),
                   label.size     = 0.2,
                   max.overlaps   = Inf,
                   show.legend    = FALSE) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values  = pal) +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    plot.margin     = margin(1, 1, 1, 1, "mm")
  )

save_eps(p, file.path(out_dir, "fig2C_umap_atlas.eps"))
save_png(p, file.path(out_dir, "fig2C_umap_atlas.png"))
