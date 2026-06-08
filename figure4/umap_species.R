#!/usr/bin/env Rscript
# Figure 4 — Species UMAPs (B.terrestris, A.mellifera, M.pharaonis)
# Same formatting as fig2C umap_atlas.R: manual geom_point + geom_label_repel
# Output: EPS + PNG at 90x60mm per species
# Author: David Bellini

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

options(Seurat.object.assay.version = "v3")

out_dir  <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure4/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 90
panel_height_mm <- 60
mm_to_in        <- 1 / 25.4

save_eps <- function(p, path) {
  postscript(path,
             width  = panel_width_mm * mm_to_in,
             height = panel_height_mm * mm_to_in,
             family = "Helvetica", paper = "special",
             onefile = FALSE, horizontal = FALSE)
  print(p); dev.off()
}

save_png <- function(p, path) {
  ggsave(path, p,
         width = panel_width_mm, height = panel_height_mm,
         units = "mm", dpi = 300, bg = "white")
}

make_umap <- function(so, title) {
  umap_coords <- as.data.frame(so@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$cell_type <- as.character(Idents(so))

  centroids <- umap_coords %>%
    group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

  cluster_levels <- names(sort(table(umap_coords$cell_type), decreasing = TRUE))
  n_clusters     <- length(cluster_levels)
  pal            <- setNames(scales::hue_pal()(n_clusters), cluster_levels)

  ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
    geom_point(size = 0.15, stroke = 0, shape = 16, show.legend = FALSE) +
    geom_label_repel(data          = centroids,
                     aes(label     = cell_type, fill = cell_type),
                     color         = "black",
                     size          = 3 / .pt,
                     fontface      = "plain",
                     family        = "Helvetica",
                     label.padding = unit(0.08, "lines"),
                     label.size    = 0.1,
                     max.overlaps  = Inf,
                     show.legend   = FALSE) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    ggtitle(title) +
    theme_void(base_family = "Helvetica") +
    theme(
      plot.title  = element_text(size = 6, hjust = 0.5, family = "Helvetica"),
      plot.margin = margin(1, 1, 1, 1, "mm")
    )
}

# ── B. terrestris ─────────────────────────────────────────────────────────────
cat("Loading bumblebee...\n")
so_bt <- readRDS("/home/davidbellini/OneDrive/gabbiani/PROJECTS/bumblebee_atlas/bumblebee_AC_v1.rds")
Idents(so_bt) <- "seurat_clusters"
cluster_names_bt <- c(
  "0"  = "PR+",           "1"  = "A/B-KC",        "2"  = "A/B-KC",
  "3"  = "ENS/AST",       "4"  = "N35/N33",        "5"  = "G-KC",
  "6"  = "N55/Tm3",       "7"  = "A/B-KC/Tm9a",   "8"  = "A/B-KC/N102",
  "9"  = "T4",            "10" = "T5",             "11" = "N55",
  "12" = "N69/TmY3",      "13" = "N1",             "14" = "N108",
  "15" = "A/B-KC",        "16" = "fl_123/C3",      "17" = "G-KC/N1",
  "18" = "TmY5a/Tm29",    "19" = "TmY5a/N116",     "20" = "Mi9",
  "21" = "T2/T2a",        "22" = "bt_22",          "23" = "Dm3",
  "24" = "Mi15",          "25" = "Dm3",            "26" = "Tm2",
  "27" = "Mi1/Tm2",       "28" = "PNG",            "29" = "Tm1/Tm2/Tm4",
  "30" = "L1/L2/L3",      "31" = "Dm8/Dm11",       "32" = "Mi4",
  "33" = "N34",           "34" = "C2",             "35" = "N42",
  "36" = "T1",            "37" = "L5",             "38" = "bt_38",
  "39" = "PLM",           "40" = "Tm9/TmY3",       "41" = "bt_41",
  "42" = "Dm12"
)
so_bt <- RenameIdents(so_bt, cluster_names_bt)
p_bt <- make_umap(so_bt, "B. terrestris")
save_eps(p_bt, file.path(out_dir, "figA_bumblebee_umap.eps"))
save_png(p_bt, file.path(out_dir, "figA_bumblebee_umap.png"))
cat("Saved figA_bumblebee_umap\n")
rm(so_bt)

# ── A. mellifera ──────────────────────────────────────────────────────────────
cat("Loading honeybee...\n")
so_hb <- readRDS("/home/davidbellini/OneDrive/gabbiani/PROJECTS/honeybee_atlas/ol_honeybee_leiden1_SAMAP_labeled.rds")
Idents(so_hb) <- "seurat_clusters"
cluster_names_hb <- c(
  "0"  = "Mi1_I",         "1"  = "AST/ENS",       "2"  = "Dm9/R7.8",
  "3"  = "Tm2",           "4"  = "Dm8/Dm11",       "5"  = "am_5",
  "6"  = "fl_26/Dm2/Poxn","7"  = "T4/T5_I",       "8"  = "Mi1_II",
  "9"  = "Tm2",           "10" = "Dm8",            "11" = "OPN",
  "12" = "am_12",         "13" = "T2",             "14" = "Dm3/N35",
  "15" = "T4/T5_II",      "16" = "Dm2",            "17" = "fl_N88",
  "18" = "Tm4/T4/T5",     "19" = "T4-5/TmY3",      "20" = "Dm2",
  "21" = "L1/T1/Lai",     "22" = "L2",             "23" = "am_23",
  "24" = "L4/L3",         "25" = "Tm4/T4/T5_II",  "26" = "am_26",
  "27" = "am_27",         "28" = "Mi1_III",        "29" = "am_29",
  "30" = "L3",            "31" = "Tm2/36",         "32" = "am_32",
  "33" = "am_33",         "34" = "am_34",          "35" = "am_35",
  "36" = "am_36",         "37" = "Mi9",            "38" = "am_38",
  "39" = "am_39",         "40" = "Tm20/N76",       "41" = "am_41",
  "42" = "am_42",         "43" = "am_43",          "44" = "am_44",
  "45" = "am_45"
)
so_hb <- RenameIdents(so_hb, cluster_names_hb)
p_hb <- make_umap(so_hb, "A. mellifera")
save_eps(p_hb, file.path(out_dir, "figB_honeybee_umap.eps"))
save_png(p_hb, file.path(out_dir, "figB_honeybee_umap.png"))
cat("Saved figB_honeybee_umap\n")
rm(so_hb)

# ── M. pharaonis ──────────────────────────────────────────────────────────────
cat("Loading pharaonis...\n")
so_mp <- readRDS("/home/davidbellini/OneDrive/gabbiani/PROJECTS/pharaonis/pharaonis_maleonly_ol_relabeled.rds")
p_mp <- make_umap(so_mp, "M. pharaonis")
save_eps(p_mp, file.path(out_dir, "figC_pharaonis_umap.eps"))
save_png(p_mp, file.path(out_dir, "figC_pharaonis_umap.png"))
cat("Saved figC_pharaonis_umap\n")
rm(so_mp)

cat("\nDone. Outputs in:", out_dir, "\n")
