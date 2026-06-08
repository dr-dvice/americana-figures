# Author: David Bellini
# Honeybee OL atlas (leiden1 SAMap-labeled) — cluster inspection & relabeling helper
# Run interactively: /opt/R/4.2.3/bin/R
library(tidyverse)
library(Seurat)
library(patchwork)

atlas_dir  <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/honeybee_atlas"
homolog_fp <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/all_atlases/SCENIC/TF_gene_find/honeybee_to_fly_homologs.tsv"
rds_fp     <- file.path(atlas_dir, "ol_honeybee_leiden1_SAMAP_labeled.rds")
out_rds    <- file.path(atlas_dir, "ol_honeybee_leiden1_relabeled.rds")

so <- readRDS(rds_fp)

cat("Assays:", paste(Assays(so), collapse = ", "), "\n")
cat("Default assay:", DefaultAssay(so), "\n")
cat("Cells:", ncol(so), " Genes:", nrow(so), "\n")
cat("Reductions:", paste(Reductions(so), collapse = ", "), "\n")
cat("Meta cols:", paste(colnames(so@meta.data), collapse = ", "), "\n")
print(table(Idents(so)))

homologs <- read_tsv(homolog_fp, show_col_types = FALSE)

# ---- Optic-lobe TF panel (fly names) -> honeybee LOC IDs ----
ol_tfs_fly <- c(
  "dac","hth","tj","vvl","bsh","D","toy","oc","dfr","Optix",
  "ato","so","eya","sv","drf","fkh","ey","Lim1","Lim3","acj6",
  "tll","sim","scro","exex","run","slp1","slp2","ap","disco",
  "br","pdm3","Hr51","Eip93F","fd59A","onecut","Ets65A","CG32532",
  "TfAP-2"
)

map_fly_to_bee <- function(fly_genes) {
  homologs %>%
    filter(fly_gene %in% fly_genes) %>%
    group_by(fly_gene) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(fly_gene, bee_gene = sp_transcript, product = fly_product, pident, bitscore)
}

tf_map <- map_fly_to_bee(ol_tfs_fly)
tf_map$in_object <- tf_map$bee_gene %in% rownames(so)
cat("\nTF mapping (fly -> honeybee LOC):\n"); print(tf_map, n = Inf)

tf_features <- tf_map %>% filter(in_object) %>% pull(bee_gene)
names(tf_features) <- tf_map %>% filter(in_object) %>% pull(fly_gene)

# ---- Convenience plotters ----
plot_tfs_feature <- function(features = tf_features, reduction = "umap", ncol = 4) {
  ps <- lapply(seq_along(features), function(i) {
    FeaturePlot(so, features = features[i], reduction = reduction, order = TRUE) +
      ggtitle(sprintf("%s (%s)", names(features)[i], features[i]))
  })
  wrap_plots(ps, ncol = ncol)
}

# Save FeaturePlots in batches (default 4 per PNG). group_size = 1 -> one TF per PNG.
save_tfs_feature <- function(out_dir, features = tf_features, group_size = 4,
                             reduction = "umap", panel_w = 6, panel_h = 4) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  idx <- split(seq_along(features), ceiling(seq_along(features) / group_size))
  for (chunk in idx) {
    sub <- features[chunk]
    ps <- lapply(seq_along(sub), function(i) {
      FeaturePlot(so, features = sub[i], reduction = reduction, order = TRUE) +
        ggtitle(sprintf("%s (%s)", names(sub)[i], sub[i]))
    })
    ncol <- if (group_size == 1) 1 else min(2, length(sub))
    nrow <- ceiling(length(sub) / ncol)
    p <- wrap_plots(ps, ncol = ncol)
    fname <- if (group_size == 1) {
      sprintf("TF_%s_%s.png", names(sub)[1], sub[1])
    } else {
      sprintf("TFs_%s.png", paste(names(sub), collapse = "-"))
    }
    ggsave(file.path(out_dir, fname), p,
           width = panel_w * ncol, height = panel_h * nrow, dpi = 150)
  }
  cat("Saved", length(idx), "PNG(s) to", out_dir, "\n")
}

plot_tfs_dot <- function(features = tf_features, group.by = NULL) {
  # Collapse duplicate LOCs (e.g. slp1 + slp2 -> same LOC), joining fly names
  by_loc <- split(names(features), unname(features))
  feats  <- names(by_loc)
  lbl    <- setNames(vapply(by_loc, paste, character(1), collapse = "/"), feats)
  p <- DotPlot(so, features = feats, group.by = group.by) + RotatedAxis()
  p + scale_x_discrete(labels = function(x) ifelse(x %in% names(lbl), paste0(lbl[x], "\n", x), x))
}

lookup_fly <- function(fly_gene) {
  hit <- homologs %>% filter(fly_gene == !!fly_gene) %>% arrange(desc(bitscore))
  if (nrow(hit) == 0) { cat("No homolog for", fly_gene, "\n"); return(invisible(NULL)) }
  hit %>% select(fly_gene, bee_gene = sp_transcript, fly_product, pident, bitscore)
}

lookup_loc <- function(loc) {
  homologs %>% filter(sp_transcript == loc) %>%
    arrange(desc(bitscore)) %>%
    select(bee_gene = sp_transcript, fly_gene, fly_product, pident, bitscore)
}

# ---- Inspect ----
DimPlot(so, label = TRUE, repel = TRUE) + NoLegend()
plot_tfs_feature()
plot_tfs_dot()
FeaturePlot(so, features = "LOC410893")  # TfAP-2 per original Rmd

# Dump per-TF featureplots into outputs/figure4/honeybee_TF_featureplots/
# save_tfs_feature("../../outputs/figure4/honeybee_TF_featureplots", group_size = 4)
# save_tfs_feature("../../outputs/figure4/honeybee_TF_featureplots_individual", group_size = 1)

# ---- Starter rename map (from honeybee_seurat_reconstruction.Rmd, leiden1) ----
Idents(so) <- "seurat_clusters"
cluster_names <- c(
  "0"  = "Mi1_I",
  "1"  = "AST/ENS",
  "2"  = "Dm9/R7.8",
  "3"  = "Tm2",
  "4"  = "Dm8/Dm11",
  "5"  = "am_5",
  "6"  = "fl_26/Dm2/Poxn",
  "7"  = "T4/T5_I",
  "8"  = "Mi1_II",
  "9"  = "Tm2",
  "10" = "Dm8",
  "11" = "OPN",
  "12" = "am_12",
  "13" = "T2",
  "14" = "Dm3/N35",
  "15" = "T4/T5_II",
  "16" = "Dm2",
  "17" = "fl_N88",
  "18" = "Tm4/T4/T5",
  "19" = "T4-5/TmY3",
  "20" = "Dm2",
  "21" = "L1/T1/Lai",
  "22" = "L2",
  "23" = "am_23",
  "24" = "L4/L3",
  "25" = "Tm4/T4/T5_II",
  "26" = "am_26",
  "27" = "am_27",
  "28" = "Mi1_III",
  "29" = "am_29",
  "30" = "L3",
  "31" = "Tm2/36",
  "32" = "am_32",
  "33" = "am_33",
  "34" = "am_34",
  "35" = "am_35",
  "36" = "am_36",
  "37" = "Mi9",
  "38" = "am_38",
  "39" = "am_39",
  "40" = "Tm20/N76",
  "41" = "am_41",
  "42" = "am_42",
  "43" = "am_43",
  "44" = "am_44",
  "45" = "am_45"
)

so2 <- RenameIdents(so, cluster_names)
p <- DimPlot(so2, reduction = "umap", label = TRUE, label.size = 2.5,
              label.box = TRUE, repel = TRUE) +
   NoLegend() + NoAxes() + ggtitle("Honeybee Optic Lobe Atlas")
p
ggsave("../../outputs/figure4/UMAP_honeybee_OL_leiden1.pdf", p, width = 10, height = 8)
ggsave("../../outputs/figure4/UMAP_honeybee_OL_leiden1.png", p, width = 10, height = 8)
# saveRDS(so2, out_rds)

save_tfs_feature("../../outputs/figure4/honeybee_TF_featureplots", group_size = 4)

cat("\nReady. Object is `so`. Helpers: plot_tfs_feature(), plot_tfs_dot(), lookup_fly(), lookup_loc().\n")
