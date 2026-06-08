# Author: David Bellini
# Bumblebee AC_v1 atlas — cluster inspection & relabeling helper
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(patchwork)

atlas_dir  <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/bumblebee_atlas"
homolog_fp <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/all_atlases/SCENIC/TF_gene_find/bumblebee_to_fly_homologs.tsv"
gff_fp     <- file.path(atlas_dir, "bumblebee_gene_annotations.tsv")
out_rds    <- file.path(atlas_dir, "bumblebee_AC_v1_relabeled.rds")

# ---- Load object (cache .rds next to .h5seurat to avoid re-converting) ----
rds_cache <- file.path(atlas_dir, "bumblebee_AC_v1.rds")
if (file.exists(rds_cache)) {
  so <- readRDS(rds_cache)
} else {
  so <- LoadH5Seurat(file.path(atlas_dir, "bumblebee_AC_v1.h5seurat"))
  saveRDS(so, rds_cache)
}

cat("Assays:", paste(Assays(so), collapse = ", "), "\n")
cat("Default assay:", DefaultAssay(so), "\n")
cat("Cells:", ncol(so), " Genes:", nrow(so), "\n")
cat("Reductions:", paste(Reductions(so), collapse = ", "), "\n")
cat("Meta cols:", paste(colnames(so@meta.data), collapse = ", "), "\n")
print(table(Idents(so)))

gffkey   <- read_tsv(gff_fp,     show_col_types = FALSE)
homologs <- read_tsv(homolog_fp, show_col_types = FALSE)

# ---- Optic-lobe TF panel (fly names) -> bumblebee LOC IDs ----
ol_tfs_fly <- c(
  "dac","hth","tj","vvl","bsh","D","toy","oc","dfr","Optix",
  "ato","so","eya","sv","drf","fkh","ey","Lim1","Lim3","acj6",
  "tll","sim","scro","exex","run","slp1","slp2","ap","disco",
  "br","pdm3","Hr51","Eip93F","fd59A","onecut","Ets65A","CG32532"
)

map_fly_to_bee <- function(fly_genes) {
  homologs %>%
    filter(fly_gene %in% fly_genes) %>%
    group_by(fly_gene) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(fly_gene, bee_gene = sp_gene_name, product = sp_product, pident, bitscore)
}

tf_map <- map_fly_to_bee(ol_tfs_fly)
tf_map$in_object <- tf_map$bee_gene %in% rownames(so)
cat("\nTF mapping (fly -> bumblebee LOC):\n"); print(tf_map, n = Inf)

tf_features <- tf_map %>% filter(in_object) %>% pull(bee_gene)
names(tf_features) <- tf_map %>% filter(in_object) %>% pull(fly_gene)

# ---- Convenience plotters ----
# Use plot.title to show fly name above the LOC ID
plot_tfs_feature <- function(features = tf_features, reduction = "umap", ncol = 4) {
  ps <- lapply(seq_along(features), function(i) {
    FeaturePlot(so, features = features[i], reduction = reduction, order = TRUE) +
      ggtitle(sprintf("%s (%s)", names(features)[i], features[i]))
  })
  wrap_plots(ps, ncol = ncol)
}

# Save FeaturePlots in batches (default 4 per PNG). group_size = 1 -> one TF per PNG.
save_tfs_feature <- function(out_dir, features = tf_features, group_size = 4,
                             reduction = "umap", panel_w = 4, panel_h = 4) {
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

# Lookup helper: fly gene name -> bumblebee LOC (with product)
lookup_fly <- function(fly_gene) {
  hit <- homologs %>% filter(fly_gene == !!fly_gene) %>% arrange(desc(bitscore))
  if (nrow(hit) == 0) { cat("No homolog for", fly_gene, "\n"); return(invisible(NULL)) }
  hit %>% select(fly_gene, bee_gene = sp_gene_name, product = sp_product, pident, bitscore)
}

# Reverse: LOC -> fly gene (for interpreting markers tables)
lookup_loc <- function(loc) {
  homologs %>% filter(sp_gene_name == loc) %>%
    arrange(desc(bitscore)) %>%
    select(bee_gene = sp_gene_name, fly_gene, fly_product, pident, bitscore)
}

# ---- Relabeling workflow (edit & rerun the block below as needed) ----
# Inspect:
  DimPlot(so, label = TRUE, repel = TRUE)
  plot_tfs_feature()
  plot_tfs_dot()
  lookup_fly("tfap2"); lookup_loc("LOC100646599")
  plot_tfs_feature(tf_features[c("dac","hth","tj","Optix")], ncol = 2)
  FeaturePlot(so, features="LOC100645993")
  
#
# Then build a rename map (old_label = "new_label") and apply.
# Starter map adapted from sa3_relabel.R — adjust cluster IDs/labels to bumblebee AC_v1.

# Idents(so) <- "seurat_clusters"

cluster_names <- c(
  "0"  = "PR+",
  "1"  = "A/B-KC",
  "2"  = "A/B-KC",
  "3"  = "ENS/AST",
  "4"  = "N35/N33",
  "5"  = "G-KC",
  "6"  = "N55/Tm3",
  "7"  = "A/B-KC/Tm9a",
  "8"  = "A/B-KC/N102",
  "9"  = "T4",  # manual override: clear TfAP-2 expression; SAMap calls T4-5
  "10" = "T5",  # manual override: no TfAP-2 expression; SAMap calls T4-5
  "11" = "N55",
  "12" = "N69/TmY3",
  "13" = "N1",
  "14" = "N108",
  "15" = "A/B-KC",
  "16" = "fl_123/C3",
  "17" = "G-KC/N1",
  "18" = "TmY5a/Tm29",
  "19" = "TmY5a/N116",
  "20" = "Mi9",
  "21" = "T2/T2a",
  "22" = "bt_22",
  "23" = "Dm3",
  "24" = "Mi15",
  "25" = "Dm3",
  "26" = "Tm2",
  "27" = "Mi1/Tm2",
  "28" = "PNG",
  "29" = "Tm1/Tm2/Tm4",
  "30" = "L1/L2/L3",
  "31" = "Dm8/Dm11",
  "32" = "Mi4",
  "33" = "N34",
  "34" = "C2",
  "35" = "N42",
  "36" = "T1",
  "37" = "L5",
  "38" = "bt_38",
  "39" = "PLM",
  "40" = "Tm9/TmY3",
  "41" = "bt_41",
  "42" = "Dm12"
)

so2 <- RenameIdents(so, cluster_names)
DimPlot(so2, label=TRUE) + NoLegend() 
so[["cell_type"]]       <- Idents(so)
so[["seurat_clusters"]] <- Idents(so)
#
# dir.create("plots", showWarnings = FALSE)
p <- DimPlot(so2, reduction = "umap", label = TRUE, label.size = 2.5,
              label.box = TRUE, repel = TRUE) +
   NoLegend() + NoAxes() + ggtitle("Bumblebee Young Adult Whole Brain Atlas")
p
ggsave("../../outputs/figure4/UMAP_bumblebee_AC_v1.pdf", p, width = 10, height = 8)
ggsave("../../outputs/figure4/UMAP_bumblebee_AC_v1.png", p, width = 10, height = 8)
# saveRDS(so, out_rds)
