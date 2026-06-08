# Author: David Bellini
# M. pharaonis male OL atlas — cluster inspection & relabeling helper
# Run interactively: /opt/R/4.2.3/bin/R
# NOTE: rownames in this object are *fly gene symbols* (not LOC IDs),
#       so no homolog mapping is needed for TF lookups.
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(patchwork)

atlas_dir  <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/pharaonis"
homolog_fp <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/all_atlases/SCENIC/TF_gene_find/pharaonis_to_fly_homologs.tsv"
out_rds    <- file.path(atlas_dir, "pharaonis_maleonly_ol_relabeled.rds")

# ---- Load object (cache .rds next to .h5seurat to avoid re-converting) ----
rds_cache <- file.path(atlas_dir, "pharaonis_male_subclusters.rds")
if (file.exists(rds_cache)) {
  so <- readRDS(rds_cache)
} else {
  so <- LoadH5Seurat(file.path(atlas_dir, "pharaonis_male_subclusters.h5seurat"))
  saveRDS(so, rds_cache)
}

cat("Assays:", paste(Assays(so), collapse = ", "), "\n")
cat("Default assay:", DefaultAssay(so), "\n")
cat("Cells:", ncol(so), " Genes:", nrow(so), "\n")
cat("Reductions:", paste(Reductions(so), collapse = ", "), "\n")
cat("Meta cols:", paste(colnames(so@meta.data), collapse = ", "), "\n")
print(table(Idents(so)))

homologs <- read_tsv(homolog_fp, show_col_types = FALSE)

# ---- Optic-lobe TF panel (fly names — used directly as gene names) ----
ol_tfs_fly <- c(
  "dac","hth","tj","vvl","bsh","D","toy","oc","dfr","Optix",
  "ato","so","eya","sv","drf","fkh","ey","Lim1","Lim3","acj6",
  "tll","sim","scro","exex","run","slp1","slp2","ap","disco",
  "br","pdm3","Hr51","Eip93F","fd59A","onecut","Ets65A","CG32532",
  "TfAP-2"
)

# Mixed namespace: rownames are mostly fly-style symbols, but some genes are
# stored under a different (pharaonis-specific) symbol. Strategy:
#   1. If fly symbol is in rownames -> use it directly.
#   2. Otherwise, look up best homolog by sp_gene_name and use that.
map_fly_to_ant <- function(fly_genes) {
  rn <- rownames(so)
  best <- homologs %>%
    filter(fly_gene %in% fly_genes, !is.na(sp_gene_name)) %>%
    group_by(fly_gene) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(fly_gene, ant_symbol = sp_gene_name, transcript = sp_transcript, bitscore)

  tibble(fly_gene = fly_genes) %>%
    left_join(best, by = "fly_gene") %>%
    mutate(
      ant_gene = case_when(
        fly_gene %in% rn   ~ fly_gene,
        ant_symbol %in% rn ~ ant_symbol,
        TRUE ~ NA_character_
      ),
      via = case_when(
        fly_gene %in% rn   ~ "direct",
        ant_symbol %in% rn ~ "homolog",
        TRUE ~ "missing"
      )
    )
}

tf_map <- map_fly_to_ant(ol_tfs_fly)
cat("\nTF mapping (fly -> pharaonis symbol):\n"); print(tf_map, n = Inf)

tf_features <- setNames(tf_map$ant_gene[!is.na(tf_map$ant_gene)],
                        tf_map$fly_gene[!is.na(tf_map$ant_gene)])

# ---- Convenience plotters ----
plot_tfs_feature <- function(features = tf_features, reduction = "umap", ncol = 4) {
  ps <- lapply(seq_along(features), function(i) {
    title <- if (names(features)[i] == features[i]) names(features)[i]
             else sprintf("%s (%s)", names(features)[i], features[i])
    FeaturePlot(so, features = features[i], reduction = reduction, order = TRUE) +
      ggtitle(title)
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
      title <- if (names(sub)[i] == sub[i]) names(sub)[i]
               else sprintf("%s (%s)", names(sub)[i], sub[i])
      FeaturePlot(so, features = sub[i], reduction = reduction, order = TRUE) +
        ggtitle(title)
    })
    ncol <- if (group_size == 1) 1 else min(2, length(sub))
    nrow <- ceiling(length(sub) / ncol)
    p <- wrap_plots(ps, ncol = ncol)
    safe <- gsub("[/ ]", "-", names(sub))
    fname <- if (group_size == 1) sprintf("TF_%s.png", safe[1])
             else sprintf("TFs_%s.png", paste(safe, collapse = "-"))
    ggsave(file.path(out_dir, fname), p,
           width = panel_w * ncol, height = panel_h * nrow, dpi = 150)
  }
  cat("Saved", length(idx), "PNG(s) to", out_dir, "\n")
}

plot_tfs_dot <- function(features = tf_features, group.by = NULL) {
  DotPlot(so, features = unname(features), group.by = group.by) + RotatedAxis()
}

# Lookup helpers (still useful for non-TF genes that might be stored as LOC/transcript IDs)
lookup_fly <- function(fly_gene) {
  hit <- homologs %>% filter(fly_gene == !!fly_gene) %>% arrange(desc(bitscore))
  if (nrow(hit) == 0) { cat("No homolog for", fly_gene, "\n"); return(invisible(NULL)) }
  hit %>% select(fly_gene, ant_gene = sp_gene_name, transcript = sp_transcript, pident, bitscore)
}

lookup_ant <- function(name) {
  homologs %>% filter(sp_gene_name == name | sp_transcript == name) %>%
    arrange(desc(bitscore)) %>%
    select(ant_gene = sp_gene_name, transcript = sp_transcript, fly_gene, fly_product, pident, bitscore)
}

# ---- Inspect ----
# DimPlot(so, label = TRUE, repel = TRUE) + NoLegend()
# plot_tfs_feature()
# plot_tfs_dot()
# FeaturePlot(so, features = "TfAP-2")

# ---- Starter rename map (from pharaonisRDS.Rmd) ----
# Idents(so) <- "seurat_clusters"
cluster_names <- c(
  "0"  = "Tm20",
  "1"  = "Tm3",
  "2"  = "102",
  "3"  = "mp_3",
  "4"  = "mp_4",
  "5"  = "fl_30",
  "6"  = "ENS1",
  "7"  = "44",
  "8"  = "T4/T5_I",
  "9"  = "Mi9",
  "10" = "L1/L2/L2/L3",
  "11" = "G-KC1",
  "12" = "G-KC2",
  "13" = "T4/T5_II",
  "14" = "mp_14",
  "15" = "G-KC",
  "16" = "fl_26",
  "17" = "fl_102",
  "18" = "mp_18",
  "19" = "ENS2",
  "20" = "G-KC",
  "21" = "102",
  "22" = "N100/Mi4/14",
  "23" = "mp_23",
  "24" = "G-KC",
  "25" = "mp_25",
  "26" = "mp_26",
  "27" = "G-KC/TmY8",
  "28" = "L2",
  "29" = "Mi1/Tm1/TmY8",
  "30" = "PNG",
  "31" = "Tm1/Tm2",
  "32" = "N134/Poxn",
  "33" = "AST",
  "34" = "mp_34",
  "35" = "Mi9/TmY5a/178",
  "36" = "Mi15",
  "37" = "ENS3",
  "38" = "T4-5",
  "39" = "Tm1",
  "40" = "L1",
  "41" = "Tm9",
  "42" = "mp_42",
  "43" = "mp_43",
  "44" = "AST",
  "45" = "mp_45",
  "46" = "mp_46",
  "47" = "SER/OCTY",
  "48" = "fl_G74",
  "49" = "fl_C46",
  "50" = "mp_50",
  "51" = "CTX",
  "52" = "fl_90/C86",
  "53" = "mp_53",
  "54" = "ENS4",
  "55" = "fl_C68",
  "56" = "ENS5",
  "57" = "mp_57",
  "58" = "mp_58",
  "59" = "PR",
  "60" = "mp_60",
  "61" = "mp_61",
  "62" = "mp_62",
  "63" = "mp_63"
)

so2 <- RenameIdents(so, cluster_names)
p <- DimPlot(so2, reduction = "umap", label = TRUE, label.size = 2.5,
              label.box = TRUE, repel = TRUE) +
   NoLegend() + NoAxes() + ggtitle("M. pharaonis Male Optic Lobe Atlas")
p
ggsave("../../outputs/figure4/UMAP_pharaonis_male_OL.pdf", p, width = 12, height = 9)
ggsave("../../outputs/figure4/UMAP_pharaonis_male_OL.png", p, width = 12, height = 9)
saveRDS(so2, out_rds)
#
save_tfs_feature("../../outputs/figure4/pharaonis_TF_featureplots", group_size = 4)
