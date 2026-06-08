#!/usr/bin/env Rscript
# Neurotransmitter Usage Figure - S. americana Atlas v3
# Panels: A (feature plots 22x22mm), C (UMAP 67x40mm no legend + legend EPS), D (bar 44x40mm)
# Output: EPS + PNG per panel
# Author: David Bellini

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
})

options(Seurat.object.assay.version = "v3")

# ── Paths ──────────────────────────────────────────────────────────────────────
seurat_path   <- "/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.rds"
homologs_path <- "/home/davidbellini/OneDrive/gabbiani/PROJECTS/all_atlases/SCENIC/TF_gene_find/americana_to_fly_homologs.tsv"
out_dir       <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure3/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mm_to_in <- 1 / 25.4

save_eps <- function(plot, filename, width_mm, height_mm) {
  postscript(file.path(out_dir, paste0(filename, ".eps")),
             width = width_mm * mm_to_in, height = height_mm * mm_to_in,
             family = "Helvetica", paper = "special",
             onefile = FALSE, horizontal = FALSE)
  print(plot)
  dev.off()
}

save_png <- function(plot, filename, width_mm, height_mm) {
  ggsave(file.path(out_dir, paste0(filename, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 300, bg = "white")
}

save_legend_eps <- function(plot_with_legend, filename, width_mm = 30, height_mm = 40) {
  leg_plot <- plot_grid(get_legend(plot_with_legend))
  postscript(file.path(out_dir, paste0(filename, ".eps")),
             width = width_mm * mm_to_in, height = height_mm * mm_to_in,
             family = "Helvetica", paper = "special",
             onefile = FALSE, horizontal = FALSE)
  print(leg_plot)
  dev.off()
}

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Loading Seurat object...\n")
so <- readRDS(seurat_path)
DefaultAssay(so) <- "RNA"

homologs <- read_tsv(homologs_path, show_col_types = FALSE)

# ── NT gene mapping ────────────────────────────────────────────────────────────
nt_fly_genes <- c("VGlut", "Gad1", "VGAT")

nt_homologs <- homologs %>%
  filter(fly_gene %in% nt_fly_genes) %>%
  group_by(fly_gene) %>%
  slice_max(order_by = pident, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(seurat_id = gsub("_", "-", sp_transcript))

nt_homologs <- bind_rows(
  nt_homologs,
  tibble(fly_gene = "ChAT",  sp_transcript = "ChAT",  seurat_id = "ChAT"),
  tibble(fly_gene = "VAChT", sp_transcript = "VAChT", seurat_id = "VAChT")
)

rna_genes   <- rownames(so[["RNA"]])
nt_homologs <- nt_homologs %>% mutate(in_seurat = seurat_id %in% rna_genes)

cat("Neurotransmitter gene mapping:\n")
print(nt_homologs %>% select(fly_gene, seurat_id, in_seurat))

nt_present <- nt_homologs %>% filter(in_seurat)
if (nrow(nt_present) == 0) stop("None of the NT genes found in Seurat object.")

# ── NT type assignment ─────────────────────────────────────────────────────────
nt_groups <- list(
  Cholinergic   = nt_present %>% filter(fly_gene %in% c("VAChT", "ChAT")) %>% pull(seurat_id),
  GABAergic     = nt_present %>% filter(fly_gene %in% c("Gad1", "VGAT"))  %>% pull(seurat_id),
  Glutamatergic = nt_present %>% filter(fly_gene == "VGlut")              %>% pull(seurat_id)
)
print(nt_groups)

expr_list <- lapply(names(nt_groups), function(nt_type) {
  genes <- nt_groups[[nt_type]]
  if (length(genes) == 0) return(NULL)
  mat <- GetAssayData(so, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  apply(mat, 2, max) %>% enframe(name = "cell", value = nt_type)
})
expr_list <- Filter(Negate(is.null), expr_list)

nt_expr <- reduce(expr_list, full_join, by = "cell") %>% replace(is.na(.), 0)
nt_cols  <- intersect(names(nt_groups), colnames(nt_expr))

nt_assignment <- nt_expr %>%
  select(cell, all_of(nt_cols)) %>%
  rowwise() %>%
  mutate(
    vals    = list(sort(c_across(all_of(nt_cols)), decreasing = TRUE)),
    top     = vals[[1]][1],
    second  = vals[[1]][2],
    nt_type = case_when(
      top == 0                    ~ "undefined",
      top / (top + second) < 0.55 ~ "undefined",
      TRUE                        ~ nt_cols[which.max(c_across(all_of(nt_cols)))]
    )
  ) %>%
  ungroup() %>%
  select(cell, nt_type)

so$nt_type <- nt_assignment$nt_type[match(colnames(so), nt_assignment$cell)]
so$nt_type[is.na(so$nt_type)] <- "undefined"
print(table(so$nt_type))

# ── Display label mapping ──────────────────────────────────────────────────────
label_map <- c(
  "Cholinergic"   = "Acetylcholine",
  "GABAergic"     = "GABA",
  "Glutamatergic" = "Glutamate",
  "undefined"     = "Undefined"
)
so$nt_label <- label_map[so$nt_type]

label_levels <- c("Acetylcholine", "GABA", "Glutamate", "Undefined")
so$nt_label  <- factor(so$nt_label, levels = intersect(label_levels, unique(so$nt_label)))

# ── Color palette ──────────────────────────────────────────────────────────────
nt_colors_raw <- c(
  "Cholinergic"   = "#E07070",
  "GABAergic"     = "#6BBF6B",
  "Glutamatergic" = "#5BB8D4",
  "undefined"     = "#B89BE0"
)
nt_colors <- setNames(nt_colors_raw, label_map[names(nt_colors_raw)])
nt_colors <- nt_colors[names(nt_colors) %in% levels(so$nt_label)]

# ── UMAP base coords ───────────────────────────────────────────────────────────
umap_base <- as.data.frame(Embeddings(so, reduction = "umap"))
colnames(umap_base) <- c("UMAP_1", "UMAP_2")

# ── Panel A: Feature plots (22x22mm, manual ggplot2) ──────────────────────────
# FeaturePlot pt.size is unreliable for EPS — using geom_point directly
FEATURE_PT_SIZE <- 0.1

feature_theme <- theme_void(base_family = "Helvetica", base_size = 6) +
  theme(
    plot.title    = element_text(face = "bold.italic", size = 6, hjust = 0.5),
    legend.text   = element_text(size = 5, family = "Helvetica"),
    legend.title  = element_text(size = 5, family = "Helvetica"),
    plot.margin   = margin(1, 1, 1, 1, "mm")
  )

thin_colorbar <- guide_colorbar(
  barwidth       = unit(0.15, "cm"),
  barheight      = unit(0.8,  "cm"),
  ticks          = FALSE,
  title.position = "top"
)

panel_a_genes <- c("VAChT", "ChAT", "VGlut", "Gad1", "VGAT")

for (fly_name in panel_a_genes) {
  gene_row <- nt_homologs %>% filter(fly_gene == fly_name)

  if (nrow(gene_row) == 0 || !gene_row$in_seurat[1]) {
    p <- ggplot(umap_base, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(color = "lightgrey", size = FEATURE_PT_SIZE, stroke = 0) +
      ggtitle(fly_name) +
      feature_theme
  } else {
    gene_id  <- gene_row$seurat_id[1]
    expr_vec <- GetAssayData(so, assay = "RNA", layer = "data")[gene_id, ]
    plot_df  <- umap_base
    plot_df$expr <- as.numeric(expr_vec)
    plot_df  <- plot_df[order(plot_df$expr), ]

    p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
      geom_point(size = FEATURE_PT_SIZE, stroke = 0) +
      scale_color_gradientn(
        colors = c("lightgrey", "#FEE08B", "#F46D43", "#A50026"),
        name   = "expr",
        guide  = thin_colorbar
      ) +
      ggtitle(fly_name) +
      feature_theme
  }

  save_eps(p, paste0("figA_", fly_name), 40, 30)
  save_png(p, paste0("figA_", fly_name), 40, 30)
  cat("Saved figA_", fly_name, "\n")
}

# ── Panel C: NT UMAP (67x40mm, no legend) + legend as separate EPS ────────────
umap_df <- umap_base
umap_df$nt_label <- so$nt_label
umap_df <- umap_df %>%
  arrange(factor(nt_label, levels = c("Undefined", "Glutamate", "GABA", "Acetylcholine")))

p_c_with_legend <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = nt_label)) +
  geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(values = nt_colors) +
  labs(color = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.text     = element_text(size = 6, family = "Helvetica"),
    legend.key.size = unit(0.3, "cm"),
    plot.margin     = margin(1, 1, 1, 1, "mm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_c_no_legend <- p_c_with_legend + theme(legend.position = "none")

save_eps(p_c_no_legend, "figC_nt_umap", 80, 60)
save_png(p_c_no_legend, "figC_nt_umap", 80, 60)
save_legend_eps(p_c_with_legend, "figC_nt_umap_legend", width_mm = 20, height_mm = 60)
cat("Saved figC_nt_umap + legend\n")

# ── Panel D: Bar chart (44x40mm) ──────────────────────────────────────────────
nt_pct <- so@meta.data %>%
  count(nt_label) %>%
  mutate(
    pct      = n / sum(n) * 100,
    nt_label = factor(nt_label, levels = label_levels)
  )

p_d <- ggplot(nt_pct, aes(x = nt_label, y = pct, fill = nt_label)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = nt_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Neurotransmitter",
    y = "Dominant neurotransmitter per cell (%)"
  ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 5, family = "Helvetica")
  )

save_eps(p_d, "nt_barplot", 35, 60)
save_png(p_d, "nt_barplot", 35, 60)
cat("Saved nt_barplot\n")

cat("\nNT figures done. Outputs in:", out_dir, "\n")
