#!/usr/bin/env Rscript
# Neuromodulator Figure - S. americana Atlas v3
# Panels: A (feature plots 22x22mm), C (UMAP 67x40mm no legend + legend EPS), D (bar 44x40mm)
# Output: EPS + PNG per panel
# Author: David Bellini

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(RColorBrewer)
})

options(Seurat.object.assay.version = "v3")

# ── Paths ──────────────────────────────────────────────────────────────────────
seurat_path <- "/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.rds"
out_dir     <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure3/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

EXPR_THRESHOLD <- 0.5
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
rna_genes <- rownames(so[["RNA"]])

# ── Gene definitions ───────────────────────────────────────────────────────────
nm_genes <- tribble(
  ~label,                        ~fly_ref, ~seurat_id,
  "Serotonergic",                "Trh",    "Trh",
  "Dopaminergic",                "pale",   "egapxtmp-011035",
  "Tyraminergic/Octopaminergic", "Tdc2",   "egapxtmp-007839",
  "Histaminergic",               "Hdc",    "egapxtmp-024608",
  "Nitrergic",                   "Nos",    "egapxtmp-001450",
  "Vmat",                        "Vmat",   "egapxtmp-000921"
) %>%
  mutate(in_seurat = seurat_id %in% rna_genes)

cat("Neuromodulator gene mapping:\n")
print(nm_genes)

missing <- nm_genes %>% filter(!in_seurat)
if (nrow(missing) > 0) {
  cat("\nWARNING - not found in RNA assay:\n")
  print(missing)
}

# ── Expression extraction ──────────────────────────────────────────────────────
nm_present <- nm_genes %>% filter(in_seurat)

expr_df <- map(nm_present$seurat_id, function(g) {
  GetAssayData(so, assay = "RNA", layer = "data")[g, , drop = TRUE] %>%
    enframe(name = "cell", value = g)
}) %>%
  reduce(full_join, by = "cell") %>%
  replace(is.na(.), 0)

# ── Type assignment ────────────────────────────────────────────────────────────
# Vmat excluded (shared marker); Nitrergic shown individually but not in UMAP assignment
category_genes <- nm_present %>%
  filter(!label %in% c("Vmat", "Nitrergic")) %>%
  select(label, seurat_id)

cat(sprintf("\nExpression threshold: %.2f (log-normalized)\n", EXPR_THRESHOLD))

assignment <- expr_df %>%
  select(cell, all_of(category_genes$seurat_id)) %>%
  rowwise() %>%
  mutate(
    vals    = list(c_across(all_of(category_genes$seurat_id))),
    max_val = max(unlist(vals)),
    nm_type = if (max_val < EXPR_THRESHOLD) "undefined" else
                category_genes$label[which.max(unlist(vals))]
  ) %>%
  ungroup() %>%
  select(cell, nm_type)

so$nm_type <- assignment$nm_type[match(colnames(so), assignment$cell)]
so$nm_type[is.na(so$nm_type)] <- "undefined"

cat("\nCell counts per neuromodulator type:\n")
print(table(so$nm_type))

# ── Display label mapping ──────────────────────────────────────────────────────
nm_label_map <- c(
  "Serotonergic"               = "Serotonin",
  "Dopaminergic"               = "Dopamine",
  "Tyraminergic/Octopaminergic"= "Tyramine/Octopamine",
  "Histaminergic"              = "Histamine",
  "undefined"                  = "Undefined"
)
so$nm_label <- nm_label_map[so$nm_type]

nm_label_levels <- c("Serotonin", "Dopamine", "Tyramine/Octopamine", "Histamine", "Undefined")
so$nm_label <- factor(so$nm_label, levels = intersect(nm_label_levels, unique(so$nm_label)))

# ── Color palette ──────────────────────────────────────────────────────────────
s2 <- brewer.pal(8, "Set2")
nm_colors_raw <- c(
  "Serotonergic"               = s2[2],
  "Dopaminergic"               = s2[1],
  "Tyraminergic/Octopaminergic"= s2[3],
  "Histaminergic"              = s2[4],
  "undefined"                  = "#CCCCCC"
)
nm_colors <- setNames(nm_colors_raw, nm_label_map[names(nm_colors_raw)])
nm_colors <- nm_colors[names(nm_colors) %in% levels(so$nm_label)]

# ── UMAP base coords ───────────────────────────────────────────────────────────
umap_base <- as.data.frame(Embeddings(so, reduction = "umap"))
colnames(umap_base) <- c("UMAP_1", "UMAP_2")

# ── Panel A: Feature plots (22x22mm, manual ggplot2) ──────────────────────────
# DimPlot/FeaturePlot pt.size unreliable in EPS — using geom_point directly
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

for (i in seq_len(nrow(nm_genes))) {
  row        <- nm_genes[i, ]
  plot_title <- sprintf("%s (%s)", row$fly_ref, row$label)

  if (!row$in_seurat) {
    p <- ggplot(umap_base, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(color = "lightgrey", size = FEATURE_PT_SIZE, stroke = 0) +
      ggtitle(plot_title) +
      labs(caption = "not detected") +
      feature_theme
  } else {
    expr_vec     <- GetAssayData(so, assay = "RNA", layer = "data")[row$seurat_id, ]
    plot_df      <- umap_base
    plot_df$expr <- as.numeric(expr_vec)
    plot_df      <- plot_df[order(plot_df$expr), ]

    p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
      geom_point(size = FEATURE_PT_SIZE, stroke = 0) +
      scale_color_gradientn(
        colors = c("lightgrey", "#FEE08B", "#F46D43", "#A50026"),
        name   = "expr",
        guide  = thin_colorbar
      ) +
      ggtitle(plot_title) +
      feature_theme
  }

  safe_name <- gsub("[^A-Za-z0-9_]", "_", row$label)
  save_eps(p, paste0("figA_", safe_name), 40, 30)
  save_png(p, paste0("figA_", safe_name), 40, 30)
  cat("Saved figA_", safe_name, "\n")
}

# ── Panel C: NM UMAP (67x40mm, no legend) + legend as separate EPS ────────────
umap_df <- umap_base
umap_df$nm_label <- so$nm_label
render_order <- c("Undefined", "Serotonin", "Dopamine", "Tyramine/Octopamine", "Histamine")
render_order <- intersect(render_order, levels(so$nm_label))
umap_df <- umap_df %>% arrange(factor(nm_label, levels = render_order))

p_c_with_legend <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = nm_label)) +
  geom_point(size = 0.15, stroke = 0) +
  scale_color_manual(values = nm_colors) +
  labs(color = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.text     = element_text(size = 6, family = "Helvetica"),
    legend.key.size = unit(0.3, "cm"),
    plot.margin     = margin(1, 1, 1, 1, "mm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_c_no_legend <- p_c_with_legend + theme(legend.position = "none")

save_eps(p_c_no_legend, "figC_nm_umap", 80, 60)
save_png(p_c_no_legend, "figC_nm_umap", 80, 60)
save_legend_eps(p_c_with_legend, "figC_nm_umap_legend", width_mm = 20, height_mm = 60)
save_png(p_c_with_legend, "figC_nm_umap_with_legend_ref", 100, 60)
cat("Saved figC_nm_umap + legend\n")

# ── Panel D: Bar chart (44x40mm) ──────────────────────────────────────────────
nm_pct <- so@meta.data %>%
  count(nm_label) %>%
  mutate(
    pct      = n / sum(n) * 100,
    nm_label = factor(nm_label, levels = nm_label_levels)
  ) %>%
  filter(nm_label != "Undefined") %>%
  mutate(nm_label = droplevels(nm_label))

p_d <- ggplot(nm_pct, aes(x = nm_label, y = pct, fill = nm_label)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = nm_colors) +
  scale_y_continuous(limits = c(0, 20), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Neuromodulator",
    y = "Dominant neuromodulator per cell (%)"
  ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 5, family = "Helvetica")
  )

save_eps(p_d, "nm_barplot", 35, 60)
save_png(p_d, "nm_barplot", 35, 60)
cat("Saved nm_barplot\n")

cat("\nNM figures done. Outputs in:", out_dir, "\n")
