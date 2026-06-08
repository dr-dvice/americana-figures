#!/usr/bin/env Rscript
# Figure 2D — TF expression violin plots in S. americana T neurons
# Nature Ecology and Evolution formatting

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

options(Seurat.object.assay.version = "v3")

rds_path <- "/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.rds"
out_dir  <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

panel_width_mm  <- 70
panel_height_mm <- 70
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

t_neuron_markers <- c(
  "egapxtmp-015520",  # TfAP-2
  "egapxtmp-011504",  # acj6
  "egapxtmp-020439",  # Lim1
  "egapxtmp-015484",  # dac
  "egapxtmp-010971",  # oc
  "egapxtmp-011685",  # ap
  "egapxtmp-006463",  # Drgx
  "egapxtmp-018167",  # Sox102F
  "egapxtmp-020700",  # bi
  "egapxtmp-001029",  # grn
  "egapxtmp-007865",  # pdm3
  "egapxtmp-011077"   # scro
)

t_marker_names <- c(
  "egapxtmp-015520" = "TfAP-2",
  "egapxtmp-011504" = "acj6",
  "egapxtmp-020439" = "Lim1",
  "egapxtmp-015484" = "dac",
  "egapxtmp-010971" = "oc",
  "egapxtmp-011685" = "ap",
  "egapxtmp-006463" = "Drgx",
  "egapxtmp-018167" = "Sox102F",
  "egapxtmp-020700" = "bi",
  "egapxtmp-001029" = "grn",
  "egapxtmp-007865" = "pdm3",
  "egapxtmp-011077" = "scro"
)

t_cluster_order <- c("T1", "T2a", "T2/T3", "T4", "T5")
stopifnot(all(t_cluster_order %in% levels(Idents(so))))

so_t <- subset(so, idents = t_cluster_order)
Idents(so_t) <- factor(Idents(so_t), levels = t_cluster_order)
DefaultAssay(so_t) <- "RNA"

ncol_grid <- 4
n_genes   <- length(t_neuron_markers)
nrow_grid <- ceiling(n_genes / ncol_grid)  # 3 rows

vln_plots <- lapply(seq_along(t_neuron_markers), function(i) {
  gene      <- t_neuron_markers[i]
  row_idx   <- ceiling(i / ncol_grid)
  is_bottom <- row_idx == nrow_grid

  p <- VlnPlot(so_t, features = gene, pt.size = 0) +
    ggtitle(t_marker_names[gene]) +
    NoLegend() +
    labs(x = NULL, y = NULL) +
    theme(
      text              = element_text(family = "Helvetica"),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(size = 6, family = "Helvetica", face = "italic",
                                       hjust = 0.5, margin = margin(b = 0.5, unit = "mm")),
      axis.text.y       = element_text(size = 5, family = "Helvetica"),
      axis.ticks.y      = element_line(linewidth = 0.2),
      axis.text.x       = if (is_bottom)
                            element_text(size = 5, family = "Helvetica", angle = 45, hjust = 1)
                          else
                            element_blank(),
      axis.ticks.x      = if (is_bottom) element_line(linewidth = 0.2) else element_blank(),
      plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "mm")
    )
  # Thin the violin outlines — VlnPlot doesn't expose linewidth directly
  for (j in seq_along(p$layers)) {
    if (inherits(p$layers[[j]]$geom, "GeomViolin")) {
      p$layers[[j]]$aes_params$linewidth <- 0.2
    }
  }
  p
})

# Shared y-axis label as a narrow spacer panel
y_label <- ggplot() +
  annotate("text", x = 0, y = 0,
           label  = "Expression ln(counts + 1)",
           angle  = 90,
           size   = 6 / .pt,
           family = "Helvetica") +
  theme_void() +
  coord_cartesian(clip = "off")

p_grid <- wrap_plots(vln_plots, ncol = ncol_grid)

px <- y_label + p_grid +
  plot_layout(widths = c(0.04, 1)) +
  plot_annotation(
    title = expression("TF expression in " ~ italic("Schistocerca americana")),
    theme = theme(
      plot.title = element_text(size = 7, family = "Helvetica", hjust = 0.5,
                                margin = margin(b = 1, unit = "mm"))
    )
  )

save_eps(px, file.path(out_dir, "fig2D_Tneuron_TFs.eps"))
save_png(px, file.path(out_dir, "fig2D_Tneuron_TFs.png"))
cat("Done: fig2D_Tneuron_TFs.eps/.png\n")
