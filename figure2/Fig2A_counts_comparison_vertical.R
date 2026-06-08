#!/usr/bin/env Rscript
# Figure 2A (vertical) — featureCounts assignment: Bellini2026 vs NCBI R100
# Nature Ecology and Evolution formatting — vertical bar orientation

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

bellini_summary <- "/home/davidbellini/files/workspace/americana_atlas3/results/feature_counts/feature_counts.txt.summary"
ncbi_summary    <- "/home/davidbellini/files/workspace/americana_atlas3/results_ncbi_original/feature_counts/feature_counts.txt.summary"
outdir          <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 70
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

read_summary <- function(path, label) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   colClasses = c("character", "character"))
  colnames(df) <- c("Status", "count")
  df$count <- as.numeric(gsub(",", "", df$count))
  df$annotation <- label
  df
}

bel <- read_summary(bellini_summary, "Bellini2026")
ncb <- read_summary(ncbi_summary,    "NCBI R100")

summ <- bind_rows(bel, ncb) %>%
  filter(Status %in% c("Assigned", "Unassigned_NoFeatures")) %>%
  mutate(Status = recode(Status,
                         "Assigned"             = "Assigned",
                         "Unassigned_NoFeatures" = "No feature"))

totals <- summ %>% group_by(annotation) %>% summarise(total = sum(count), .groups = "drop")
summ <- summ %>%
  left_join(totals, by = "annotation") %>%
  mutate(pct   = count / total,
         label = sprintf("%.1fM\n(%.1f%%)", count / 1e6, 100 * pct))

assigned <- summ %>%
  filter(Status == "Assigned") %>%
  select(annotation, count) %>%
  pivot_wider(names_from = annotation, values_from = count)
delta_abs <- assigned$Bellini2026 - assigned$`NCBI R100`
delta_pct <- 100 * delta_abs / assigned$`NCBI R100`

summ$annotation <- factor(summ$annotation, levels = c("NCBI R100", "Bellini2026"))
summ$Status     <- factor(summ$Status,     levels = c("No feature", "Assigned"))

pal <- c("Assigned" = "#2E8B8B", "No feature" = "#BEBEBE")

p <- ggplot(summ, aes(x = annotation, y = count / 1e6, fill = Status)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5 / .pt, color = "white", fontface = "bold", lineheight = 0.9,
            family = "Helvetica") +
  scale_fill_manual(values = pal, breaks = c("Assigned", "No feature")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Mapped Reads (Million)", x = NULL, fill = NULL,
       title    = "featureCounts assignment: Bellini2026 vs NCBI Release 100",
       subtitle = sprintf("+%.1fM assigned reads (+%.1f%%) with Bellini2026 annotation",
                          delta_abs / 1e6, delta_pct)) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(size = 5, family = "Helvetica"),
    axis.text.y        = element_text(size = 5, family = "Helvetica"),
    axis.title.y       = element_text(size = 6, family = "Helvetica"),
    legend.position    = "right",
    legend.text        = element_text(size = 5, family = "Helvetica"),
    legend.key.size    = unit(3, "mm"),
    legend.spacing.y   = unit(0.5, "mm"),
    legend.box.spacing = unit(1, "mm"),
    plot.title         = element_text(size = 6, family = "Helvetica", face = "bold"),
    plot.subtitle      = element_text(size = 5, family = "Helvetica", color = "grey30"),
    plot.margin        = margin(1, 1, 1, 1, "mm")
  )

save_eps(p, file.path(outdir, "fig2A_counts_comparison_vertical.eps"))
save_png(p, file.path(outdir, "fig2A_counts_comparison_vertical.png"))

cat(sprintf("\nDone: fig2A_counts_comparison_vertical.eps/.png\n"))
