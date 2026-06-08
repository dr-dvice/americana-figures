#!/usr/bin/env Rscript
# Figure 2E — Optic lobe cell type composition: S. americana vs D. melanogaster
# Nature Ecology and Evolution formatting

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

options(Seurat.object.assay.version = "v3")

output_dir <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 90
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

so2 <- readRDS("/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.rds")

cluster_categories <- c(
  "PNG"  = "Glial", "ENS" = "Glial", "AST" = "Glial",
  "CTX"  = "Glial", "G12" = "Glial", "G57" = "Glial",
  "G74"  = "Glial", "G125" = "Glial",
  "PR"   = "Photoreceptor",
  "C2/C3" = "Centrifugal",
  "L1"   = "Lamina", "L3" = "Lamina", "L2-L4" = "Lamina", "N158/L2" = "Lamina",
  "Mi1"  = "Medulla intrinsic", "Mi15" = "Medulla intrinsic",
  "Transmedullary_neurons1" = "Transmedullary",
  "Transmedullary_neurons2" = "Transmedullary",
  "Tm9"  = "Transmedullary", "Tm3" = "Transmedullary",
  "Tm1/TmY8" = "Transmedullary", "TmY3" = "Transmedullary",
  "TmY5a" = "Transmedullary", "PLM" = "Transmedullary",
  "Dm8/Dm11" = "Distal medulla", "Dm10/Mi4" = "Distal medulla",
  "T1"   = "T neurons", "T2a" = "T neurons", "T2/T3" = "T neurons",
  "T4"   = "T neurons", "T5"  = "T neurons",
  "N164/N115" = "Uncharacterized", "N115/N119" = "Uncharacterized",
  "N96"  = "Uncharacterized", "N34" = "Uncharacterized",
  "N183" = "Uncharacterized", "C86" = "Uncharacterized",
  "C16"  = "Uncharacterized", "sa_40" = "Uncharacterized"
)

unmapped <- setdiff(levels(so2@active.ident), names(cluster_categories))
if (length(unmapped) > 0) warning("Unmapped clusters: ", paste(unmapped, collapse = ", "))

cell_counts_df <- as.data.frame(table(so2@active.ident))
colnames(cell_counts_df) <- c("Cluster", "Count")
cell_counts_df$Percentage <- cell_counts_df$Count / sum(cell_counts_df$Count) * 100
cell_counts_df$Category   <- cluster_categories[as.character(cell_counts_df$Cluster)]

category_counts <- cell_counts_df %>%
  group_by(Category) %>%
  summarise(Count = sum(Count), Percentage = sum(Percentage), .groups = "drop") %>%
  arrange(desc(Percentage))

fly_stats_dir    <- "/home/davidbellini/OneDrive/gabbiani/SAMap_fly_datasets/fly_atlas_statistics"
zipursky_summary <- read.csv(file.path(fly_stats_dir,
                    "zipursky_2020_96h_projection_summary_HARMONIZED.csv"))
zipursky_summary$projection_category <- gsub("_", " ", zipursky_summary$projection_category)
zipursky_summary$Species <- "D. melanogaster"

schistocerca_summary <- category_counts %>%
  rename(projection_category = Category) %>%
  select(projection_category, percentage = Percentage) %>%
  mutate(Species = "S. americana")

combined_data <- bind_rows(
  zipursky_summary %>% select(projection_category, percentage, Species),
  schistocerca_summary
) %>% filter(!is.na(projection_category))

all_categories <- unique(combined_data$projection_category)
combined_data_complete <- combined_data %>%
  complete(Species, projection_category = all_categories, fill = list(percentage = 0))

category_order <- combined_data_complete %>%
  group_by(projection_category) %>%
  summarise(total = sum(percentage), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(projection_category)

combined_data_complete$projection_category <- factor(
  combined_data_complete$projection_category, levels = category_order)

p <- ggplot(combined_data_complete,
            aes(x = projection_category, y = percentage, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c("S. americana"    = "#E69F00",
               "D. melanogaster" = "#56B4E9"),
    labels = c(
      expression(italic("S. americana")),
      expression(italic("D. melanogaster"))
    )
  ) +
  labs(x = "Cell type category", y = "% of total cells", fill = NULL) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 5, family = "Helvetica"),
    axis.text.y        = element_text(size = 5, family = "Helvetica"),
    axis.title.x       = element_text(size = 6, family = "Helvetica"),
    axis.title.y       = element_text(size = 6, family = "Helvetica"),
    legend.position    = "right",
    legend.text        = element_text(size = 5, face = "italic", family = "Helvetica"),
    legend.key.size    = unit(3, "mm"),
    legend.spacing.y   = unit(0.5, "mm"),
    plot.margin        = margin(1, 1, 1, 1, "mm")
  )

save_eps(p, file.path(output_dir, "fig2E_celltype_composition.eps"))
save_png(p, file.path(output_dir, "fig2E_celltype_composition.png"))
cat("Done: fig2E_celltype_composition.eps/.png\n")
