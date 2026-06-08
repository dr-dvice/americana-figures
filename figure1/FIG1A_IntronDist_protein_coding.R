library(dplyr)
library(ggplot2)
library(scales)
library(purrr)

out_path <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure1/FIG1A_intron_length_distribution_protein_coding.eps"

save_png <- function(plot, eps_path, width_mm, height_mm) {
  png_path <- sub("\\.eps$", ".png", eps_path)
  ggsave(png_path, plot, width = width_mm, height = height_mm, units = "mm", dpi = 300)
  cat("Saved:", png_path, "\n")
}

# === 1. Set file paths ===
data_dir <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/results/intermediate/intron_sizes_protein_coding_only/"
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

# Create species names from filenames
species_names <- gsub("_introns\\.csv$", "", basename(files))

# Map species names to proper scientific names
species_mapping <- c(
  "BomTerr"                     = "Bombus terrestris",
  "fly"                         = "Drosophila melanogaster",
  "mellifera"                   = "Apis mellifera",
  "schistocercaAmericana"       = "Schistocerca americana",
  "schistocercaCancellata"      = "Schistocerca cancellata",
  "schistocercaGregaria"        = "Schistocerca gregaria",
  "schistocercaNitens"          = "Schistocerca nitens",
  "schistocercaPiceifrons"      = "Schistocerca piceifrons",
  "schistocercaSerialisCubense" = "Schistocerca serialis cubense",
  "lekang_locusta"              = "Locusta migratoria",
  "locusta_migratoria"          = "Locusta migratoria",
  "pharaonis"                   = "Monomorium pharaonis",
  "saltator"                    = "Harpegnathos saltator"
)

# Apply mapping (if a name is not in the mapping, keep original with underscores replaced by spaces)
species_names <- sapply(species_names, function(x) {
  if (x %in% names(species_mapping)) {
    species_mapping[x]
  } else {
    gsub("_", " ", x)
  }
})

# Filter out sheina2025 (the in-silico annotation)
exclude_pattern <- "sheina2025"
keep_indices <- !grepl(exclude_pattern, basename(files))
files <- files[keep_indices]
species_names <- species_names[keep_indices]

# ---- 2. Read and combine data ----
all_species_data <- map2_dfr(
  files,
  species_names,
  ~ read.csv(.x) %>%
    select(intron_size) %>%
    filter(!is.na(intron_size), intron_size > 0) %>%
    mutate(
      intron_length = as.numeric(intron_size),
      Species = .y
    )
)

# ---- 3. Define custom color palette ----
unique_species <- unique(all_species_data$Species)

color_palette <- c()
schistocerca_colors <- c("#1a5490", "#2970b8", "#3d8cd0", "#5ca8e8", "#85c1f0", "#aed6f8", "#d4ebfc")
schistocerca_count <- 1

for (sp in unique_species) {
  if (grepl("^Schistocerca", sp)) {
    color_palette[sp] <- schistocerca_colors[schistocerca_count]  # Shades of blue for Schistocerca species
    schistocerca_count <- schistocerca_count + 1
  } else if (grepl("Bombus", sp)) {
    color_palette[sp] <- "#FDD835"  # Bright yellow for bumblebee
  } else if (grepl("Apis", sp)) {
    color_palette[sp] <- "#F4511E"  # Deep orange for honeybee
  } else if (grepl("Drosophila", sp)) {
    color_palette[sp] <- "#D32F2F"  # Red for fly
  } else if (grepl("Locusta", sp)) {
    color_palette[sp] <- "#388E3C"  # Green for locust
  } else if (grepl("Monomorium", sp)) {
    color_palette[sp] <- "#8D6E63"  # Light brown for Monomorium ant
  } else if (grepl("Harpegnathos", sp)) {
    color_palette[sp] <- "#4E342E"  # Very dark brown for Harpegnathos ant
  } else {
    color_palette[sp] <- "#757575"  # Gray for others
  }
}

# ---- 4. Frequency polygon ----
# Pre-compute per-species proportions explicitly to avoid after_stat(count/sum(count))
# summing globally across all species instead of per species.
all_species_data$log10_length <- log10(all_species_data$intron_length)

data_range <- range(all_species_data$log10_length)
binwidth    <- diff(data_range) / 60
breaks      <- seq(data_range[1], data_range[2], by = binwidth)

binned_data <- all_species_data %>%
  mutate(bin = cut(log10_length, breaks = breaks, include.lowest = TRUE),
         bin_mid = breaks[findInterval(log10_length, breaks)]) %>%
  filter(!is.na(bin), log10_length >= 1) %>%
  group_by(Species, bin_mid) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Species) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_freqpoly <- ggplot(binned_data, aes(x = bin_mid, y = prop, color = Species)) +
  geom_line(linewidth = 0.35) +
  scale_x_continuous(
    limits = c(1, 6),
    breaks = 1:6,
    labels = math_format(10^.x),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = color_palette) +
  labs(
    x     = "Intron length (bp)",
    y     = "Frequency (proportion of introns)",
    color = NULL
  ) +
  theme_minimal(base_size = 6, base_family = "Helvetica") +
  theme(
    text             = element_text(family = "Helvetica", size = 6),
    axis.title       = element_text(size = 6),
    axis.text        = element_text(size = 5, color = "black"),
    axis.text.x      = element_text(size = 5, color = "black"),
    legend.text      = element_text(size = 5, face = "italic"),
    legend.key.size  = unit(3, "mm"),
    legend.spacing.y = unit(0.5, "mm"),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm")
  )

# ---- 5. Save ----
# width = 88mm, height = 65mm per panel in a 2x3 layout (180mm total figure width)
# Helvetica is a native PostScript font so no embedding issues
mm_to_in        <- 1 / 25.4
panel_width_mm  <- 90
panel_height_mm <- 60

postscript(
  file       = out_path,
  width      = panel_width_mm  * mm_to_in,
  height     = panel_height_mm * mm_to_in,
  family     = "Helvetica",
  horizontal = FALSE,
  paper      = "special",
  onefile    = FALSE
)
print(p_freqpoly)
dev.off()

cat("Saved:", out_path, "\n")
save_png(p_freqpoly, out_path, panel_width_mm, panel_height_mm)
