#!/usr/bin/env Rscript
# Figure 4D — Cross-species conservation confidence tiles (publication EPS)
# Axes inverted vs prototype: canonical types on x, species on y, family facets as columns
# Output: EPS + PNG at 90x60mm
# Author: David Bellini

library(tidyverse)

harm_path <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/scripts/figure4/harmonization.csv"
out_dir   <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure4/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_width_mm  <- 90
panel_height_mm <- 80
mm_to_in        <- 1 / 25.4

species_order <- c("S.americana", "B.terrestris", "A.mellifera", "M.pharaonis")

harm <- read_csv(harm_path, show_col_types = FALSE,
                 col_types = cols(
                   family = "c", canonical = "c",
                   sa = "i", sa_lab = "c", bt = "i", bt_lab = "c",
                   hb = "i", hb_lab = "c", mp = "i", mp_lab = "c"
                 )) %>%
  mutate(across(ends_with("_lab"), ~replace_na(., "")))

conf_long <- harm %>%
  pivot_longer(c(sa, bt, hb, mp), names_to = "species_code", values_to = "confidence") %>%
  left_join(
    harm %>%
      select(family, canonical, sa_lab, bt_lab, hb_lab, mp_lab) %>%
      pivot_longer(c(sa_lab, bt_lab, hb_lab, mp_lab),
                   names_to = "species_code", values_to = "local_label") %>%
      mutate(species_code = sub("_lab$", "", species_code)),
    by = c("family", "canonical", "species_code")
  ) %>%
  mutate(
    species = recode(species_code,
                     sa = "S.americana", bt = "B.terrestris",
                     hb = "A.mellifera",  mp = "M.pharaonis"),
    species = factor(species, levels = species_order),
    conf_lbl = factor(case_when(
      confidence == 2 ~ "confident",
      confidence == 1 ~ "tentative",
      TRUE            ~ "absent"
    ), levels = c("absent", "tentative", "confident"))
  )

# Order canonical types within each family by n_confident desc, then n_any desc
row_order <- conf_long %>%
  group_by(family, canonical) %>%
  summarise(n_conf = sum(confidence == 2), n_any = sum(confidence > 0), .groups = "drop") %>%
  arrange(family, -n_conf, -n_any)

conf_long <- conf_long %>%
  mutate(canonical = factor(canonical, levels = rev(row_order$canonical)))

pal <- c(absent = "grey92", tentative = "#fdbb84", confident = "#2b8cbe")

p <- ggplot(conf_long, aes(species, canonical, fill = conf_lbl)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(confidence > 0, local_label, "")),
            size = 4 / .pt, color = "black", family = "Helvetica") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_discrete(position = "top") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Helvetica", base_size = 6) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, hjust = 1, face = "bold",
                                      size = 5, family = "Helvetica"),
    panel.grid         = element_blank(),
    axis.text.x.top    = element_text(face = "italic", size = 6,
                                      family = "Helvetica"),
    axis.text.y        = element_text(size = 5, family = "Helvetica"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 5, family = "Helvetica"),
    legend.key.size    = unit(0.25, "cm"),
    plot.margin        = margin(1, 1, 1, 1, "mm")
  ) +
  guides(fill = guide_legend(nrow = 1))

postscript(file.path(out_dir, "figD_conservation_tiles.eps"),
           width = panel_width_mm * mm_to_in, height = panel_height_mm * mm_to_in,
           family = "Helvetica", paper = "special",
           onefile = FALSE, horizontal = FALSE)
print(p)
dev.off()

ggsave(file.path(out_dir, "figD_conservation_tiles.png"), p,
       width = panel_width_mm, height = panel_height_mm,
       units = "mm", dpi = 300, bg = "white")

cat("Saved figD_conservation_tiles.eps/.png to", out_dir, "\n")
