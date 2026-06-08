#!/usr/bin/env Rscript
# GO Enrichment Analysis for Successfully Merged Genes

library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(GO.db)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# Read eggNOG-mapper output
eggnog <- read_tsv("/home/davidbellini/OneDrive/gabbiani/PROJECTS/GO/GO_americana_sc/eggnog_files/out.emapper.annotations",
                   skip = 4,
                   comment = "##",
                   show_col_types = FALSE)

# Clean eggNOG query IDs - extract gene ID from full identifier
# Format: gnl|WGS:ZZZZ|egapxtmp_001307-P1 -> egapxtmp_001307
eggnog <- eggnog %>%
  mutate(gene_id = str_extract(`#query`, "egapxtmp_\\d+"))

# Create gene ID to description mapping for plot labels
gene_descriptions <- eggnog %>%
  filter(!is.na(gene_id)) %>%
  dplyr::select(gene_id, Preferred_name, Description) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(
    label = case_when(
      !is.na(Preferred_name) & Preferred_name != "-" ~ Preferred_name,
      !is.na(Description) & Description != "-" ~ str_trunc(Description, 30),
      TRUE ~ gene_id
    )
  ) %>%
  dplyr::select(gene_id, label) %>%
  # Make labels unique by appending gene_id to duplicates
  group_by(label) %>%
  mutate(
    n = n(),
    label = if_else(n > 1, paste0(label, " (", gene_id, ")"), label)
  ) %>%
  ungroup() %>%
  dplyr::select(gene_id, label)

# Extract GO annotations
go_annotations <- eggnog %>%
  dplyr::select(gene_id, GOs) %>%
  filter(!is.na(GOs), GOs != "", !is.na(gene_id)) %>%
  separate_rows(GOs, sep = ",") %>%
  dplyr::rename(GO = GOs) %>%
  # Remove duplicates (same gene can have same GO from different isoforms)
  distinct()


# Create TERM2GENE mapping for clusterProfiler
term2gene_go <- go_annotations %>%
  dplyr::select(GO, gene_id)

# Create TERM2NAME mapping (GO ID -> description)
go_terms <- AnnotationDbi::select(GO.db,
                                  keys = unique(term2gene_go$GO),
                                  columns = c("GOID", "TERM", "ONTOLOGY"),
                                  keytype = "GOID")

# Filter for Biological Process only
go_terms_bp <- go_terms %>%
  filter(ONTOLOGY == "BP") %>%
  dplyr::select(-ONTOLOGY)

term2name_go <- go_terms_bp %>%
  dplyr::rename(GO = GOID, Description = TERM)

# Update TERM2GENE to only include BP terms
term2gene_go <- term2gene_go %>%
  filter(GO %in% term2name_go$GO)


merged_genes_df <- read_csv("/home/davidbellini/OneDrive/gabbiani/TOOLS/Gene-Misannotation-Tool/annotation_comparison/results/merged_genes_only.csv",
                            show_col_types = FALSE)

# Extract unique gene IDs from the new_genes column
merged_gene_ids <- merged_genes_df %>%
  pull(new_genes) %>%
  str_extract_all("egapxtmp_\\d+") %>%
  unlist() %>%
  unique()

genes_with_go <- sum(merged_gene_ids %in% term2gene_go$gene_id)

# Run enrichment using clusterProfiler
enrich_go <- enricher(
  gene = merged_gene_ids,
  TERM2GENE = term2gene_go,
  TERM2NAME = term2name_go,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  minGSSize = 10,      # Minimum 10 genes per GO term
  maxGSSize = 500,     # Maximum 500 genes per GO term
  universe = unique(term2gene_go$gene_id)  # All annotated genes as background
)



results <- enrich_go@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)

OUT_DIR <- "/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure1/"
RESULTS_DIR <- "/home/davidbellini/OneDrive/gabbiani/TOOLS/Gene-Misannotation-Tool/annotation_comparison/results/"

write_csv(results, file.path(RESULTS_DIR, "go_enrichment_results.csv"))

# --- EXCLUDED TERMS (p.adjust < 0.05 but removed from figure) ---------------
# David note: For clarity and lack of space in the figure, some of the statistically significant GO TERMS were removed
# Here is a list of the ones I tossed out. 
# WARNING: THE FOLLOWING DESCRIPTIONS ARE AI-GENERATED
# REDUNDANT GO HIERARCHY (same 6 genes across multiple ontology levels):
#   Mo-molybdopterin cofactor biosynthetic process
#   Mo-molybdopterin cofactor metabolic process
#   molybdopterin cofactor biosynthetic process
#   molybdopterin cofactor metabolic process
#   prosthetic group metabolic process
#
# CROSS-SPECIES ANNOTATION ARTIFACT (not meaningful in locust brain context):
#   regulation of heart rate by cardiac conduction  (3 genes)
#
# BROAD PARENT TERM (redundant with specific GPCR terms already included):
#   G protein-coupled receptor signaling pathway
#
# MARGINAL / UNCLEAR RELEVANCE:
#   positive regulation of Rho protein signal transduction  (4 genes)
#   regulation of Rho protein signal transduction           (6 genes)
#   glycosyl compound metabolic process                     (12 genes)
#
# NOTE: "response to antineoplastic agent" retained in Lipid metabolism —
#   gene set (n=17) is almost entirely leukotriene/arachidonic acid pathway
#   genes; term reflects mammalian GO annotation context.
# NOTE: "neuron fate commitment" placed in Electrophysiology per user decision.
# ----------------------------------------------------------------------------
EXCLUDE_TERMS <- c(
  "Mo-molybdopterin cofactor biosynthetic process",
  "Mo-molybdopterin cofactor metabolic process",
  "molybdopterin cofactor biosynthetic process",
  "molybdopterin cofactor metabolic process",
  "prosthetic group metabolic process",
  "regulation of heart rate by cardiac conduction",
  "G protein-coupled receptor signaling pathway",
  "positive regulation of Rho protein signal transduction",
  "regulation of Rho protein signal transduction",
  "glycosyl compound metabolic process"
)

cat_order <- c(
  "Cuticle development",
  "Lipid metabolism",
  "GPCR / kinase signaling",
  "Synaptic signaling",
  "Planar polarity",
  "Electrophysiology"
)

results_plot <- results %>%
  filter(!Description %in% EXCLUDE_TERMS) %>%
  mutate(category = case_when(
    Description %in% c(
      "chitin-based cuticle development", "cuticle development"
    ) ~ "Cuticle development",
    Description %in% c(
      "leukotriene metabolic process", "prostanoid metabolic process",
      "prostaglandin metabolic process", "icosanoid metabolic process",
      "unsaturated fatty acid metabolic process",
      "fatty acid derivative metabolic process",
      "fatty acid metabolic process", "response to antineoplastic agent"
    ) ~ "Lipid metabolism",
    Description %in% c(
      "adenylate cyclase-inhibiting G protein-coupled receptor signaling pathway",
      "adenylate cyclase-modulating G protein-coupled receptor signaling pathway",
      "G protein-coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger",
      "regulation of protein kinase A signaling"
    ) ~ "GPCR / kinase signaling",
    Description %in% c(
      "modulation of chemical synaptic transmission",
      "regulation of trans-synaptic signaling",
      "regulation of synaptic plasticity",
      "neuromuscular process"
    ) ~ "Synaptic signaling",
    Description %in% c(
      "establishment of ommatidial planar polarity",
      "establishment of planar polarity",
      "establishment of tissue polarity",
      "morphogenesis of a polarized epithelium"
    ) ~ "Planar polarity",
    Description %in% c(
      "membrane depolarization during action potential",
      "cardiac muscle cell action potential",
      "membrane depolarization", "action potential",
      "neuron fate commitment"
    ) ~ "Electrophysiology",
    TRUE ~ "Other"
  )) %>%
  mutate(category = factor(category, levels = cat_order)) %>%
  arrange(desc(category), Count) %>%
  mutate(
    label      = str_trunc(Description, 48, ellipsis = "..."),
    label = factor(label, levels = unique(label))
  )

panel_width_mm  <- 90
panel_height_mm <- 60
mm_to_in        <- 1 / 25.4

p <- ggplot(results_plot, aes(x = Count, y = label, fill = p.adjust)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = label),
            x = 0.3, hjust = 0, size = 4 / .pt, color = "black") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low  = "#a8d5a2",  # Soft green for most significant
                      high = "#f0f9ee",  # Near-white for least significant
                      name = "p-value\n(BH correction)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Gene count", y = NULL) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = 5),
    axis.title.x       = element_text(size = 6),
    strip.text.y.left  = element_text(size = 5, angle = 0, hjust = 1),
    strip.placement    = "outside",
    legend.text        = element_text(size = 5),
    legend.title       = element_text(size = 5),
    legend.key.size    = unit(3, "mm"),
    legend.spacing.y   = unit(0.5, "mm"),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.5, "mm"),
    plot.margin        = margin(1, 1, 1, 1, "mm")
  )

eps_path <- file.path(OUT_DIR, "FIG1_go_enrichment.eps")
postscript(eps_path,
           width = panel_width_mm * mm_to_in,
           height = panel_height_mm * mm_to_in,
           horizontal = FALSE, paper = "special",
           family = "Helvetica")
print(p)
dev.off()
cat("Saved:", eps_path, "\n")

png_path <- file.path(OUT_DIR, "FIG1_go_enrichment.png")
ggsave(png_path, plot = p,
       width = panel_width_mm, height = panel_height_mm,
       units = "mm", dpi = 300)
cat("Saved:", png_path, "\n")
