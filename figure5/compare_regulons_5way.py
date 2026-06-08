"""
5-Way Cross-Species Regulon TF Comparison (Fly-Hub Approach)

Maps all species' regulon TFs to Drosophila namespace via BLAST homology
files at >=50% protein identity, then computes set intersections.

Species: D. melanogaster, A. mellifera, M. pharaonis, S. americana, B. terrestris

Usage: conda run -n bio python compare_regulons_5way.py
"""

import pandas as pd
import numpy as np
from itertools import combinations
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from upsetplot import UpSet, from_memberships
import warnings
warnings.filterwarnings('ignore')

# --- Nature Ecology & Evolution formatting (180 mm double-column, Helvetica 5-7 pt) ---
MM_TO_IN = 1 / 25.4
FIG_WIDTH_IN = 180 * MM_TO_IN
FIG_HEIGHT_IN = 70 * MM_TO_IN
rcParams['font.family'] = 'Helvetica'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.size'] = 6
rcParams['axes.labelsize'] = 7
rcParams['axes.titlesize'] = 7
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6
rcParams['legend.fontsize'] = 6
rcParams['axes.linewidth'] = 0.5
rcParams['xtick.major.width'] = 0.5
rcParams['ytick.major.width'] = 0.5
rcParams['xtick.major.size'] = 2
rcParams['ytick.major.size'] = 2

# --- Paths ---
# Inputs point at the canonical SCENIC cross-species folder (source of truth)
INPUTS = Path('/home/davidbellini/OneDrive/gabbiani/PROJECTS/all_atlases/SCENIC/cross_species_comparison/inputs')
OUTPUTS = Path(__file__).resolve().parents[2] / 'outputs' / 'figure5'
OUTPUTS.mkdir(parents=True, exist_ok=True)

SPECIES_ORDER = ['fly', 'honeybee', 'pharaonis', 'americana', 'bumblebee']
SPECIES_LABELS = {
    'fly': 'D. melanogaster',
    'honeybee': 'A. mellifera',
    'pharaonis': 'M. pharaonis',
    'americana': 'S. americana',
    'bumblebee': 'B. terrestris',
}

REGULON_FILES = {sp: INPUTS / f'{sp}_regulons.csv' for sp in SPECIES_ORDER}
HOMOLOG_FILES = {
    'honeybee': INPUTS / 'honeybee_to_fly_homologs.tsv',
    'pharaonis': INPUTS / 'pharaonis_to_fly_homologs.tsv',
    'americana': INPUTS / 'americana_to_fly_homologs.tsv',
    'bumblebee': INPUTS / 'bumblebee_to_fly_homologs.tsv',
}

PIDENT_THRESHOLD = 50.0


def load_regulon_tfs(path):
    """Extract unique TFs from a pySCENIC regulon CSV (2 junk header rows)."""
    df = pd.read_csv(path, skiprows=[0, 1])
    tfs = df['TF'].dropna().unique().tolist()
    # Seurat converts _ to - in gene names; homolog files use underscores
    tfs = [tf.replace('-', '_') for tf in tfs]
    return tfs


def load_homolog_map(path, species):
    """Build species_gene -> (fly_gene, pident, bitscore) dict from homolog TSV.

    Filters at >=50% pident, takes best hit (highest bitscore) per species gene.
    For pharaonis, builds lookups from both sp_transcript and sp_gene_name columns.
    """
    df = pd.read_csv(path, sep='\t')
    df = df[df['pident'] >= PIDENT_THRESHOLD].copy()

    if species == 'pharaonis':
        # Build two lookup dicts: one by sp_transcript (Pharaoh_ant_*), one by sp_gene_name (fly-style)
        # Both map to fly_gene via BLAST, not by name matching
        result = {}
        for col in ['sp_gene_name', 'sp_transcript']:
            sub = df.dropna(subset=[col, 'fly_gene']).copy()
            sub[col] = sub[col].astype(str).str.replace('-', '_')
            sub['fly_gene'] = sub['fly_gene'].astype(str).str.replace('-', '_')
            sub = sub.sort_values('bitscore', ascending=False).drop_duplicates(subset=[col], keep='first')
            for _, row in sub.iterrows():
                key = row[col]
                if key not in result:  # sp_gene_name takes priority (tried first)
                    result[key] = (row['fly_gene'], row['pident'], row['bitscore'])
        return result
    else:
        col = 'sp_transcript'
        df = df.dropna(subset=[col, 'fly_gene']).copy()
        df[col] = df[col].astype(str).str.replace('-', '_')
        df['fly_gene'] = df['fly_gene'].astype(str).str.replace('-', '_')
        df = df.sort_values('bitscore', ascending=False).drop_duplicates(subset=[col], keep='first')
        return {row[col]: (row['fly_gene'], row['pident'], row['bitscore'])
                for _, row in df.iterrows()}


# Paralogs that are not reliably distinguishable by protein identity are merged
# to a single label so the comparison treats them as one regulon. ey/toy are the
# two ancestral Drosophila Pax6 paralogs; the "ey-locus" hits ey and toy at near-
# equal identity across species, so the ey-vs-toy call is arbitrary (see
# pax6_zscores_qc.py). Merged for Figure 5.
PARALOG_MERGE = {'ey': 'toy/ey', 'toy': 'toy/ey'}


def _canon_fly_gene(g):
    return PARALOG_MERGE.get(g, g)


def map_tfs_to_fly(tfs, homolog_map, species):
    """Map a list of species TFs to fly namespace.

    Returns: (mapped_dict {sp_tf: (fly_gene, pident, bitscore)}, unmapped_list)
    """
    if species == 'fly':
        return {tf: (_canon_fly_gene(tf), 100.0, 0.0) for tf in tfs}, []

    mapped = {}
    unmapped = []
    for tf in tfs:
        if tf in homolog_map:
            fg, pid, bit = homolog_map[tf]
            mapped[tf] = (_canon_fly_gene(fg), pid, bit)
        else:
            unmapped.append(tf)
    return mapped, unmapped


def find_many_to_one(mapped_dict, species):
    """Find cases where multiple species TFs map to the same fly gene."""
    fly_to_sp = defaultdict(list)
    for sp_tf, (fly_gene, pident, bitscore) in mapped_dict.items():
        fly_to_sp[fly_gene].append((sp_tf, pident, bitscore))
    return {fg: tfs for fg, tfs in fly_to_sp.items() if len(tfs) > 1}


def main():
    # --- Load TFs ---
    species_tfs = {}
    for sp in SPECIES_ORDER:
        species_tfs[sp] = load_regulon_tfs(REGULON_FILES[sp])

    # --- Load homolog maps ---
    homolog_maps = {}
    for sp in SPECIES_ORDER:
        if sp == 'fly':
            continue
        homolog_maps[sp] = load_homolog_map(HOMOLOG_FILES[sp], sp)

    # --- Map TFs to fly namespace ---
    mapped = {}      # {species: {sp_tf: (fly_gene, pident, bitscore)}}
    unmapped = {}    # {species: [sp_tf, ...]}
    for sp in SPECIES_ORDER:
        hmap = homolog_maps.get(sp, {})
        mapped[sp], unmapped[sp] = map_tfs_to_fly(species_tfs[sp], hmap, sp)

    # --- Build fly-gene sets (mapped TFs only) ---
    # Unmapped TFs get unique placeholders so they count as species-unique
    fly_sets = {}
    for sp in SPECIES_ORDER:
        s = set()
        for sp_tf, (fly_gene, _, _) in mapped[sp].items():
            s.add(fly_gene)
        for sp_tf in unmapped[sp]:
            s.add(f'_unmapped_{sp}_{sp_tf}')
        fly_sets[sp] = s

    # --- Summary ---
    print("=" * 70)
    print("5-Way Cross-Species Regulon TF Comparison (Fly-Hub, >=50% pident)")
    print("=" * 70)
    for sp in SPECIES_ORDER:
        n_total = len(species_tfs[sp])
        n_mapped = len(mapped[sp])
        n_unmapped = len(unmapped[sp])
        n_fly_unique = len(set(v[0] for v in mapped[sp].values()))
        print(f"\n{SPECIES_LABELS[sp]} ({sp}):")
        print(f"  Total regulon TFs:    {n_total}")
        print(f"  Mapped to fly:        {n_mapped} -> {n_fly_unique} unique fly genes")
        print(f"  Unmapped (no hit):    {n_unmapped}")

        # Many-to-one
        m2o = find_many_to_one(mapped[sp], sp)
        if m2o:
            print(f"  Many-to-one mappings: {len(m2o)} fly genes have multiple {sp} TFs:")
            for fg, tfs in sorted(m2o.items()):
                tf_str = ', '.join(f'{t[0]} ({t[1]:.1f}%)' for t in tfs)
                print(f"    {fg} <- {tf_str}")

    # --- Save mapping table ---
    rows = []
    for sp in SPECIES_ORDER:
        for sp_tf, (fly_gene, pident, bitscore) in mapped[sp].items():
            rows.append({'species': sp, 'species_tf': sp_tf, 'fly_ortholog': fly_gene,
                         'pident': pident, 'bitscore': bitscore})
    pd.DataFrame(rows).to_csv(OUTPUTS / 'tf_mapping_table.csv', index=False)

    # --- Save unmapped TFs ---
    urows = []
    for sp in SPECIES_ORDER:
        for sp_tf in unmapped[sp]:
            urows.append({'species': sp, 'species_tf': sp_tf})
    pd.DataFrame(urows).to_csv(OUTPUTS / 'unmapped_tfs.csv', index=False)

    # --- Save many-to-one report ---
    m2o_rows = []
    for sp in SPECIES_ORDER:
        m2o = find_many_to_one(mapped[sp], sp)
        for fly_gene, tfs in m2o.items():
            for sp_tf, pident, bitscore in tfs:
                m2o_rows.append({'species': sp, 'fly_gene': fly_gene,
                                 'species_tf': sp_tf, 'pident': pident, 'bitscore': bitscore})
    pd.DataFrame(m2o_rows).to_csv(OUTPUTS / 'many_to_one_mappings.csv', index=False)

    # --- Compute intersections ---
    # For the Venn/UpSet, use only mapped fly genes (no placeholders)
    # Unmapped TFs are species-unique by definition
    mapped_sets = {}
    for sp in SPECIES_ORDER:
        mapped_sets[sp] = set(v[0] for v in mapped[sp].values())

    # All fly genes across all species
    all_fly_genes = set()
    for s in mapped_sets.values():
        all_fly_genes |= s

    # For each fly gene, determine which species have it
    gene_memberships = {}
    for gene in all_fly_genes:
        members = frozenset(sp for sp in SPECIES_ORDER if gene in mapped_sets[sp])
        gene_memberships[gene] = members

    # Group by membership pattern
    pattern_genes = defaultdict(list)
    for gene, members in gene_memberships.items():
        pattern_genes[members].append(gene)

    # Save intersection summary
    int_rows = []
    for pattern in sorted(pattern_genes.keys(), key=lambda x: (-len(x), sorted(x))):
        genes = sorted(pattern_genes[pattern])
        species_str = ', '.join(sorted(pattern))
        int_rows.append({
            'n_species': len(pattern),
            'species': species_str,
            'count': len(genes),
            'fly_genes': '; '.join(genes)
        })
    # Add unmapped counts
    for sp in SPECIES_ORDER:
        if unmapped[sp]:
            int_rows.append({
                'n_species': 1,
                'species': f'{sp} (unmapped)',
                'count': len(unmapped[sp]),
                'fly_genes': '; '.join(sorted(unmapped[sp]))
            })
    pd.DataFrame(int_rows).to_csv(OUTPUTS / 'intersection_summary.csv', index=False)

    print(f"\n{'=' * 70}")
    print("Intersection Summary (mapped TFs only):")
    print(f"{'=' * 70}")
    print(f"Total unique fly genes across all species: {len(all_fly_genes)}")
    for pattern in sorted(pattern_genes.keys(), key=lambda x: (-len(x), sorted(x))):
        genes = pattern_genes[pattern]
        sp_str = ' & '.join(sorted(pattern))
        print(f"  {sp_str}: {len(genes)}")
    for sp in SPECIES_ORDER:
        if unmapped[sp]:
            print(f"  {sp} (unmapped, species-unique): {len(unmapped[sp])}")

    # --- UpSet plot (Nature EE, 180 mm width, Helvetica 5-7 pt) ---
    categories = [SPECIES_LABELS[sp] for sp in SPECIES_ORDER]
    index_tuples = [tuple(gene in mapped_sets[sp] for sp in SPECIES_ORDER)
                    for gene in sorted(all_fly_genes)]
    mi = pd.MultiIndex.from_tuples(index_tuples, names=categories)
    upset_data = pd.Series(np.ones(len(index_tuples), dtype=int), index=mi)
    upset_data = upset_data.groupby(level=categories).sum()

    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
    # element_size=None makes UpSet honor figsize instead of forcing a fixed per-column width
    upset = UpSet(upset_data, show_counts=True, sort_by='cardinality',
                  show_percentages=False, facecolor='black',
                  element_size=None,
                  intersection_plot_elements=3, totals_plot_elements=3)
    axes = upset.plot(fig=fig)

    # Helvetica sizes: counts 7pt bold, species / tick labels 5pt regular
    for name, ax in axes.items():
        ax.tick_params(labelsize=5, width=0.5, length=2)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontsize(5)
            lbl.set_fontfamily('Helvetica')
        for txt in ax.texts:
            txt.set_fontfamily('Helvetica')
            if name == 'intersections':
                txt.set_fontsize(7)
                txt.set_fontweight('bold')
            else:
                txt.set_fontsize(6)
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

    # Shrink the matrix dots (default 200 is far too large at 180 mm) and
    # thin the connector lines. PDF supports transparency, so we keep the
    # default alpha-based inactive dots and alternating-row shading.
    if 'matrix' in axes:
        from matplotlib.collections import LineCollection, PathCollection
        for coll in axes['matrix'].collections:
            if isinstance(coll, PathCollection):
                coll.set_sizes([25])
                coll.set_linewidths(0.3)
            elif isinstance(coll, LineCollection):
                coll.set_linewidth(0.5)

    # No bbox_inches='tight' — that re-expands the canvas; we want a hard 180 mm width
    fig.savefig(OUTPUTS / 'regulon_upset.pdf', format='pdf')
    fig.savefig(OUTPUTS / 'regulon_upset.png', format='png', dpi=600)
    plt.close(fig)

    print(f"\nOutputs saved to: {OUTPUTS}/")
    print("  - tf_mapping_table.csv")
    print("  - unmapped_tfs.csv")
    print("  - many_to_one_mappings.csv")
    print("  - intersection_summary.csv")
    print("  - regulon_upset.pdf")
    print("  - regulon_upset.png")


if __name__ == '__main__':
    main()
