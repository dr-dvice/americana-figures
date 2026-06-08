#!/usr/bin/env python
"""
Figure 5B (final): per-family mean per-cell AUCell for three regulons
(TfAP-2, grn, toy) across 5 species.

One row × three panels, max width 180 mm, Helvetica 5–7 pt, EPS + PNG.

Data flow:
  1. For each species: load per-cell AUC (*_auc_matrix.loom) + cell-type labels
     (from h5ad or cluster_names.csv mapping).
  2. For each fly TF (TfAP-2, grn, toy) × species, pick the species-specific
     regulon that has the most unique target genes among those mapped to the
     fly ortholog (parsed from the ctx *_regulons.csv).
  3. Map cell-type labels → harmonised family via harmonization_fig5.csv
     (excludes honeybee cluster 21 "L1/T1/Lai" — mixed-family label).
  4. Mean AUC per (species × family) is the integral ∫ x p̂(x) dx — equivalently
     the sample mean (see build_fig5b_final.md).

Run:
    conda activate pyscenic
    python build_fig5b.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import loompy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ── Paths ──
BASE = Path(__file__).resolve().parent
PAPER = BASE.parent.parent.parent
SCENIC = Path("/home/davidbellini/OneDrive/gabbiani/PROJECTS/all_atlases/SCENIC")
OUT = PAPER / "figures/outputs/figure5"
OUT.mkdir(parents=True, exist_ok=True)

MAPPING_TABLE = SCENIC / "cross_species_comparison/outputs/tf_mapping_table.csv"
HARMONIZATION = BASE / "harmonization_fig5.csv"

# ── Figure contents ──
FIG_TFS = [("TfAP_2", "TfAP-2"), ("grn", "grn"), ("toy", "toy/ey")]

# Pax6 (toy/ey) panel: the ey-vs-toy paralog call is ambiguous (near-equal protein
# identity), so the regulon shown per species is fixed explicitly rather than auto-
# selected. Prior picks are kept, EXCEPT honeybee: its toy-locus regulon
# (LOC102654776) has only 2 targets and is flat across cell types, so it is replaced
# by its informative ey-locus regulon LOC724238.
PAX6_PANEL_REGULON = {
    "fly": "toy",
    "americana": "egapxtmp_018601",
    "bumblebee": "LOC100647828",
    "honeybee": "LOC724238",
    "pharaonis": "toy",
}

SP_ORDER = ["fly", "americana", "bumblebee", "honeybee", "pharaonis"]
SP_LABEL = {"fly": "D. melanogaster", "americana": "S. americana",
            "bumblebee": "B. terrestris", "honeybee": "A. mellifera",
            "pharaonis": "M. pharaonis"}
SPECIES_MARKERS = {
    "fly":        ("o", "#1f77b4"),
    "americana":  ("s", "#d62728"),
    "bumblebee":  ("^", "#2ca02c"),
    "honeybee":   ("D", "#ff7f0e"),
    "pharaonis":  ("v", "#9467bd"),
}
FAMILIES = ["Photoreceptors", "Lamina", "T neurons", "Medulla",
            "Transmedullary", "Distal medulla", "Centrifugal", "Glia"]

EXCLUDED_LABELS = {"honeybee": {"L1/T1/Lai"}}   # mixed-family atlas cluster

# ── Species config ──
CFG = {
    "fly": dict(
        h5ad="/home/davidbellini/OneDrive/gabbiani/SAMap_fly_datasets/fly_atlas/zipursky/fly_ol_2020.h5ad",
        loom=SCENIC/"zipursky/outputs/fly_auc_matrix.loom",
        regulons=SCENIC/"zipursky/outputs/fly_regulons.csv",
        cluster_col="subtype", cluster_names=None, harm_lab_col="fl_lab",
    ),
    "americana": dict(
        h5ad="/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.h5ad",
        loom=Path("/home/davidbellini/files/workspace/SCENIC_nextflow/nextflow/results/step8_aucell/auc_matrix.loom"),
        regulons=Path("/home/davidbellini/files/workspace/SCENIC_nextflow/nextflow/results/step6_ctx/regulons.csv"),
        cluster_col="cell_type", cluster_names=None, harm_lab_col="sa_lab",
    ),
    "bumblebee": dict(
        h5ad="/home/davidbellini/OneDrive/gabbiani/PROJECTS/bumblebee_atlas/bumblebee_AC_v1.h5ad",
        loom=SCENIC/"bumblebee/outputs/bumblebee_auc_matrix.loom",
        regulons=SCENIC/"bumblebee/outputs/bumblebee_regulons.csv",
        cluster_col="seurat_clusters",
        cluster_names=SCENIC/"bumblebee/inputs/cluster_names.csv",
        harm_lab_col="bt_lab",
    ),
    "honeybee": dict(
        h5ad="/home/davidbellini/OneDrive/gabbiani/PROJECTS/honeybee_atlas/honeybee_ol_subset_ldn1.h5ad",
        loom=SCENIC/"honeybee/outputs/honeybee_auc_matrix.loom",
        regulons=SCENIC/"honeybee/outputs/honeybee_regulons.csv",
        cluster_col="seurat_clusters",
        cluster_names=SCENIC/"honeybee/inputs/cluster_names.csv",
        harm_lab_col="hb_lab",
    ),
    "pharaonis": dict(
        h5ad="/home/davidbellini/OneDrive/gabbiani/PROJECTS/pharaonis/pharaonis_male_subclusters.h5ad",
        loom=SCENIC/"pharaonis/outputs/pharaonis_auc_matrix.loom",
        regulons=SCENIC/"pharaonis/outputs/pharaonis_regulons.csv",
        cluster_col="seurat_clusters",
        cluster_names=SCENIC/"pharaonis/inputs/cluster_names.csv",
        harm_lab_col="mp_lab",
    ),
}

# ── Publication style (Nature Ecol Evol: Helvetica, 5–7 pt, max 180 mm) ──
MM_PER_IN = 25.4
FIG_WIDTH_MM = 180
rcParams["font.family"] = "Helvetica"
rcParams["font.size"] = 7
rcParams["axes.titlesize"] = 7
rcParams["axes.labelsize"] = 7
rcParams["xtick.labelsize"] = 6
rcParams["ytick.labelsize"] = 6
rcParams["legend.fontsize"] = 6
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["axes.linewidth"] = 0.5
rcParams["xtick.major.width"] = 0.5
rcParams["ytick.major.width"] = 0.5
rcParams["xtick.major.size"] = 2
rcParams["ytick.major.size"] = 2


# ─────────────────────────────────────────────────────────────────────────
# Data-loading helpers (self-contained; archive-safe)
# ─────────────────────────────────────────────────────────────────────────
TARGET_PATTERN = re.compile(r"'([^']+)',\s*[-\d.eE+]+")


def count_targets_per_tf(regulons_csv):
    df = pd.read_csv(regulons_csv, header=[0, 1, 2])
    def pick(levels):
        for lv in levels[::-1]:
            s = str(lv)
            if s and not s.startswith("Unnamed"):
                return s
        return str(levels[-1])
    df.columns = [pick(c) for c in df.columns]
    counts = {}
    for tf, sub in df.groupby("TF"):
        genes = set()
        for tg in sub["TargetGenes"].dropna():
            genes.update(TARGET_PATTERN.findall(str(tg)))
        counts[str(tf).replace("-", "_")] = len(genes)
    return counts


def load_fly_to_sp_tfs():
    df = pd.read_csv(MAPPING_TABLE)
    out = {}
    for (sp, fg), sub in df.groupby(["species", "fly_ortholog"]):
        tfs = [str(t).replace("-", "_") for t in sub["species_tf"].tolist()]
        out[(sp, fg)] = sorted(set(tfs))
    return out


def load_auc_df(loom_path):
    with loompy.connect(str(loom_path), mode="r", validate=False) as ds:
        auc = ds.ca.RegulonsAUC
        cells = ds.ca.CellID
    cols_raw = list(auc.dtype.names)
    data = np.stack([auc[c] for c in cols_raw], axis=1)
    cols = [c[:-3].replace("-", "_") if c.endswith("(+)")
            else c.replace("-", "_") for c in cols_raw]
    return pd.DataFrame(data, index=list(cells), columns=cols)


def load_cluster_labels(sp):
    cfg = CFG[sp]
    a = ad.read_h5ad(cfg["h5ad"], backed="r")
    raw = a.obs[cfg["cluster_col"]].astype(str)
    if cfg["cluster_names"] is not None:
        names = pd.read_csv(cfg["cluster_names"])
        name_col = [c for c in names.columns if c != "cluster"][0]
        name_map = dict(zip(names["cluster"].astype(str), names[name_col].astype(str)))
        labels = raw.map(name_map).fillna(raw)
    else:
        labels = raw
    labels.index = a.obs_names
    return labels


def tokens(s):
    return set(re.findall(r"[A-Za-z][A-Za-z0-9]*", str(s)))


def build_label_to_family(harm_df, harm_lab_col):
    fam_map = {}
    for _, r in harm_df.iterrows():
        lab_str = r[harm_lab_col]
        fam = r["family"]
        if not isinstance(lab_str, str) or not lab_str.strip():
            continue
        for sl in (s.strip() for s in lab_str.split(";")):
            if sl:
                fam_map[sl] = fam
    return fam_map


def assign_family(cluster_label, fam_exact):
    if cluster_label in fam_exact:
        return fam_exact[cluster_label]
    tks = tokens(cluster_label)
    if not tks:
        return "Other"
    for lbl, fam in fam_exact.items():
        if tokens(lbl) & tks:
            return fam
    return "Other"


# ─────────────────────────────────────────────────────────────────────────
# Pipeline
# ─────────────────────────────────────────────────────────────────────────
def select_regulons(fly_tfs, fly_to_sp_tfs, target_counts, auc_regs_by_sp):
    chosen = {}
    rows = []
    for sp in SP_ORDER:
        avail = auc_regs_by_sp[sp]
        tcnt = target_counts[sp]
        for fg in fly_tfs:
            cands = [t for t in fly_to_sp_tfs.get((sp, fg), []) if t in avail]
            if not cands:
                rows.append(dict(species=sp, fly_tf=fg, chosen=None,
                                 n_targets=0, n_candidates=0))
                continue
            pick = max(cands, key=lambda t: (tcnt.get(t, 0), -len(t)))
            chosen[(sp, fg)] = pick
            rows.append(dict(species=sp, fly_tf=fg, chosen=pick,
                             n_targets=tcnt.get(pick, 0),
                             n_candidates=len(cands)))
    return chosen, pd.DataFrame(rows)


def compute_means(fly_tfs, chosen, auc_by_sp, family_by_sp, labels_by_sp):
    """Per (species × fly_tf × family): compute mean/median plus two
    species-level normalisations of those values:
      z_mean_sd   = (family_mean   - species_mean)   / species_sd
      z_med_mad   = (family_median - species_median) / species_mad

    species-level baselines are computed across ALL kept cells of that species
    for the *specific regulon* (per-regulon baseline, since each panel is one TF).
    """
    rows = []
    for sp in SP_ORDER:
        fam = family_by_sp[sp]
        labels = labels_by_sp[sp]
        excl = EXCLUDED_LABELS.get(sp, set())
        keep = ~labels.isin(excl) if excl else pd.Series(True, index=labels.index)
        fam_keep = fam[keep]
        for fg in fly_tfs:
            reg = chosen.get((sp, fg))
            if reg is None:
                continue
            auc = auc_by_sp[sp][reg]
            common = auc.index.intersection(fam_keep.index)
            a = auc.loc[common].values
            f = fam_keep.loc[common].values

            sp_mean = float(np.mean(a))
            sp_sd = float(np.std(a, ddof=1))
            sp_med = float(np.median(a))
            sp_mad = float(np.median(np.abs(a - sp_med)))  # raw (unscaled) MAD

            for fam_name in FAMILIES:
                mask = (f == fam_name)
                if mask.sum() < 5:
                    continue
                fm_mean = float(a[mask].mean())
                fm_med = float(np.median(a[mask]))
                z1 = (fm_mean - sp_mean) / sp_sd if sp_sd > 0 else np.nan
                z2 = (fm_med - sp_med) / sp_mad if sp_mad > 0 else np.nan
                rows.append(dict(species=sp, fly_tf=fg, family=fam_name,
                                 mean_auc=fm_mean, median_auc=fm_med,
                                 species_mean=sp_mean, species_sd=sp_sd,
                                 species_median=sp_med, species_mad=sp_mad,
                                 z_mean_sd=z1, z_median_mad=z2,
                                 n_cells=int(mask.sum())))
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────────
def plot_final(df, metric, ylabel, out_stem):
    n_panels = len(FIG_TFS)
    # 180 mm total, allow a bit of padding
    fig_w_in = FIG_WIDTH_MM / MM_PER_IN
    fig_h_in = 2.1  # balanced aspect for 1×3 row
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_w_in, fig_h_in),
                             sharex=True)

    n_sp = len(SP_ORDER)
    spread = 0.16 * (n_sp - 1)
    jitter = {sp: -spread/2 + i * (spread / max(n_sp - 1, 1))
              for i, sp in enumerate(SP_ORDER)}
    xticks = np.arange(len(FAMILIES))

    for ax, (fg_key, fg_disp) in zip(axes, FIG_TFS):
        # Alternating column shading + boundary lines so families are visually distinct
        for i in range(len(FAMILIES)):
            if i % 2 == 1:
                ax.axvspan(i - 0.5, i + 0.5, color="#f2f2f2",
                           zorder=0, linewidth=0)
        for i in range(1, len(FAMILIES)):
            ax.axvline(i - 0.5, color="#bbbbbb", lw=0.3, zorder=0.5)

        sub = df[df["fly_tf"] == fg_key]
        for sp in SP_ORDER:
            marker, color = SPECIES_MARKERS[sp]
            sp_sub = sub[sub["species"] == sp].set_index("family")
            xs, ys = [], []
            for i, fam in enumerate(FAMILIES):
                if fam in sp_sub.index and np.isfinite(sp_sub.loc[fam, metric]):
                    xs.append(i + jitter[sp])
                    ys.append(sp_sub.loc[fam, metric])
            if xs:
                ax.scatter(xs, ys, marker=marker, s=18, color=color,
                           edgecolors="black", linewidths=0.35, alpha=0.9,
                           zorder=3, label=SP_LABEL[sp])
        ax.set_title(fg_disp, fontsize=7, pad=2, style="italic")
        ax.set_xticks(xticks)
        ax.set_xticklabels(FAMILIES, rotation=40, ha="right")
        ax.set_xlim(-0.5, len(FAMILIES) - 0.5)
        ax.axhline(0, color="#cccccc", lw=0.4, zorder=1)
        ax.grid(axis="y", color="#eeeeee", lw=0.4, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="both", pad=1.5)

    axes[0].set_ylabel(ylabel, fontsize=7)

    # Shared legend on the far right (outside the last panel)
    handles = [plt.Line2D([0], [0], marker=SPECIES_MARKERS[sp][0],
                          color="white", markerfacecolor=SPECIES_MARKERS[sp][1],
                          markeredgecolor="black", markersize=5,
                          markeredgewidth=0.35, label=SP_LABEL[sp])
               for sp in SP_ORDER]
    axes[-1].legend(handles=handles, loc="upper left",
                    bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=6,
                    handletextpad=0.4, borderaxespad=0)

    fig.tight_layout()
    fig.subplots_adjust(right=0.82, wspace=0.22)

    # Verify width
    w_mm = fig.get_figwidth() * MM_PER_IN
    print(f"figure width: {w_mm:.2f} mm (limit 180 mm)", flush=True)

    fig.savefig(OUT / f"{out_stem}.eps", bbox_inches=None)  # preserve set width
    fig.savefig(OUT / f"{out_stem}.png", dpi=600, bbox_inches=None)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
def main():
    fly_tfs = [k for k, _ in FIG_TFS]

    print("Loading TF mapping + target counts...", flush=True)
    fly_to_sp_tfs = load_fly_to_sp_tfs()
    target_counts = {sp: count_targets_per_tf(cfg["regulons"])
                     for sp, cfg in CFG.items()}

    print("Loading AUC matrices...", flush=True)
    auc_by_sp = {sp: load_auc_df(cfg["loom"]) for sp, cfg in CFG.items()}
    auc_regs_by_sp = {sp: set(df.columns) for sp, df in auc_by_sp.items()}

    print("Loading cluster labels + family assignment...", flush=True)
    harm = pd.read_csv(HARMONIZATION)
    labels_by_sp, family_by_sp = {}, {}
    for sp, cfg in CFG.items():
        lab = load_cluster_labels(sp)
        fam_exact = build_label_to_family(harm, cfg["harm_lab_col"])
        fam = lab.map(lambda x: assign_family(x, fam_exact))
        labels_by_sp[sp] = lab
        family_by_sp[sp] = fam
        print(f"  {sp}: cells={len(lab)} harmonised={(fam != 'Other').sum()}", flush=True)

    print("Selecting regulons...", flush=True)
    chosen, sel_df = select_regulons(fly_tfs, fly_to_sp_tfs,
                                     target_counts, auc_regs_by_sp)

    # Override the Pax6 (toy/ey) panel with the explicit per-species regulon.
    sel_df = sel_df[sel_df["fly_tf"] != "toy"].copy()
    for sp in SP_ORDER:
        reg = PAX6_PANEL_REGULON[sp]
        if reg not in auc_regs_by_sp[sp]:
            raise SystemExit(f"Pax6 panel regulon {reg} not in {sp} AUCell loom")
        chosen[(sp, "toy")] = reg
        sel_df = pd.concat([sel_df, pd.DataFrame([dict(
            species=sp, fly_tf="toy", chosen=reg,
            n_targets=target_counts[sp].get(reg, 0), n_candidates=-1)])],
            ignore_index=True)
    sel_df.to_csv(OUT / "fig5b_regulon_selection.csv", index=False)
    print(sel_df.to_string(index=False), flush=True)

    print("Computing per-family means...", flush=True)
    df = compute_means(fly_tfs, chosen, auc_by_sp, family_by_sp, labels_by_sp)
    df.to_csv(OUT / "fig5b_data.csv", index=False)

    print("Plotting...", flush=True)
    # Variant 1: standard z-score, (family_mean − species_mean) / species_sd
    plot_final(df, metric="z_mean_sd",
               ylabel="z = (mean − μ_sp) / σ_sp",
               out_stem="fig5b_zscore_mean_sd")
    # Variant 2: robust z-score, (family_median − species_median) / species_MAD
    plot_final(df, metric="z_median_mad",
               ylabel="z = (median − med_sp) / MAD_sp",
               out_stem="fig5b_zscore_median_mad")

    # Drop the legacy raw-mean-AUCell outputs if they exist (replaced by z-score variants)
    for stale in ("fig5b.eps", "fig5b.png"):
        p = OUT / stale
        if p.exists():
            p.unlink()

    print(f"\nDone. Outputs in {OUT}/fig5b_zscore_*.{{eps,png}}", flush=True)


if __name__ == "__main__":
    main()
