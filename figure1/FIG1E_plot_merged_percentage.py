#!/usr/bin/env python3
"""
FIG 1C
Plot percentage of gene pairs classified as MERGED, per composite score interval.
Nature Ecology and Evolution — Helvetica, 90x60mm, EPS + PNG.
"""

import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['ps.fonttype'] = 3
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

PANEL_W_MM = 90
PANEL_H_MM = 60
MM_TO_IN = 1 / 25.4

DATA_FILE = '/home/davidbellini/OneDrive/gabbiani/TOOLS/Gene-Misannotation-Tool/annotation_comparison/results/annotation_comparison.csv'
OUT_DIR = '/home/davidbellini/OneDrive/gabbiani/All_paperwork/PAPERS/americana/figures/outputs/figure1/'
PREFIX = 'FIG1_merged_percentage'


def save_eps_png(fig, prefix):
    os.makedirs(OUT_DIR, exist_ok=True)
    eps_path = os.path.join(OUT_DIR, f'{prefix}.eps')
    png_path = os.path.join(OUT_DIR, f'{prefix}.png')
    fig.savefig(eps_path, format='eps', bbox_inches='tight')
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    print(f'Saved: {eps_path}')
    print(f'Saved: {png_path}')


def plot_merged_percentage():
    df = pd.read_csv(DATA_FILE)

    bins = [50, 60, 70, 80, 90, 101]
    labels = ['50-59', '60-69', '70-79', '80-89', '90+']
    df['score_interval'] = pd.cut(df['composite_score'], bins=bins,
                                  labels=labels, include_lowest=True)

    totals = df.groupby('score_interval').size()
    merged = df[df['status'] == 'MERGED'].groupby('score_interval').size()
    pct = (merged.reindex(labels, fill_value=0) / totals.reindex(labels).replace(0, pd.NA) * 100).fillna(0)

    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM_TO_IN, PANEL_H_MM * MM_TO_IN))
    bars = ax.bar(labels, pct.values, color='#4682B4',
                  edgecolor='black', linewidth=0.4, width=0.7)

    # Value labels above bars: pct% (merged/total)
    for i, bar in enumerate(bars):
        h = bar.get_height()
        merged_n = int(merged.reindex(labels, fill_value=0).iloc[i])
        total_n = int(totals.reindex(labels, fill_value=0).iloc[i])
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.8,
                f'{h:.1f}%\n({merged_n}/{total_n})',
                ha='center', va='bottom', fontsize=5)

    ax.set_xlabel('Composite Score Interval', fontsize=6)
    ax.set_ylabel('Merged Gene Pairs (%)', fontsize=6)
    ax.tick_params(axis='both', labelsize=5, length=2, pad=2)
    ax.set_ylim(0, max(pct.values) * 1.40 if pct.values.max() > 0 else 1)

    ax.yaxis.grid(True, alpha=0.3, linestyle='--', linewidth=0.4)
    ax.xaxis.grid(False)
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.text(0.02, 0.98, f'N = {len(df)}', transform=ax.transAxes,
            ha='left', va='top', fontsize=5)

    plt.subplots_adjust(left=0.16, right=0.97, top=0.95, bottom=0.20)

    save_eps_png(fig, PREFIX)
    plt.close(fig)

    print('\nMerge rate by interval:')
    for lab in labels:
        m = int(merged.reindex(labels, fill_value=0).loc[lab])
        t = int(totals.reindex(labels, fill_value=0).loc[lab])
        p = pct.loc[lab]
        print(f'  {lab}: {m}/{t} = {p:.1f}%')


if __name__ == '__main__':
    plot_merged_percentage()
