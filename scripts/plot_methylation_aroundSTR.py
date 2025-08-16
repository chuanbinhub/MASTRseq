#MASTR-seq STR containing gene promoter region DNA methylation plotting script <plot_methylation_aroundSTR.py>
#Written by Han-Seul Ryu (Last modified Sept 16th 2024)

#Modified 8/5/25 by Chuanbin Su to integrate STR calling scripts into snakemake container 

#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import os
import sys

def main():

    if len(sys.argv) != 7:
        print("Usage: script.py [output_plot] [chrom] [start] [end] [threshold] [input_tsv]")
        sys.exit(1)

    output = sys.argv[1]
    plot_chrom = sys.argv[2]
    plot_start = int(sys.argv[3])
    plot_end = int(sys.argv[4])
    threshold = float(sys.argv[5])
    tsv_file = sys.argv[6]

    # Load tsv
    try:
        df = pd.read_table(tsv_file, sep=r"\s+")
    except Exception as e:
        print(f"Error reading {tsv_file}: {e}")
        sys.exit(1)

    # Filter to methylated cytosines only
    df = df[df['mod_code'] == 'm'][['read_id', 'ref_position', 'chrom', 'mod_qual']]
    df_filter = df[
        (df['chrom'] == plot_chrom) &
        (df['ref_position'] >= plot_start) &
        (df['ref_position'] <= plot_end)
    ]

    if df_filter.empty:
        print(f"No mCpGs found in {plot_chrom}:{plot_start}-{plot_end}. Skipping plot.")
        sys.exit(0)

    # Build dictionary of read_id [True, False, ...] for methylation status
    store = {}
    for _, row in df_filter.iterrows():
        store.setdefault(row['read_id'], []).append(row['mod_qual'] >= threshold)

    # Calculate methylation fraction per read
    fractions = [sum(v) / len(v) for v in store.values() if len(v) > 0]

    if not fractions:
        print("No reads passed methylation filtering. Exiting.")
        sys.exit(0)

    # Plot
    fig, ax = plt.subplots(figsize=(6.5, 5))
    sns.kdeplot(fractions, ax=ax, color='darkblue', fill=True, linewidth=2)

    ax.set_xlabel(f"Fraction of CpGs methylated in {plot_chrom}:{plot_start}-{plot_end}", fontsize=12)
    ax.set_ylabel("Read density", fontsize=12)
    ax.set_xlim(0, 1)
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(output)

    print(f"\n✅ Plot saved: {output}")
    print(f"Included reads: {len(fractions)}")

if __name__ == '__main__':
    main()