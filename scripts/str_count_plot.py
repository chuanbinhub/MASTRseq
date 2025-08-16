#Writen by Chuanbin Su 
#Date: 08/01/2025

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python str_count_plot_bin_100.py <input_txt> <output_pdf>", file=sys.stderr)
        sys.exit(1)

    input_txt = sys.argv[1]
    output_pdf = sys.argv[2]

    # Derive sample name from filename
    sample_name = os.path.basename(input_txt).split("_")[0]

    # Load data
    df = pd.read_csv(input_txt, sep="\t", header=None, names=["read_name", "str_count"])
    
    # Set seaborn style
    sns.set(style="whitegrid")

    # Plot histogram
    plt.figure(figsize=(8, 4))
    ax = sns.histplot(df["str_count"], bins=100, kde=False, element="step", fill=True, color="#1f77b4")


    # Determine max repeat count for dynamic scaling
    max_count = int(df["str_count"].max())
    rounded_max = int(np.ceil(max_count / 100.0)) * 100  # round up to the nearest 100 for cleaner axis

    # Create histogram
    counts, bins = np.histogram(df["str_count"], bins=100, range=(0, rounded_max))

    # Plot
    ax.set_xlim(0, rounded_max)
    ax.set_xticks(np.arange(0, rounded_max + 1, 100))
    ax.set_xlabel("Repeat Count")
    ax.set_ylabel("Read Count")
    ax.set_title(f"Repeat Count Distribution ({sample_name})", fontsize=13)
    ax.grid(axis='x', color='lightgray')
    ax.grid(axis='y', visible=False)



    # Vertical reference lines
    ax.axvline(x=45, color='orange', linestyle='--', linewidth=1)
    ax.axvline(x=200, color='red', linestyle='--', linewidth=1)

    # Calculate and annotate max bin count
    counts, bins = np.histogram(df["str_count"], bins=100, range=(0, 900))
    #max_count = counts.max()
    #ax.text(0.02, 0.95, f"max count: {max_count}", transform=ax.transAxes,
    #        fontsize=10, verticalalignment='top', horizontalalignment='left')

    # Save
    plt.tight_layout()
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[✓] Plot saved to {output_pdf}")


if __name__ == "__main__":
    main()