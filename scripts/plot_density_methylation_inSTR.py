'''
Written by Hanseul Ryu (Last updated Sept 19 2024)
Modified by Chuanbin Su to integrate STR calling scripts into snakemake container (Aug 2025)

INPUT: 
(1) STRlength_5mC-likelihood.tsv file generated from run_methylation_inSTR.py, 
(2) methylation-likelihood threshold, 
(3) desired output name, 
(4) (optional) mutation-length threshold

OUTPUT: density plot with fraction of CpGs methylated within STR on the x-axis and density of reads on y-axis
- if mutation-length threshold is provided, 2 different density plots will be plotted: one for mutation-length alleles and one for normal-length alleles
- if mutation-length threshold is NOT provided, all filtered reads will be included in one density plot

NOTE about methylation-likelihood threshold:
    The 'methylation_likelihood' column in your tsv files store the methylation likelihood value from 0 to 1, as computed by Dorado basecalling.
    If specifying 0.6 as threshold for methylation-likelihood, that means regarding any score higher than 0.6 as methylated. Any score below this threshold will be considered "unmethylated."

'''
#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main():


    input_tsv = sys.argv[1]
    methylation_threshold = float(sys.argv[2])
    output_name = sys.argv[3]
    
    df = pd.read_table(input_tsv) 
    df = df [['read_type','name','STR_length','methylation_likelihood']]

    if len(sys.argv) == 5:
        mutation_threshold = int(sys.argv[4])
        
        percent_of_CpGs_methylated_ml = [] 
        percent_of_CpGs_methylated_nl = [] 

        for i, row in df.iterrows():
            string = row['methylation_likelihood']
            if string.strip() in ("[]", ""):
                continue  # skip empty rows
            string_arr = string.lstrip('[').rstrip(']').split(', ')
            try:
                np_arr = np.array(string_arr, dtype=float)
            except ValueError:
               continue  # skip if conversion fails
            bool_arr = np_arr > methylation_threshold
            fraction = sum(bool_arr) / len(bool_arr)


            if int(row['STR_length']) > mutation_threshold:
                percent_of_CpGs_methylated_ml.append(fraction)
            else:
                percent_of_CpGs_methylated_nl.append(fraction)

        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax = sns.kdeplot(percent_of_CpGs_methylated_ml, label = 'ML read')
        ax = sns.kdeplot(percent_of_CpGs_methylated_nl, label = 'NL read')
        ax.set_xlabel("Fraction of CpGs methylated in STR")
        ax.grid(axis='y')
        ax.grid(axis='x')
        plt.legend()
        plt.savefig(output_name)


    else:
        percent_of_CpGs_methylated_all = [] 

        for i, row in df.iterrows():
            string = row['methylation_likelihood']
            string_arr = string.lstrip('[').rstrip(']').split(', ')
            np_arr = np.array(string_arr,dtype=float)
            bool_arr = np_arr > methylation_threshold
            fraction = sum(bool_arr)/len(bool_arr)
            percent_of_CpGs_methylated_all.append(fraction)
        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax = sns.kdeplot(percent_of_CpGs_methylated_all)
        ax.set_xlabel("Fraction of CpGs methylated in STR")
        ax.grid(axis='y')
        ax.grid(axis='x')
        plt.savefig(output_name)

    print(f'\nDensity plots can be found here: {output_name}\n')  
   

    return




if __name__ == '__main__':
    main()