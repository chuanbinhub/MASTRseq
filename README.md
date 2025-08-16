# MASTR-seq Snakemake Pipeline

This Snakemake pipeline supports short tandem repeats (STR) counting, methylation profiling analysis from nanopore long-read sequencing data.

# Installation
### Conda Environment Setup
Snakemake requires Conda version 24.7.1 or later.  
Create a conda environment using the following `env.yaml`:

### Pre-installation
Currently, the ont-modkit package (version 0.5.0), a bioinformatics tool for analyzing modified bases from Oxford Nanopore, is distributed via Bioconda for the following platforms: linux-64, osx-64 (Intel-based macOS), and linux-aarch64. There is no osx-arm64 (Apple Silicon) build available at this time.
For macOS Apple Silicon users, modkit can be installed natively using the Cargo tool:

```bash
cargo install --git https://github.com/nanoporetech/modkit.git
```
```yaml
name: mastrseq
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.12
  - pandas=2.2.2
  - matplotlib=3.9.2
  - pysam=0.22.1
  - seaborn=0.13.2
  - samtools=1.20
  - numpy
  - snakemake
  - graphviz
  - ont-modkit=0.5
```

### Run the conda command to create the environment:

```bash
conda env create -f env.yaml
```

### Activate the environment:

```bash
conda activate mastrseq
```


# Usage
### General Configuration

```python
# Load pipeline configuration
configfile: "config.yaml"

SAMPLES = config["samples"]
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
STR_type = config["str_type"]
REF_DIR = config["ref_dir"]
ENV = "env.yaml"
```
**Note**: The `methylation_inSTR` and `methylation_inSTR_plot` rules are excluded for the HTT STR type.

Set up your `config.yaml` by adjusting:
sample names, 
input paths,
output paths, 
reference directory: path to the reference genome, e.g., ``` reference/hg39.fa ``` 
STR type (e.g., `HTT`, `FMR1`, or `C9orf72`), 
STR seq motif,
methylation threshold, 
mutation length threshold. Then run the snakemake pipeline using the following commands:

```bash
snakemake --cores <NUM_CORES>
```
<NUM_CORES> specifies the number of CPU cores for alignment.

If you run sucessfully, you will see output files in <results>. look like this:
```
{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt
{OUTPUT_DIR}/str_plot/{{sample}}_str_plot.pdf
{OUTPUT_DIR}/methylation_aroundSTR/{{sample}}_filtered_extractCpG.tsv
{OUTPUT_DIR}/methylation_aroundSTR/{{sample}}_filtered.bam
{OUTPUT_DIR}/methylation_aroundSTR/{{sample}}_filtered.bam.bai
{OUTPUT_DIR}/methylation_plot/{{sample}}_methylation_density_{STR type}.pdf

{OUTPUT_DIR}/methylation_inSTR/{{sample}}_methylation_inSTR_density_plot.pdf
{OUTPUT_DIR}/methylation_inSTR/{{sample}}/STRlength_5mC-likelihood.tsv
{OUTPUT_DIR}/methylation_inSTR/{{sample}}/plots_for_each_read/*.png

{OUTPUT_DIR}/logs/str_count/{{sample}}.log
{OUTPUT_DIR}/log/str_count_plot/{{sample}}.log
{OUTPUT_DIR}/logs/methylation_aroundSTR/{{sample}}.log
{OUTPUT_DIR}/logs/methylation_density_plot/{{sample}}.log
{OUTPUT_DIR}/logs/methylation_inSTR/{{sample}}.log
{OUTPUT_DIR}/logs/methylation_inSTR_plot/{{sample}}.log
```









