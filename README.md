# MASTR-seq Snakemake Pipeline

This Snakemake pipeline supports STR counting, methylation profiling analysis for the target nanopore dataset.

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

### Rule: all

```python
rule all:
    input:
        expand(f"{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/str_plot/{{sample}}_str_plot.pdf", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/methylation_out/{{sample}}_flank_methyl.tsv", sample=SAMPLES)
```

### Rule: get_STRbam

```python
rule get_STRbam:
    input:
        bam = f"{INPUT_DIR}/{{sample}}.bam"
    output:
        f"{OUTPUT_DIR}/bam_by_STRtype/{{sample}}_{{STR_type}}.bam"
    log:
        f"{OUTPUT_DIR}/bam_by_STRtype/logs/{{sample}}_{{STR_type}}_extract.log"
    conda:
        ENV
    shell:
        "samtools view -h {input.bam} | "
        "awk '{{if ($1 ~ /^@/ || $0 ~ /STR_MARK/) print $0}}' | "
        "samtools view -b -o {output} - ; "
        "echo 'STR BAM extraction done' > {log}"
```

### Rule: count_STR

```python
rule count_STR:
    input:
        bam = f"{OUTPUT_DIR}/bam_by_STRtype/{{sample}}_{{STR_type}}.bam"
    output:
        f"{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt"
    conda:
        ENV
    script:
        "scripts/str_counter.py"
```

### Rule: plot_STR

```python
rule plot_STR:
    input:
        f"{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt"
    output:
        f"{OUTPUT_DIR}/str_plot/{{sample}}_str_plot.pdf"
    conda:
        ENV
    script:
        "scripts/plot_STR_distribution.py"
```

### Rule: methylation_flank

```python
rule methylation_flank:
    input:
        bam = f"{OUTPUT_DIR}/bam_by_STRtype/{{sample}}_{{STR_type}}.bam"
    output:
        f"{OUTPUT_DIR}/methylation_out/{{sample}}_flank_methyl.tsv"
    conda:
        ENV
    script:
        "scripts/plot_methylation_aroundSTR.py"
```

**Note**: The `methylation_inSTR` and `methylation_inSTR_plot` rules are excluded for the HTT STR type.

# Running the Pipeline

Make sure your `config.yaml` is set up correctly with sample names, input/output paths, reference directory, STR type (e.g., `HTT`, `FMR1`, or `C9orf72`), STR seq motif,methylation threshold and mutation length threshold. Then run:

```bash
snakemake --cores <NUM_CORES>
```
