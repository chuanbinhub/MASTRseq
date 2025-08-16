# Load pipeline configuration
configfile: "config.yaml"

SAMPLES = config["samples"]
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
STR_type = config["str_type"] 
REF_DIR = config["ref_dir"]
ENV = "env.yaml"
methylation_threshold = config["methylation_threshold"]
mutation_len_threshold = config["mutation_len_threshold"]


rule all:
    input:
        expand("{output_dir}/str_count/{sample}_allcounts.txt",
               output_dir=OUTPUT_DIR, sample=SAMPLES),
        expand("{output_dir}/str_count/{sample}_fwdcounts.txt",
               output_dir=OUTPUT_DIR, sample=SAMPLES),
        expand("{output_dir}/str_count/{sample}_revcounts.txt",
               output_dir=OUTPUT_DIR, sample=SAMPLES),
        expand("{output_dir}/str_plot/{sample}_str_plot.pdf",
               output_dir=OUTPUT_DIR, sample=SAMPLES),
        expand("{output_dir}/methylation_aroundSTR/{sample}_filtered_extractCpG.tsv",
               output_dir=OUTPUT_DIR, sample=SAMPLES),
        expand("{output_dir}/methylation_plot/{sample}_methylation_density_{str_type}.pdf",
               output_dir=OUTPUT_DIR, str_type=STR_type, sample=SAMPLES),
        *([] if STR_type == "HTT" else
          expand("{output_dir}/methylation_inSTR/{sample}/STRlength_5mC-likelihood.tsv",
                 output_dir=OUTPUT_DIR, sample=SAMPLES)),
        *([] if STR_type == "HTT" else
          expand("{output_dir}/methylation_inSTR/{sample}/plots_for_each_read",
                 output_dir=OUTPUT_DIR, sample=SAMPLES)),
        *([] if STR_type == "HTT" else
          expand("{output_dir}/methylation_inSTR/{sample}_methylation_inSTR_density_plot.pdf",
                 output_dir=OUTPUT_DIR, sample=SAMPLES))


rule str_count:
    input:
        bam = lambda wildcards: f"{INPUT_DIR}/{wildcards.sample}.sorted.bam"
    output:
        counts = f"{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt",
        fwd = f"{OUTPUT_DIR}/str_count/{{sample}}_fwdcounts.txt",
        rev = f"{OUTPUT_DIR}/str_count/{{sample}}_revcounts.txt"
    params:
        outdir = f"{OUTPUT_DIR}/str_count",
        str_type = STR_type
    log:
        f"{OUTPUT_DIR}/logs/str_count/{{sample}}.log"
    conda:
        ENV
    shell:
        """
        mkdir -p "{params.outdir}"
        python scripts/run_str_counter.py {params.str_type} "{input.bam}" "{params.outdir}" > "{log}" 2>&1
        """

rule str_count_plot:
    input:
        txt = f"{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt"
    output:
        plot = f"{OUTPUT_DIR}/str_plot/{{sample}}_str_plot.pdf"
    log:
        f"{OUTPUT_DIR}/logs/str_count_plot/{{sample}}.log"
    conda:
        ENV
    shell:
        """
        python scripts/str_count_plot.py "{input.txt}" "{output.plot}" > "{log}" 2>&1
        """

rule methylation_aroundSTR:
    input:
        bam = lambda wildcards: f"{INPUT_DIR}/{wildcards.sample}.sorted.bam",
        counts = f"{OUTPUT_DIR}/str_count/{{sample}}_allcounts.txt"
    output:
        tsv = f"{OUTPUT_DIR}/methylation_aroundSTR/{{sample}}_filtered_extractCpG.tsv"
    params:
        outdir = f"{OUTPUT_DIR}/methylation_aroundSTR",
        ref = f"{REF_DIR}/hg38.fa",
        modkit = "modkit"
    log:
        f"{OUTPUT_DIR}/logs/methylation_aroundSTR/{{sample}}.log"
    conda:
        ENV
    shell:
        """
        mkdir -p "{params.outdir}"
        samtools view -N "{input.counts}" -o "{params.outdir}/{wildcards.sample}_filtered.unsorted.bam" "{input.bam}"
        samtools sort -o "{params.outdir}/{wildcards.sample}_filtered.bam" "{params.outdir}/{wildcards.sample}_filtered.unsorted.bam"
        samtools index "{params.outdir}/{wildcards.sample}_filtered.bam"
        "{params.modkit}" extract full --cpg --reference "{params.ref}" \
            "{params.outdir}/{wildcards.sample}_filtered.bam" \
            "{output.tsv}" > "{log}" 2>&1
        """

# Define region of gene promoter per STR type
def get_region(wildcards):
    str_regions = {
        "FMR1":  {"chrom": "chrX",  "start": 147911419, "end": 147911919 },
        "HTT":   {"chrom": "chr4",  "start": 3074181,   "end": 3074681},
        "C9ORF72": {"chrom": "chr9", "start": 27573866,  "end": 27574866},
    }
    return str_regions[STR_type]  # STR_type is from config.yaml

rule methylation_density_plot:
    input:
        tsv = f"{OUTPUT_DIR}/methylation_aroundSTR/{{sample}}_filtered_extractCpG.tsv"
    output:
        plot = f"{OUTPUT_DIR}/methylation_plot/{{sample}}_methylation_density_{STR_type}.pdf"
    params:
        region = get_region
    conda:
        ENV
    log:
        f"{OUTPUT_DIR}/logs/methylation_density_plot/{{sample}}.log"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/methylation_plot
        python scripts/plot_methylation_aroundSTR.py \
            {output.plot} \
            {params.region[chrom]} {params.region[start]} {params.region[end]} \
            {methylation_threshold} \
            {input.tsv} > {log} 2>&1
        """


def get_bam_by_strtype(wildcards):
    if STR_type == "FMR1":
        return f"{INPUT_DIR}/FXS_R10.4.1_basecalled_sup_wCpG_unaligned.bam"
    elif STR_type == "C9ORF72":
        return f"{INPUT_DIR}/ALS_R10.4.1_basecalled_sup_wCpG_unaligned.bam"
    else:
        raise ValueError(f"Unrecognized STR_type: {STR_type}")


if STR_type != "HTT":

    rule methylation_inSTR:
        input:
            bam = get_bam_by_strtype,
            fwd = f"{OUTPUT_DIR}/str_count/{{sample}}_fwdcounts.txt",
            rev = f"{OUTPUT_DIR}/str_count/{{sample}}_revcounts.txt"
        output:
            tsv = f"{OUTPUT_DIR}/methylation_inSTR/{{sample}}/STRlength_5mC-likelihood.tsv",
            plots_dir = directory(f"{OUTPUT_DIR}/methylation_inSTR/{{sample}}/plots_for_each_read")
        params:
            outdir = f"{OUTPUT_DIR}/methylation_inSTR/{{sample}}",
            str_seq = config["str_seq"]
        log:
            f"{OUTPUT_DIR}/logs/methylation_inSTR/{{sample}}.log"
        conda:
            ENV
        shell:
            """
            mkdir -p {params.outdir}/plots_for_each_read
            python scripts/run_methylation_inSTR.py \
                {input.bam} {input.fwd} {input.rev} \
                {params.str_seq} {params.outdir} > {log} 2>&1
            """

    rule methylation_inSTR_plot:
        input:
            tsv = f"{OUTPUT_DIR}/methylation_inSTR/{{sample}}/STRlength_5mC-likelihood.tsv"
        output:
            plot = f"{OUTPUT_DIR}/methylation_inSTR/{{sample}}_methylation_inSTR_density_plot.pdf"
        params:
            threshold = config["methylation_threshold"],
            mutation_len = config["mutation_len_threshold"]
        log:
            f"{OUTPUT_DIR}/logs/methylation_inSTR_plot/{{sample}}.log"
        conda:
            ENV
        run:
            import os
            os.makedirs(f"{OUTPUT_DIR}/methylation_inSTR", exist_ok=True) 
            cmd = f"python scripts/plot_density_methylation_inSTR.py {input.tsv} {params.threshold} {output.plot} {params.mutation_len} > {log} 2>&1"
            shell(cmd)
