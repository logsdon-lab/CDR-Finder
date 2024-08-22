from os.path import join


configfile: "config.yaml"


MANIFEST = config["samples"]
WINDOW = config["window_size"]
SLIDE = config["slide"]
REPORT_THRESHOLD = config.get("report_threshold", 25000)
ALR_THRESHOLD = config.get("alr_threshold", 100000)
EDGE_SEARCH = config.get("edge_search", 50000)
LOW_THRESHOLD = config.get("low_threshold", 39)
OUTPUT_DIR = "results"
LOG_DIR = "logs"
BMK_DIR = "benchmarks"


wildcard_constraints:
    sample="|".join(MANIFEST.keys()),


rule all:
    input:
        expand("results/{sample}_CDR.bed", sample=MANIFEST.keys()),


###
# Subset fasta: Extracts subsequence from target bed
###


rule subset_fasta:
    input:
        fasta=lambda wc: MANIFEST[wildcards.sample]["fasta"],
        bed=lambda wc: MANIFEST[wildcards.sample]["regions"],
    output:
        rm_fasta=join(OUTPUT_DIR, "fasta", "{sample}_subset.fasta"),
    log:
        join(LOG_DIR, "subset_fasta_{sample}.log"),
    resources:
        mem=8,
        hrs=1,
    threads: 1
    conda:
        "envs/tools.yaml"
    shell:
        """
        bedtools getfasta -fi {input.fasta} -bed {input.bed} > {output.rm_fasta} 2> {log}
        """


###
# Get methylation information.
###


rule get_methylation_bed:
    input:
        ref=lambda wc: MANIFEST[wildcards.sample]["fasta"],
        bam=lambda wc: MANIFEST[wildcards.sample]["bam"],
    output:
        join(OUTPUT_DIR, "bed", "{sample}_methyl.bed"),
    log:
        join(LOG_DIR, "get_methylation_bed_{sample}.log"),
    benchmark:
        join(BMK_DIR, "get_methylation_bed_{sample}.tsv")
    params:
        output_file_prefix="{sample}",
        methylation="5mC",
    threads: 1
    resources:
        mem=8,
        hrs=1,
    conda:
        "envs/tools.yaml"
    shell:
        """
        modbam2bed \
        -p {params.output_file_prefix} \
        -m {params.methylation} \
        -t {threads} \
        -e --aggregate --cpg \
        {input.ref} {input.bam} > {output.bed} 2> {log}
        """


###
# Subset meth: Intersects target bed with methylation bed
###


rule subset_meth:
    input:
        methylation_tsv=rules.get_methylation_bed.output,
        target_bed=lambda wc: MANIFEST[wildcards.sample]["regions"],
    output:
        subset_bed=join(OUTPUT_DIR, "bed", "{sample}_subset.bed"),
    log:
        join(LOG_DIR, "subset_meth_{sample}.log"),
    benchmark:
        join(BMK_DIR, "subset_meth_{sample}.tsv")
    resources:
        mem=8,
        hrs=12,
    threads: 1
    conda:
        "envs/tools.yaml"
    shell:
        """
        bedtools intersect -a {input.methylation_tsv} -b {input.target_bed} -wa > {output.subset_bed} 2> {log}
        """


###
# run_rm: Runs repeatmasker on the extracted subsequence fasta from rule subset_fasta
###


rule run_rm:
    input:
        fasta=rules.subset_fasta.output.rm_fasta,
    output:
        rm_out=join(
            OUTPUT_DIR,
            "rm",
            "{sample}_subset.fasta.out",
        ),
    log:
        join(LOG_DIR, "run_rm_{sample}.log"),
    benchmark:
        join(BMK_DIR, "run_rm_{sample}.tsv")
    resources:
        mem=8,
        hrs=12,
    threads: 12
    conda:
        "envs/tools.yaml"
    shell:
        """
        RepeatMasker -species human -dir $( dirname {output.rm_out} ) -qq -pa {threads} {input.fasta} 2> {log}
        """


###
# Calculates mean frequency in windows/bins of the methylation tsv over the target region
###


rule calc_windows:
    input:
        script="workflow/scripts/calculate_windows.py",
        methylation_tsv=rules.subset_meth.output.subset_bed,
        target_bed=lambda wc: MANIFEST[wildcards.sample]["regions"],
    output:
        binned_freq=join(OUTPUT_DIR, "bed", "{sample}_binned_freq.bed"),
    log:
        join(LOG_DIR, "calc_windows_{sample}.log"),
    benchmark:
        join(BMK_DIR, "calc_windows_{sample}.tsv")
    resources:
        mem=8,
        hrs=1,
    threads: 1
    conda:
        "envs/python.yaml"
    params:
        window_size=WINDOW,
        slide=SLIDE,
    shell:
        """
        python {input.script} \
        --target_bed {input.target_bed} \
        --methylation_tsv {input.methylation_tsv} \
        --window_size {params.window_size} \
        --slide {params.slide} > {output} 2> {log}
        """


###
# Format RM: Converts Repeatmasker output to bedfile and extracts only ALR annotations
###


rule format_filter_RM:
    input:
        rm_out=rules.run_rm.output.rm_out,
    output:
        rm_bed=join(OUTPUT_DIR, "bed", "{sample}_rm_ALR.bed"),
    log:
        join(LOG_DIR, "format_filter_RM_{sample}.log"),
    resources:
        mem=8,
        hrs=1,
    conda:
        "envs/tools.yaml"
    threads: 1
    shell:
        """
        {{
            awk -v OFS="\\t" '{{
                split($5, ":", chr_names)
                print chr_names[1], $6, $7, $10, $9
            }}' {input.rm_out} | \
            grep "ALR" ;}} > {output.rm_bed} 2> {log}
        """


# ###
# # Filter RM: Returns Repeatmasker bed contig name to original fasta name
# ###


# rule filter_RM:
#     input:
#         repeat_masker=rules.format_RM.output.repeat_bed,
#     output:
#         rm_bed="results/{sample}_rm_ALR.bed",
#     resources:
#         mem=8,
#         hrs=1,
#     threads: 1
#     run:
#         RM_df = pd.read_csv(
#             input.repeat_masker,
#             header=None,
#             sep="\t",
#             names=["chr", "start", "stop", "ALR", "or"],
#         )
#         RM_df["chr"] = RM_df["chr"].str.split(pat=":").str[0]
#         RM_df.to_csv(output.rm_bed, header=None, index=None, sep="\t")


###
# Intersect RM: Merges repeatmasker bed file with 500bp slop, and only keeps those above the ALR Threshold
#               Intersects merged ALR bed file with binned methylation frequency bedfile
###


rule intersect_RM:
    input:
        repeat_masker=rules.format_filter_RM.output.rm_bed,
        binned_freq=rules.calc_windows.output.binned_freq,
    output:
        intersect_bed=join(OUTPUT_DIR, "bed", "{sample}_intersect.bed"),
        merged=join(OUTPUT_DIR, "bed", "{sample}_rm_merged.bed"),
    log:
        join(LOG_DIR, "intersect_RM_{sample}.log"),
    benchmark:
        join(BMK_DIR, "intersect_RM_{sample}.tsv")
    resources:
        mem=8,
        hrs=1,
    params:
        alr_threshold=ALR_THRESHOLD,
    threads: 1
    conda:
        "envs/tools.yaml"
    shell:
        """
        {{ bedtools merge -i {input.repeat_masker} -d 500 | \
        awk '{{if ($3-$2 > {params.alr_threshold}) print}}' ;}} > {output.merged} 2> {log}
        bedtools intersect -a {input.binned_freq} -b {output.merged} -f 1 -wa -u > {output.intersect_bed} 2> {log}
        """


###
# in threshold: Evaluates bed file from intersected bed file and surrounding regions for signatures of CDR
# Initial signature is mean methylation frequency below dynamically calculated threshold
# Confidence is evaluated based on whether or not these low methylation frequency windows clear a size threshold
# + have flanking peaks of high methylation
###


rule in_threshold:
    input:
        script="workflow/scripts/in_threshold.py",
        intersect_bed=rules.intersect_RM.output.intersect_bed,
        target_bed=lambda wc: MANIFEST[wildcards.sample]["regions"],
    output:
        final_call=join(OUTPUT_DIR, "bed", "{sample}_CDR.bed"),
        windows_call=join(OUTPUT_DIR, "bed", "{sample}_windows_CDR.bed"),
    log:
        join(LOG_DIR, "in_threshold_{sample}.log"),
    benchmark:
        join(BMK_DIR, "in_threshold_{sample}.tsv")
    resources:
        mem=8,
        hrs=1,
    params:
        low_threshold=LOW_THRESHOLD,
        report_threshold=REPORT_THRESHOLD,
        edge_search=EDGE_SEARCH,
    conda:
        "envs/python.yaml"
    threads: 1
    shell:
        """
        python {input.script} \
        --intersect_bed {input.intersect_bed} \
        --target_bed {input.target_bed} \
        --final_call {output.final_call} \
        --windows_call {output.windows_call} \
        --low_threshold {params.low_threshold} \
        --report_threshold {params.report_threshold} \
        --edge_search {params.edge_search} 2> {log}
        """
