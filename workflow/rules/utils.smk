
rule add_target_bed_coords:
    input:
        target_bed="",
        bed="",
    output:
        bed="",
    log:
        "",
    resources:
        mem=8,
        hrs=1,
    params:
        col_st="$2",
        col_end="$3",
        col_add="$4",
        col_other="",
    conda:
        "envs/tools.yaml"
    shell:
        """
        {{ join -1 1 -2 1 <(sort -k 1 {input.bed}) <(sort -k 1 {input.target_bed}) | \
        awk -v OFS="\\t" '{{
            print $1, $2 + {params.col_add}, $3 + {params.col_add} {params.col_other}
        }}' | \
        sort -k 2n -k 1 ;}} > {output} 2> {log}
        """
