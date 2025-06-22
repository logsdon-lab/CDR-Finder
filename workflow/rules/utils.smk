
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
        col_other="0",
    conda:
        "envs/tools.yaml"
    shell:
        """
        {{ join -1 1 -2 1 <(sort -k 1,1 {input.bed}) <(sort -k 1,1 {input.target_bed}) | \
        awk -v OFS="\\t" '{{
            col_add=({params.col_add} == 0) ? 0 : {params.col_add};
            st = $2 + col_add;
            end = $3 + col_add;
            print $1, st, end {params.col_other}
        }}' | \
        sort -k 1,1 -k 2,2n | \
        uniq ;}} > {output} 2> {log}
        """
