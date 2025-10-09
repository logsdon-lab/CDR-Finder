# CDR-Finder
[![CI](https://github.com/koisland/CDR-Finder/actions/workflows/main.yaml/badge.svg)](https://github.com/koisland/CDR-Finder/actions/workflows/main.yaml)

This repository contains a Snakemake workflow to identify and quantify hypomethylated regions within centromeres, or Centromere Dip Regions (CDRs; Altemose et al., Science, 2022).

![](docs/chr8.png)

Original pipeline constructed by @arozanski97 with help from Glennis Logsdon

Adapted and distributed by @fkmastrorosa, @wharvey31, and @koisland.

This is done by:
- Extracting the sequence of interest.
- Intersect the sequence with its methylation data.
- Bin the region of interest into 5 kbp windows and calculate the mean methylation percentage.
- Run RepeatMasker on the sequence of interest to identify regions containing alpha-satellite repeats (ALR/Alpha)
- For each alpha-satellite containing region:
    * Merge consecutive bins.
    * Detect valleys in average methylation percentage based on relative prominence.
    * Check sides of each found CDR, filtering other CDRs, to calculate its relative height.
    * Optionally, extend edges to mean signal and some number of standard deviations to get missed peaks.
    * Optionally, merge adjacent detected CDRs.


## Getting Started
```bash
git clone https://github.com/EichlerLab/CDR-Finder.git
cd CDR-Finder
```

To setup with `conda`:
```bash
conda env create --name cdr_finder -f env.yaml
conda activate cdr_finder
```

If using `singularity`, only `snakemake` is required. This will create a python virtual environment and install `snakemake`.
```bash
python -m venv venv && source venv/bin/activate
pip install snakemake
```

## Usage
To run with `conda`.
```bash
snakemake -np --sdm conda -c 4
```

To run with `singularity`. This will use [`logsdonlab/cdr-finder`](https://hub.docker.com/r/logsdonlab/cdr-finder) to run the workflow.
```bash
snakemake -np --sdm apptainer conda  -c 4
```

Or alternatively to add as a workflow to an existing Snakefile.
```bash
git submodule add https://github.com/koisland/CDR-Finder
```

```python
# Pass CDR config here.
CDR_CONFIG = {}

module CDRFinder:
    snakefile:
        "CDR-Finder/workflow/Snakefile"
    config:
        CDR_CONFIG

use rule * from CDRFinder as cdr_*

rule all:
    input:
        rules.cdr_all.input
```

## Input
Multiple samples can be provided via the configfile. Each sample should contain the following:
- `fasta`
    * Genome assembly.
- `bamfile`
    * Alignment BAM file with methylation tags.
    * Requires aligned reads with the Mm and Ml tags (MM and ML also supported), and the reference sequence used for alignment.
    * PacBio data requires that the alignment does not hard-clip supplementary alignment using `pbmm2` or -Y flag for `minimap2`/`winnowmap`.
    * https://github.com/epi2me-labs/modbam2bed?tab=readme-ov-file#usage
- `regions`
    * BED file of target region coordinates.
- `override_chrom_params`
    * Optional
    * JSON file with chrom specific parameters that override defaults specified. See `python workflow/scripts/cdr_finder.py -h` for parameter list.
- `repeatmasker`
    * Optional
    * Repeatmasker BED4 file in absolute coordinates where 4th column in repeat type (ex. `ALR/Alpha`).
    * If provided, avoids running RepeatMasker.

## Output
- `cdr_bed`
    * CDR regions.
    * `{config.output_dir}/bed/{sample}_CDR.bed`
- `cdr_plot`
    * CDR regions plotted with RepeatMasker annotations.
    * `{config.output_dir}/plot/{sample}/{ctg}.png`
    * To change the layout or text of the plot, modify the [`cenplot`](https://github.com/logsdon-lab/cenplot) track TOML file.:
        * `{config.output_dir}/plot/{sample}/{ctg}_plot_layout.toml`
        * `cenplot draw -t results/plot/${sm}/${ctg}_plot_layout.toml -c <(echo ${ctg}) -d results/plot/${sm}/`

## Parameters
|parameter|description|default|
|-|-|-|
|`species`|Species to use for RepeatMasker Dfam database. **CDR-Finder was designed for human centromeres so performance in other species is untested. Use with caution.** See https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html.|human|
|`restrict_alr`|Restrict CDR search to sequence annotated as `ALR/Alpha` by RepeatMasker. **Disabling this will cause false positives and should only be set to `false` if the input region is well curated.**|true|
|`window_size`|Size of the methylation windows to average over.|5000|
|`alr_threshold`|Size of ALR repeat stretches to include in search of CDR.|100,000|
|`bp_merge`|Distance in bases to merge adjacent CDRs. Can be omitted.|1|
|`bp_alr_merge`|Distance in bases to merge adjacent alpha-satellite regions.|1,000|
|`bp_edge`|Distance in bases to check CDR edges. Used to determine height of dip. Large values give better estimates of true height.|500,000|
|`height_perc_valley_threshold`|Threshold percent of the median methylation percentage needed as the minimal height of a valley from the median. Larger values filter for deeper valleys.|0.34|
|`prom_perc_valley_threshold`|Threshold percent of the median methylation percentage needed as the minimal [prominence](https://en.wikipedia.org/wiki/Topographic_prominence) of a valley from the median. Larger values filter for more prominent valleys. Helps in removing low-confidence CDRs.|0.3|
|`edge_height_heuristic`|Heuristic used when determining edge height of CDR. Either min, max, or avg.|min|
|`extend_edges_std`|Extend edges of CDR until the mean average methylation signal and some number of standard deviations is reached. May fix smaller, less prominent CDRs being missed. A value of 0 is the mean while +1/-1 is one stdev above/below the mean.|-1|
|`baseline_avg_methyl`|Baseline average methylation per region. Adjusts `height_perc_valley_threshold` and `prom_perc_valley_threshold` based on the average methyl percent relative to baseline. Will only increase thresholds and never reduce them. Reduces false positives in centromeres with low methylation relative to coverage.|0.4|

See [`docs/ISSUES.md`](docs/ISSUES.md) for edge-cases and effects of parameter tuning.

## Testing
Set up the conda environment and pull test data with [`git-lfs`](https://git-lfs.com/).
```bash
conda env create --name cdr_finder -f env.yaml
conda activate cdr_finder
git lfs install && git lfs pull
```

To run the test case on chr8 and chr21.
```bash
snakemake -c 1 -p --sdm conda --configfile test/config/config.yaml
```

To run integration tests.
```bash
pytest -vvv
```

To replicate the README image, run the test case. And then run the following:
```bash
# Replace title
sed -i 's/title = "{chrom}"/title = "CHM13 chromosome 8 centromere"/g' results/plot/CHM13/chr8_plot_layout.toml
# Get cenplot conda env
conda activate $(grep -l "cenplot" .snakemake/conda/*.yaml | xargs echo | sed 's/.yaml//g')
# Draw image again.
cenplot draw -t results/plot/CHM13/chr8_plot_layout.toml -c <(echo chr8) -d results/plot/CHM13/readme
```

## Containerization
Requires `root` user and `docker`.
```bash
make docker
```

```bash
make singularity
```

## Cite
* Francesco Kumara Mastrorosa, Keisuke K Oshima, Allison N Rozanski, William T Harvey, Evan E Eichler, Glennis A Logsdon, Identification and annotation of centromeric hypomethylated regions with CDR-Finder, Bioinformatics, Volume 40, Issue 12, December 2024, btae733, https://doi.org/10.1093/bioinformatics/btae733
