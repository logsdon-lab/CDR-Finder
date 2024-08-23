library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(stringr)
library(argparser)
library(data.table)
library(ggnewscale)


reformat_repeatmasker_bed <- function(df) {
    df_intermediate <- df
    # Get mask for Satellite repeats and rename repeat types
    mask <- (df_intermediate$rClass == "Satellite/centr") | (df_intermediate$rClass == "Satellite")
    df_intermediate$rClass[mask] <- df_intermediate$type[mask]
    df_intermediate$rClass <- sub("/ERVK", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/ERVL", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/ERV1", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/CR1", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/L1", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/L2", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/RTE-X", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/RTE-BovB", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/Gypsy", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("-MaLR", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/Alu", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/Deu", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/MIR", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("?", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/hAT", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/hAT-Blackjack", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/hAT-Charlie", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/MULE-MuDR", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/PiggyBac", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/TcMar-Mariner", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/TcMar", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/TcMar?", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/hAT-Tip100", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/TcMar-Tigger", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/Dong-R4", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("/tRNA", "", df_intermediate$rClass)
    df_intermediate$rClass <- sub("DNA-Tc2", "DNA", df_intermediate$rClass)
    df_intermediate$rClass <- sub("DNA?", "DNA", df_intermediate$rClass)
    df_intermediate$rClass <- sub("DNA-Blackjack", "DNA", df_intermediate$rClass)
    df_intermediate$rClass <- sub("DNA-Charlie", "DNA", df_intermediate$rClass)
    df_intermediate$rClass <- sub("DNA-Tigger", "DNA", df_intermediate$rClass)
    df_intermediate$rClass <- sub("DNA-Tip100", "DNA", df_intermediate$rClass)
    df_intermediate$rClass <- sub("GSATX", "GSAT", df_intermediate$rClass)
    df_intermediate$rClass <- sub("LTR\\S", "LTR", df_intermediate$rClass)
    df_intermediate$type <- as.factor(df_intermediate$type)

    # Rename classes.
    df_intermediate$rClass <- sub("SAR", "HSat1A", df_intermediate$rClass)
    df_intermediate$rClass <- sub("HSAT", "HSat1B", df_intermediate$rClass)
    df_intermediate$rClass <- sub("HSATII", "HSat2", df_intermediate$rClass)
    df_intermediate$rClass <- sub("(CATTC)n", "HSat2", df_intermediate$rClass)
    df_intermediate$rClass <- sub("(GAATG)n", "HSat2", df_intermediate$rClass)

    # Adjust for reverse complement.
    df_intermediate$start2 <- df_intermediate$start
    df_intermediate$end2 <- df_intermediate$end
    mask <- df_intermediate$C == "C"
    df_intermediate$start2[mask] <- df_intermediate$end[mask]
    df_intermediate$end2[mask] <- df_intermediate$start[mask]

    return(df_intermediate)
}

read_repeatmasker_bed <- function(input_file) {
    # read in BED file
    df <- fread(
        input_file,
        select = c(1:6),
        stringsAsFactors = TRUE,
        fill = TRUE,
        sep = "\t",
        quote = "",
        header = FALSE,
        col.names = c("chr", "start", "end", "rClass", "type", "strand")
    )

    # start2 and end2
    df_reformatted <- reformat_repeatmasker_bed(df)
    return(df_reformatted)
}

get_colors <- function() {
    default_colors <- c(
        "#522758",
        "#84bac6", "#84bac6",
        "#f2b80b", "#ad8c2a"
    )
    names(default_colors) <- levels(as.factor(c("ALR/Alpha", "HSat1A", "HSat1B", "HSat2", "HSat3")))
    return(default_colors)
}

p <- arg_parser("Plot cumulative centromere HOR array lengths.")
p <- add_argument(
  p, "--input_methyl",
  help = "Input binned average methylation frequency bed.", type = "character"
)
p <- add_argument(
  p, "--input_cdr",
  help = "Input CDR regions.", type = "character"
)
p <- add_argument(
  p, "--input_rm",
  help = "Input RepeatMasker annotations.", type = "character"
)
p <- add_argument(
  p, "--input_target_regions",
  help = "Input target regions.", type = "character"
)
p <- add_argument(
  p, "--output",
  help = "Output binned average methylation frequency plot PNG.", type = "character",
  default = "plot.png"
)
argv <- parse_args(p)

df_target_regions <- fread(
    argv$input_target_regions,
    header = FALSE,
    col.names = c("chr", "start", "end")
)
df_cdr <- fread(
    argv$input_cdr,
    header = FALSE,
    col.names = c("chr", "start", "end", "desc")
)
df_cdr <- df_cdr %>% filter(desc == "high_confidence")

# Load binned freq file
df_methyl_binned <- fread(
    argv$input_methyl,
    header = FALSE,
    col.names = c("chr", "start", "end", "meth_prob", "index")
)
# Adjust coordinates
df_methyl_binned <- df_methyl_binned %>%
    left_join(df_target_regions, join_by(chr)) %>%
    mutate(
        start = start.x + start.y,
        end = end.x + start.y
    )

# Adjust coordinates
df_rm_out <- read_repeatmasker_bed(argv$input_rm) %>%
    left_join(df_target_regions, join_by(chr)) %>%
    mutate(
        start = start2 + start.y,
        end = end2 + start.y,
    )

# Plot bins
plt <- ggplot() +
    geom_area(
        data = df_methyl_binned,
        aes(x = start, y = as.numeric(meth_prob)),
        fill = "#008080",
        colour = "black",
        alpha = 0.5
    ) +
    geom_rect(
        data = df_cdr,
        aes(
            xmin = start,
            ymin = 0,
            xmax = end,
            ymax = 100,
        ),
        fill = "red",
        alpha = 0.5
    ) +
    new_scale_color() +
    geom_segment(
        data = df_rm_out,
        aes(
            x = start,
            y = 110,
            xend = end,
            yend = 110,
            colour = rClass
        ),
        linewidth = 10
    ) +
    scale_color_manual(values = get_colors()) +
    facet_wrap(vars(chr), ncol = 1, scales = "free_x") +
    ylim(0, 110) +
    scale_x_continuous(labels = unit_format(scale = 1e-6, accuracy=0.1, unit="")) +
    xlab("Position (Mbp)") +
    ylab("Average Methylation Percent") +
    labs(colour = "Repeats") +
    theme_classic() +
    theme(
        # Place labels outside and remove bg.
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black'),
        strip.placement = "outside"
    )

ggsave(
    argv$output,
    plt,
    width = 12,
    height = nrow(df_target_regions) * 2.5,
    limitsize = FALSE
)
