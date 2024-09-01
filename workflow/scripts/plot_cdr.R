library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(stringr)
library(argparser)
library(data.table)
library(ggnewscale)
library(patchwork)


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
  p, "--output_dir",
  help = "Output directory of binned average methylation frequency plot PNG.", type = "character",
  default = "plot"
)
argv <- parse_args(p)

df_cdr <- fread(
    argv$input_cdr,
    header = FALSE,
    select = c(1:3),
    col.names = c("chr", "start", "end")
)

# Load binned freq file
df_methyl_binned <- fread(
    argv$input_methyl,
    header = FALSE,
    select = c(1:5),
    col.names = c("chr", "start", "end", "meth_prob", "cov")
)

df_rm_out <- read_repeatmasker_bed(argv$input_rm)

# Make directory
dir.create(argv$output_dir, showWarnings = FALSE)

# 1st coverage, 2nd chr, 3rd average methylation freq
# Plot bins
for (chr_name in unique(df_methyl_binned$chr)) {
    plt_cov <- ggplot(
        data = df_methyl_binned %>% filter(chr == chr_name) %>% select(start, cov),
        aes(x = start, y = cov),
    ) +
    geom_area(position = "identity") +
    ylab("Read Depth") +
    theme_classic() +
    theme(
        axis.title.x =  element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
    )

    plt_methyl <- ggplot() +
        geom_segment(
            data = df_cdr %>% filter(chr == chr_name),
            aes(x = start, y = 130, xend = end, yend = 130),
        ) +
        geom_segment(
            data = df_rm_out %>% filter(chr == chr_name),
            aes(
                x = start,
                y = 115,
                xend = end,
                yend = 115,
                colour = rClass
            ),
            linewidth = 10
        ) +
        scale_color_manual(values = get_colors()) +
        labs(colour = "Repeats") +
        geom_area(
            data = df_methyl_binned %>% filter(chr == chr_name),
            aes(x = start, y = as.numeric(meth_prob)),
            fill = "black",
        ) +
        # Annotate CDR with with red overlap rectangle and horizontal bar.
        geom_rect(
            data = df_cdr %>% filter(chr == chr_name),
            aes(
                xmin = start,
                ymin = 0,
                xmax = end,
                ymax = 100,
            ),
            fill = "red",
            alpha = 0.5
        ) +
        scale_y_continuous(
            labels = unit_format(unit="%"),
            breaks = seq(0, 100, by = 20)
        ) +
        ylab("Average Methylation Percent") +
        theme_classic() +
        theme(
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
        )

    plt_ht_prop <- c(1.5, 3)
    plt_final <- plt_cov / plt_methyl +
        plot_layout(
            ncol = 1,
            heights = unit(plt_ht_prop, c('null', 'null'))
        ) +
        plot_annotation(
            title = chr_name,
            theme = theme(plot.title = element_text(size = 18))
        ) &
        scale_x_continuous(labels = unit_format(scale = 1e-6, accuracy=0.1, unit="")) &
        xlab("Position (Mbp)")

    outfile_pdf <- file.path(argv$output, paste0(chr_name, ".pdf"))
    outfile_png <- file.path(argv$output, paste0(chr_name, ".png"))
    ggsave(
        outfile_pdf,
        device = "pdf",
        plt_final,
        width = 12,
        height = sum(plt_ht_prop),
    )
    ggsave(
        outfile_png,
        device = "png",
        plt_final,
        width = 12,
        height = sum(plt_ht_prop),
    )
}
