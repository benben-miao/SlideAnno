#' @title Plot DEGs up/down along chromosomes
#' @description Plot \bold{\emph{DEGs up/down}} along chromosomes.
#' @author benben-miao
#'
#' @return A \bold{\emph{ggplot object}} visualizing DEGs along chromosomes.
#' @param deg_file DEG table from \bold{\emph{DESeq2}} analysis.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param id_col Gene IDs column name. (\bold{\emph{"GeneID"}}).
#' @param fc_col Log2(fold change) column name. (\bold{\emph{"log2FoldChange"}}).
#' @param orientation Coordinate orientation. (\bold{\emph{"horizontal"}}, "vertical").
#' @param chrom_alpha Chromosome bar alpha. (\bold{\emph{0.1}}).
#' @param chrom_color Chromosome bar color. (\bold{\emph{"#008888"}}).
#' @param bar_height Chromosome bar thickness in y units. (\bold{\emph{0.8}}).
#' @param point_size Point size. (\bold{\emph{2}}).
#' @param point_alpha Point alpha. (\bold{\emph{0.3}}).
#' @param up_color Color for up-regulated genes. (\bold{\emph{"#ff0000"}}).
#' @param down_color Color for down-regulated genes. (\bold{\emph{"#008800"}}).
#' @param mark_style Marker style for DEGs. (\bold{\emph{"point"}}, "line").
#' @param line_width Line width. (\bold{\emph{0.6}}).
#' @param line_height Line height relative to bar radius. (\bold{\emph{0.8}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # DEG results from DESeq2
#' deg_file <- system.file(
#'     "extdata",
#'     "example.deg",
#'     package = "GAnnoViz")
#'
#' # Genomic structure annotation
#' gff_file <- system.file(
#'     "extdata",
#'     "example.gff3.gz",
#'     package = "GAnnoViz")
#'
#' # Plot DEGs along chromosomes
#' plot_deg_exp(
#'   deg_file = deg_file,
#'   gff_file = gff_file,
#'   format = "auto",
#'   id_col = "GeneID",
#'   fc_col = "log2FoldChange",
#'   orientation = "horizontal",
#'   chrom_alpha = 0.1,
#'   chrom_color = "#008888",
#'   bar_height = 0.8,
#'   point_size = 2,
#'   point_alpha = 0.3,
#'   up_color = "#ff0000",
#'   down_color = "#008800",
#'   mark_style = "point",
#'   line_width = 0.6,
#'   line_height = 0.8)
#'
plot_deg_exp <- function(deg_file,
                         gff_file,
                         format = "auto",
                         id_col = "GeneID",
                         fc_col = "log2FoldChange",
                         orientation = "horizontal",
                         chrom_alpha = 0.1,
                         chrom_color = "#008888",
                         bar_height = 0.8,
                         point_size = 2,
                         point_alpha = 0.3,
                         up_color = "#ff0000",
                         down_color = "#008800",
                         mark_style = "point",
                         line_width = 0.6,
                         line_height = 0.8) {

  # DEG annotation
  mark_style <- match.arg(mark_style)
  df <- anno_deg_chrom(
    deg_file = deg_file,
    gff_file = gff_file,
    format = format,
    id_col = id_col,
    fc_col = fc_col,
    use_strand = FALSE,
    drop_unmapped = TRUE
  )
  df <- df %>%
    dplyr::rename(chr = chrom,
                  log2fc = score,
                  geneid = gene_id)
  chrom_raw <- as.character(df$chr)
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  df <- df[keep_mask, , drop = FALSE]
  df <- df[stats::complete.cases(df$chr, df$start, df$end, df$log2fc), , drop = FALSE]
  if (nrow(df) == 0)
    stop("No DEGs with valid positions found")

  # Chromosomes
  df$pos <- (df$start + df$end) / 2
  chrom_chr <- as.character(df$chr)
  chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
  suppressWarnings({
    chrom_num <- as.numeric(chrom_num)
  })
  order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
  chrom_levels <- unique(chrom_chr[order_index])
  df$chrom <- factor(df$chr, levels = chrom_levels)

  length_map <- stats::aggregate(end ~ chr, data = df, FUN = max)
  length_map <- length_map[match(chrom_levels, length_map$chr), , drop = FALSE]
  length_map$y <- seq_along(chrom_levels)
  length_map$ymin <- length_map$y - bar_height / 2
  length_map$ymax <- length_map$y + bar_height / 2
  length_map$radius <- bar_height / 2

  df$y <- length_map$y[match(df$chr, length_map$chr)]
  df$state <- ifelse(df$log2fc >= 0, "up", "down")
  col_map <- c(up = up_color, down = down_color)

  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = length_map,
      ggplot2::aes(
        xmin = 0,
        xmax = end,
        ymin = ymin,
        ymax = ymax
      ),
      fill = chrom_color,
      color = NA,
      alpha = chrom_alpha
    ) +
    {
      if (identical(mark_style, "point")) {
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = pos, y = y, color = state),
          size = point_size,
          alpha = point_alpha
        )
      } else {
        # line style: vertical ticks with rounded ends
        df_seg <- df
        r <- length_map$radius[match(df_seg$chr, length_map$chr)]
        h <- pmax(0.1, as.numeric(r) * line_height)
        df_seg$y1 <- df_seg$y - h
        df_seg$y2 <- df_seg$y + h
        ggplot2::geom_segment(
          data = df_seg,
          ggplot2::aes(
            x = pos,
            xend = pos,
            y = y1,
            yend = y2,
            color = state
          ),
          linewidth = line_width,
          lineend = "round",
          alpha = point_alpha
        )
      }
    } +
    ggplot2::scale_color_manual(
      values = col_map,
      drop = FALSE,
      name = NULL,
      labels = c(down = "Down", up = "Up")
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::scale_y_continuous(
      breaks = length_map$y,
      labels = chrom_levels,
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(x = "Genomic position", y = "Chromosome") +
    my_theme()
  if (orientation == "horizontal")
    p <- p + ggplot2::coord_flip()
  p
}
