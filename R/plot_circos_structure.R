#' @title Circos-style chromosome gene density ring
#' @description Plot \bold{\emph{circos-style ring}} of windowed \bold{\emph{gene density}} across all chromosomes, with optional gene highlights.
#'
#' @return A \bold{\emph{ggplot object}} using polar coordinates.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param bin_size Window size (bp) for density calculation. (\bold{\emph{1e6}}).
#' @param palette Continuous color palette for density. (\bold{\emph{c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5")}}).
#' @param ring_height Baseline ring height (arbitrary units). (\bold{\emph{1}}).
#' @param highlight_genes Optional character vector of gene identifiers to highlight. (\bold{\emph{NULL}}).
#' @param annotate Annotation mode for highlights. (\bold{\emph{"id"}}, "name").
#' @param highlight_color Color for highlighted genes. (\bold{\emph{"#ff0000"}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example GFF3 file in GAnnoViz
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff",
#'   package = "GAnnoViz")
#'
#' # Circos ring with optional highlights
#' plot_circos_structure(
#'   gff_file = gff_file,
#'   format = "auto",
#'   bin_size = 1e6,
#'   palette = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5"),
#'   ring_height = 1,
#'   highlight_genes = c("HdF029609", "HdF029610"),
#'   annotate = "id",
#'   highlight_color = "#ff0000"
#' )
plot_circos_structure <- function(gff_file,
                                  format = "auto",
                                  bin_size = 1e6,
                                  palette = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5"),
                                  ring_height = 1,
                                  highlight_genes = NULL,
                                  annotate = "id",
                                  highlight_color = "#ff0000") {
  genes <- extract_genes(gff_file = gff_file, format = format, gene_info = "all")
  chrom_raw <- as.character(GenomicRanges::seqnames(genes))
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix) stringr::str_detect(chrom_raw, "^(?i)chr") else rep(TRUE, length(chrom_raw))
  genes <- genes[keep_mask]

  df <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(genes)),
    start = BiocGenerics::start(genes),
    end = BiocGenerics::end(genes),
    gene_id = genes$gene_id,
    stringsAsFactors = FALSE
  )
  suppressWarnings({
    df$chrom_num <- as.numeric(stringr::str_extract(df$chrom, "\\d+"))
  })
  order_index <- order(is.na(df$chrom_num), df$chrom_num, df$chrom)
  chrom_levels <- unique(df$chrom[order_index])

  seqlens <- stats::aggregate(end ~ chrom, data = df, FUN = max)
  seqlens <- setNames(seqlens$end, seqlens$chrom)
  tiles <- GenomicRanges::tileGenome(seqlengths = seqlens, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
  counts <- GenomicRanges::countOverlaps(tiles, genes)
  df_bin <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(tiles)),
    start = BiocGenerics::start(tiles),
    end = BiocGenerics::end(tiles),
    count = as.integer(counts),
    stringsAsFactors = FALSE
  )

  df_bin <- df_bin %>% dplyr::group_by(chrom) %>% dplyr::arrange(start, .by_group = TRUE)
  df_bin$bin_index <- dplyr::row_number()
  span_map <- df_bin %>% dplyr::group_by(chrom) %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::mutate(start_idx = cumsum(dplyr::lag(n, default = 0)) + 1)
  df_bin <- dplyr::left_join(df_bin, span_map, by = "chrom")
  df_bin$theta <- df_bin$start_idx + df_bin$bin_index - 1

  chr_starts <- span_map$start_idx
  chr_starts_df <- data.frame(x = chr_starts)

  p <- ggplot2::ggplot(df_bin, ggplot2::aes(x = theta, y = ring_height)) +
    ggplot2::geom_tile(ggplot2::aes(fill = count), width = 1, height = 0.2) +
    ggplot2::geom_vline(data = chr_starts_df, ggplot2::aes(xintercept = x), linewidth = 0.2, color = "#dddddd") +
    ggplot2::scale_fill_gradientn(colors = palette) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Gene density") +
    ggplot2::coord_polar() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")

  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    df_h <- df
    if (annotate == "id") {
      df_h <- df_h[df_h$gene_id %in% highlight_genes, , drop = FALSE]
    } else {
      if (!is.null(genes$gene_name)) {
        df_h$gene_name <- genes$gene_name
        df_h <- df_h[df_h$gene_name %in% highlight_genes, , drop = FALSE]
      } else {
        df_h <- df_h[0, ]
      }
    }
    if (nrow(df_h) > 0) {
      pos_map <- df_bin %>% dplyr::group_by(chrom) %>% dplyr::summarise(xmin = min(theta), xmax = max(theta))
      df_h <- dplyr::left_join(df_h, pos_map, by = "chrom")
      df_h$theta <- df_h$xmin + floor((df_h$start / (df_h$end)) * (df_h$xmax - df_h$xmin))
      p <- p + ggplot2::geom_point(data = df_h, ggplot2::aes(x = theta, y = ring_height + 0.12), color = highlight_color, size = 1.5)
    }
  }
  p
}

