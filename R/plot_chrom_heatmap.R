#' @title Plot genomic feature density heatmap
#' @description Plot genomic feature density \bold{\emph{heatmap}}.
#' @author benben-miao
#'
#' @return A \bold{\emph{ggplot object}} showing genomic feature density heatmap.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param feature Genomic feature to quantify. (\bold{\emph{"gene"}}, "exon", "CDS", "promoter").
#' @param bin_size Window size (bp) for density calculation. (\bold{\emph{1e6}}).
#' @param orientation Coordinate orientation. (\bold{\emph{"horizontal"}}, "vertical").
#' @param palette Continuous color palette for density. (\bold{\emph{c("#ffffff", "#0055aa")}}).
#' @param alpha Tile alpha for density heatmap. (\bold{\emph{0.9}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example GFF3 file in GAnnoViz
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff3.gz",
#'   package = "GAnnoViz")
#'
#' # Gene density heatmap
#' plot_chrom_heatmap(
#'   gff_file = gff_file,
#'   format = "auto",
#'   feature = "gene",
#'   bin_size = 1e6,
#'   orientation = "horizontal",
#'   palette = c("#ffffff", "#0055aa"),
#'   alpha = 0.9
#' )
#'
plot_chrom_heatmap <- function(gff_file,
                                       format = "auto",
                                       feature = "gene",
                                       bin_size = 1e6,
                                       orientation = "horizontal",
                                       palette = c("#ffffff", "#0055aa"),
                                       alpha = 0.9) {
  # Parameters check
  if (!(feature %in% c("gene", "exon", "CDS", "promoter")))
    stop("feature must be one of: 'gene','exon','CDS','promoter'")

  # Genomic ranges
  gr <- if (feature == "gene") {
    extract_genes(gff_file = gff_file,
                  format = format,
                  gene_info = "all")
  } else if (feature == "exon") {
    extract_exons(gff_file = gff_file,
                  format = format,
                  exon_info = "all")
  } else if (feature == "CDS") {
    extract_cds(gff_file = gff_file,
                format = format,
                cds_info = "all")
  } else {
    extract_promoters(gff_file = gff_file,
                      format = format,
                      promoter_info = "all")
  }

  # Chromosomes
  chrom_raw <- as.character(GenomicRanges::seqnames(gr))
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  gr <- gr[keep_mask]

  # Data
  df <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = BiocGenerics::start(gr),
    end = BiocGenerics::end(gr),
    stringsAsFactors = FALSE
  )
  suppressWarnings({
    df$chrom_num <- as.numeric(stringr::str_extract(df$chrom, "\\d+"))
  })
  order_index <- order(is.na(df$chrom_num), df$chrom_num, df$chrom)
  chrom_levels <- unique(df$chrom[order_index])

  seqlens <- stats::aggregate(end ~ chrom, data = df, FUN = max)
  seqlens <- setNames(seqlens$end, seqlens$chrom)
  tiles <- GenomicRanges::tileGenome(
    seqlengths = seqlens,
    tilewidth = bin_size,
    cut.last.tile.in.chrom = TRUE
  )
  counts <- GenomicRanges::countOverlaps(tiles, gr)
  df_bin <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(tiles)),
    start = BiocGenerics::start(tiles),
    end = BiocGenerics::end(tiles),
    count = as.integer(counts),
    stringsAsFactors = FALSE
  )
  df_bin$center <- (df_bin$start + df_bin$end) / 2

  length_map <- stats::aggregate(end ~ chrom, data = df_bin, FUN = max)
  names(length_map) <- c("chrom", "length")
  y_max <- max(length_map$length, na.rm = TRUE)

  # Plot
  p <- ggplot2::ggplot(df_bin, ggplot2::aes(x = center, y = factor(chrom, levels = chrom_levels))) +
    ggplot2::geom_tile(ggplot2::aes(
      width = end - start,
      height = 0.5,
      fill = count
    ),
    alpha = alpha) +
    ggplot2::scale_fill_gradientn(colors = palette) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_number(),
                                breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::labs(x = "Genomic position",
                  y = "Chromosome",
                  fill = paste0(feature, " density")) +
    ggplot2::coord_cartesian(xlim = c(0, y_max * 1.05)) +
    my_theme()
  if (orientation == "horizontal")
    p <- p + ggplot2::coord_flip()
  p
}
