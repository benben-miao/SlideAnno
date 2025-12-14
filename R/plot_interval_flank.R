#' @title Gene neighborhood architecture around a focal gene
#' @description Plot \bold{\emph{gene neighborhood}} around a focal \bold{\emph{gene}} with distances, arrows for strand, and optional promoter track.
#' @author benben-miao
#'
#' @return A \bold{\emph{ggplot object}} of genes within a flank window around the focal gene.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param gene_id Focal gene id consistent with GFF/GTF. (\bold{\emph{necessary}}).
#' @param flank_upstream Upstream flank window (bp). (\bold{\emph{200000}}).
#' @param flank_downstream Downstream flank window (bp). (\bold{\emph{200000}}).
#' @param show_promoters Whether to draw a promoter track. (\bold{\emph{TRUE}}).
#' @param upstream Promoter upstream (bp). (\bold{\emph{2000}}).
#' @param downstream Promoter downstream (bp). (\bold{\emph{200}}).
#' @param arrow_length Length of strand arrows (pt). (\bold{\emph{5}}).
#' @param arrow_unit Arrow length unit. (\bold{\emph{"pt"}}, "mm").
#' @param gene_color Gene bar color. (\bold{\emph{"#0088ff"}}).
#' @param promoter_color Promoter bar color. (\bold{\emph{"#ff8800"}}).
#' @param label_size Label text size. (\bold{\emph{3}}).
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
#' # Neighborhood around a focal gene on its chromosome
#' plot_interval_flank(
#'   gff_file = gff_file,
#'   format = "auto",
#'   gene_id = "HdF029609",
#'   flank_upstream = 200000,
#'   flank_downstream = 200000,
#'   show_promoters = TRUE,
#'   upstream = 2000,
#'   downstream = 200,
#'   arrow_length = 5,
#'   arrow_unit = "pt",
#'   gene_color = "#0088ff",
#'   promoter_color = "#ff8800",
#'   label_size = 3
#' )
#'
plot_interval_flank <- function(gff_file,
                                   format = "auto",
                                   gene_id,
                                   flank_upstream = 200000,
                                   flank_downstream = 200000,
                                   show_promoters = TRUE,
                                   upstream = 2000,
                                   downstream = 200,
                                   arrow_length = 5,
                                   arrow_unit = "pt",
                                   gene_color = "#0088ff",
                                   promoter_color = "#ff8800",
                                   label_size = 3) {
  genes <- extract_genes(gff_file = gff_file, format = format, gene_info = "all")
  if (missing(gene_id)) stop("gene_id is required")
  if (!gene_id %in% genes$gene_id) stop("Input gene_id not found in GFF/GTF: ", gene_id)

  focal <- genes[genes$gene_id == gene_id]
  chrom_id <- as.character(GenomicRanges::seqnames(focal))
  g_start <- BiocGenerics::start(focal)
  g_end <- BiocGenerics::end(focal)
  region_start <- max(1, g_start - flank_upstream)
  region_end <- g_end + flank_downstream
  region <- GenomicRanges::GRanges(seqnames = chrom_id, ranges = IRanges::IRanges(start = region_start, end = region_end))

  chrom_raw <- as.character(GenomicRanges::seqnames(genes))
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix) stringr::str_detect(chrom_raw, "^(?i)chr") else rep(TRUE, length(chrom_raw))
  genes <- genes[keep_mask]
  hits_g <- GenomicRanges::findOverlaps(genes, region, type = "any", ignore.strand = TRUE)
  genes <- genes[S4Vectors::queryHits(hits_g)]
  if (length(genes) == 0) stop("No neighboring genes in the specified flank window")

  df_genes <- data.frame(
    gene_id = genes$gene_id,
    start = BiocGenerics::start(genes),
    end = BiocGenerics::end(genes),
    strand = as.character(GenomicRanges::strand(genes)),
    stringsAsFactors = FALSE
  )
  df_genes$y <- seq_len(nrow(df_genes))

  df_gene_rect <- data.frame(
    xmin = df_genes$start,
    xmax = df_genes$end,
    ymin = df_genes$y - 0.35,
    ymax = df_genes$y + 0.35,
    gene_id = df_genes$gene_id,
    stringsAsFactors = FALSE
  )

  # Strand arrows (towards transcription direction)
  df_arrows <- df_genes
  df_arrows$x <- ifelse(df_arrows$strand == "+", df_arrows$end, df_arrows$start)
  df_arrows$xend <- df_arrows$x + ifelse(df_arrows$strand == "+", 1, -1) * (region_end - region_start) * 0.01
  df_arrows$y <- df_arrows$y
  df_arrows$yend <- df_arrows$y

  # Promoter track (transcript-level)
  df_prom <- data.frame()
  if (isTRUE(show_promoters)) {
    promoters <- extract_promoters(gff_file = gff_file, format = format, upstream = upstream, downstream = downstream, promoter_info = "all")
    hits_p <- GenomicRanges::findOverlaps(promoters, region, type = "any", ignore.strand = TRUE)
    promoters <- promoters[S4Vectors::queryHits(hits_p)]
    df_prom <- data.frame(
      start = BiocGenerics::start(promoters),
      end = BiocGenerics::end(promoters),
      y = 0,
      stringsAsFactors = FALSE
    )
  }

  x_breaks <- 10
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = df_gene_rect,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = gene_color,
      color = NA,
      alpha = 0.9
    ) +
    ggplot2::geom_segment(
      data = df_arrows,
      ggplot2::aes(x = x, xend = xend, y = y, yend = y),
      color = "#333333",
      linewidth = 0.6,
      arrow = grid::arrow(length = grid::unit(arrow_length, arrow_unit), type = "closed")
    ) +
    ggplot2::geom_text(
      data = df_genes,
      ggplot2::aes(x = pmax(df_genes$end, df_genes$start) + (region_end - region_start) * 0.02, y = y, label = gene_id),
      size = label_size,
      color = "#000000",
      hjust = 0,
      vjust = 0.5
    )

  if (isTRUE(show_promoters) && nrow(df_prom) > 0) {
    p <- p + ggplot2::geom_rect(
      data = df_prom,
      ggplot2::aes(xmin = start, xmax = end, ymin = -0.3, ymax = 0.3),
      fill = promoter_color,
      color = NA,
      alpha = 0.8
    )
  }

  p <- p +
    ggplot2::scale_y_continuous(
      breaks = c(0, df_genes$y),
      labels = c("Promoter", df_genes$gene_id),
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(
      x = "Genomic position",
      y = "Gene / Promoter",
      title = paste0("Chrom: ", chrom_id, "; Focal: ", gene_id)
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = x_breaks)) +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    my_theme()
  p
}
