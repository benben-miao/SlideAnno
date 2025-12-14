#' @title Plot genomic FST with Top-N gene annotations
#' @description Plot genomic FST with \bold{\emph{Top-N gene annotations}}.
#' @author benben-miao
#'
#' @return A \bold{\emph{ggplot object}} representing genomic FST with Top-N gene annotations.
#' @param fst_file FST sliding window results. (\bold{\emph{CHROM, BIN_START, BIN_END, WEIGHTED_FST, N_VARIANTS}}.
#' @param gff_file Path to \bold{\file{GFF3/GTF}} file as input.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param chrom_id Chromosome identifier (\bold{\emph{"chr1"}}).
#' @param top_n Number of top genes to annotate. (\bold{\emph{20}}).
#' @param orientation Coordinate orientation. (\bold{\emph{"vertical"}}, "horizontal").
#' @param smooth_span Span for local regression smoothing. (\bold{\emph{0.1}}).
#' @param fst_color Point color for FST. (\bold{\emph{"#0088ff"}}).
#' @param point_size Point size. (\bold{\emph{1}}).
#' @param point_alpha Point alpha. (\bold{\emph{0.3}}).
#' @param label_size Text size for gene labels. (\bold{\emph{3}}).
#' @param connector_dx1 First connector horizontal offset (bp). (\bold{\emph{2e4}}).
#' @param connector_dx2 Second connector horizontal offset (bp). (\bold{\emph{4e4}}).
#' @param gap_frac Minimum vertical gap between labels (fraction of FST range). (\bold{\emph{0.05}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example data in GAnnoViz
#' fst_file <- system.file(
#'   "extdata",
#'   "example.fst",
#'   package = "GAnnoViz")
#'
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff3.gz",
#'   package = "GAnnoViz")
#'
#' # Chromosome FST with Top-20 gene annotations on chr11
#' plot_snp_anno(
#'   fst_file = fst_file,
#'   gff_file = gff_file,
#'   format = "auto",
#'   chrom_id = "chr2",
#'   top_n = 20,
#'   orientation = "vertical",
#'   smooth_span = 0.5,
#'   fst_color = "#0088ff",
#'   point_size = 1,
#'   point_alpha = 0.3,
#'   label_size = 3,
#'   connector_dx1 = 2e4,
#'   connector_dx2 = 4e4,
#'   gap_frac = 0.05
#' )
#'
plot_snp_anno <- function(fst_file,
                          gff_file,
                          format = "auto",
                          chrom_id,
                          top_n = 20,
                          orientation = "vertical",
                          smooth_span = 0.5,
                          fst_color = "#0088ff",
                          point_size = 1,
                          point_alpha = 0.3,
                          label_size = 3,
                          connector_dx1 = 2e4,
                          connector_dx2 = 4e4,
                          gap_frac = 0.05) {
  # FST results
  fst <- utils::read.table(
    fst_file,
    header = TRUE,
    sep = "\t",
    na.strings = "NA",
    stringsAsFactors = FALSE
  )
  fst <- fst %>% dplyr::filter(CHROM == chrom_id)
  if (nrow(fst) == 0)
    stop("No FST windows on selected chromosome")

  # Genes
  genes <- extract_genes(gff_file = gff_file,
                         format = format,
                         gene_info = "all")
  chrom_raw <- as.character(GenomicRanges::seqnames(genes))
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  genes <- genes[keep_mask]
  genes <- genes[as.character(GenomicRanges::seqnames(genes)) == chrom_id]

  df_genes <- data.frame(
    gene_id = genes$gene_id,
    start = BiocGenerics::start(genes),
    end = BiocGenerics::end(genes),
    stringsAsFactors = FALSE
  )
  df_genes$mid <- (df_genes$start + df_genes$end) / 2

  agg_fst_gene <- function(s, e, windows) {
    mask <- (windows$BIN_START <= e) & (windows$BIN_END >= s)
    if (!any(mask))
      return(NA_real_)
    w <- windows[mask, , drop = FALSE]
    ov <- pmin(e, w$BIN_END) - pmax(s, w$BIN_START) + 1
    ov[ov < 0] <- 0
    wsum <- sum(ov, na.rm = TRUE)
    if (wsum <= 0)
      return(NA_real_)
    sum(w$WEIGHTED_FST * ov, na.rm = TRUE) / wsum
  }

  df_genes$fst_wmean <- vapply(seq_len(nrow(df_genes)), function(i) {
    agg_fst_gene(df_genes$start[i], df_genes$end[i], fst)
  }, numeric(1))

  df_top <- df_genes[stats::complete.cases(df_genes$fst_wmean), , drop = FALSE]
  df_top <- df_top[order(-df_top$fst_wmean), ]
  df_top <- utils::head(df_top, n = top_n)

  p_base <- ggplot2::ggplot(fst, ggplot2::aes(x = BIN_START, y = WEIGHTED_FST)) +
    # ggplot2::geom_line(color = fst_color, linewidth = 0.6) +
    ggplot2::geom_point(color = fst_color,
                        size = point_size,
                        alpha = point_alpha) +
    ggplot2::geom_smooth(span = smooth_span, se = FALSE, color = fst_color) +
    ggplot2::labs(x = "Genomic position", y = "WEIGHTED FST") +
    my_theme()

  if (nrow(df_top) == 0) {
    if (orientation == "horizontal")
      p_base <- p_base + ggplot2::coord_flip()
    return(p_base)
  }

  y_min <- min(fst$WEIGHTED_FST, na.rm = TRUE)
  y_max <- max(fst$WEIGHTED_FST, na.rm = TRUE)
  y_range <- y_max - y_min
  base_y <- y_max + y_range * 0.08
  min_gap <- y_range * gap_frac

  df_top <- df_top[order(df_top$mid), , drop = FALSE]
  df_top$y_label <- NA_real_
  cur <- base_y
  for (i in seq_len(nrow(df_top))) {
    df_top$y_label[i] <- cur
    cur <- cur + min_gap
  }

  df_top$x1 <- df_top$mid
  df_top$y1 <- y_max
  df_top$x2 <- df_top$x1 + connector_dx1
  df_top$y2 <- df_top$y_label - min_gap * 0.3
  df_top$x3 <- df_top$x2 + connector_dx2
  df_top$y3 <- df_top$y_label

  df_seg1 <- data.frame(
    x = df_top$x1,
    y = df_top$y1,
    xend = df_top$x2,
    yend = df_top$y2
  )
  df_seg2 <- data.frame(
    x = df_top$x2,
    y = df_top$y2,
    xend = df_top$x3,
    yend = df_top$y3
  )

  p <- p_base +
    ggplot2::geom_segment(
      data = df_seg1,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = "#333333",
      linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = df_seg2,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = "#333333",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = df_top,
      ggplot2::aes(x = x3, y = y3, label = gene_id),
      size = label_size,
      color = "#000000",
      hjust = 0,
      vjust = 0.5
    )

  if (orientation == "horizontal")
    p <- p + ggplot2::coord_flip()
  p
}
