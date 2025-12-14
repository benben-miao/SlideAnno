#' @title Plot gene structures for a genomic interval
#' @description Plot gene structures for a \bold{\emph{genomic interval}}.
#'
#' @return Plot of genomic interval structure.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param chrom_id Chromosome id (must match the annotation).
#' @param start Window start coordinate (bp, genomic absolute coordinate).
#' @param end Window end coordinate (bp, genomic absolute coordinate).
#' @param x_breaks X axis breaks number. (\bold{\emph{10}}).
#' @param upstream Promoter upstream (bp). (\bold{\emph{2000}}).
#' @param downstream Promoter downstream (bp). (\bold{\emph{200}}).
#' @param feature_alpha Elements alpha. (\bold{\emph{0.8}}).
#' @param intron_width Intron line width. (\bold{\emph{1}}).
#' @param arrow_count Intron arrow number bold. (\bold{\emph{1}}).
#' @param arrow_length Intron arrows length（pt). (\bold{\emph{1}}).
#' @param arrow_unit Intron arrow length unit. (\bold{\emph{"pt"}}, "mm").
#' @param promoter_color Promoter color. (\bold{\emph{"#ff8800"}}).
#' @param utr5_color 5'UTR color. (\bold{\emph{"#008833"}}).
#' @param utr3_color 3'UTR color. (\bold{\emph{"#ff0033"}}).
#' @param exon_color Exon color. (\bold{\emph{"#0033ff"}}).
#' @param intron_color Intron color. (\bold{\emph{"#333333"}}).
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
#' plot_interval_structure(
#'   gff_file = gff_file,
#'   format = "auto",
#'   chrom_id = "chr1",
#'   start = 950000,
#'   end = 1180000,
#'   x_breaks = 10,
#'   upstream = 2000,
#'   downstream = 200,
#'   feature_alpha = 0.8,
#'   intron_width = 1,
#'   arrow_count = 1,
#'   arrow_length = 5,
#'   arrow_unit = "pt",
#'   promoter_color = "#ff8800",
#'   utr5_color = "#008833",
#'   utr3_color = "#ff0033",
#'   exon_color = "#0033ff",
#'   intron_color = "#333333"
#' )
#'
plot_interval_structure <- function(gff_file,
                              format = "auto",
                              chrom_id,
                              start,
                              end,
                              x_breaks = 10,
                              upstream = 2000,
                              downstream = 200,
                              feature_alpha = 0.8,
                              intron_width = 1,
                              arrow_count = 1,
                              arrow_length = 5,
                              arrow_unit = "pt",
                              promoter_color = "#ff8800",
                              utr5_color = "#008833",
                              utr3_color = "#ff0033",
                              exon_color = "#0033ff",
                              intron_color = "#333333") {
  # Check parameters
  if (missing(chrom_id) ||
      missing(start) ||
      missing(end))
    stop("Please provide chrom_id, start and end")
  if (!is.numeric(start) ||
      !is.numeric(end) ||
      length(start) != 1 ||
      length(end) != 1)
    stop("start/end must be numeric bp coordinates")
  if (end < start)
    stop("end must be >= start")

  # TXDB
  fmt <- resolve_gff_format(gff_file, format)
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(file = gff_file, format = fmt))

  # Genes in interval
  genes_all <- suppressWarnings(GenomicFeatures::genes(txdb))
  if (!(chrom_id %in% levels(genes_all@seqnames@values)))
    stop("chrom_id not found in annotation")
  chrom_vec <- as.character(GenomicRanges::seqnames(genes_all))
  keep_idx <- chrom_vec == chrom_id &
    BiocGenerics::end(genes_all) >= start &
    BiocGenerics::start(genes_all) <= end
  genes_sel <- genes_all[keep_idx]
  if (length(genes_sel) == 0)
    stop("No genes found in the specified interval")

  # Promoters, Exons, Introns
  tx_by_gene <- suppressWarnings(GenomicFeatures::transcriptsBy(txdb, by = "gene"))
  exons_by_gene <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "gene"))
  utr5_by_tx <- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE))
  utr3_by_tx <- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE))
  to_df <- function(gr, feature, gene_label, strand_label) {
    if (length(gr) == 0)
      return(data.frame())
    data.frame(
      gene_id = gene_label,
      chrom = as.character(GenomicRanges::seqnames(gr)),
      start = BiocGenerics::start(gr),
      end = BiocGenerics::end(gr),
      strand = strand_label,
      feature = feature,
      stringsAsFactors = FALSE
    )
  }
  df_all <- data.frame()
  gene_ids <- genes_sel$gene_id
  gene_starts <- BiocGenerics::start(genes_sel)
  ord <- order(gene_starts)
  gene_ids <- gene_ids[ord]
  for (gid in gene_ids) {
    gene_gr <- genes_all[genes_all$gene_id == gid]
    gene_strand <- as.character(GenomicRanges::strand(gene_gr))
    exons_gr <- if (gid %in% names(exons_by_gene))
      exons_by_gene[[gid]]
    else
      GenomicRanges::GRanges()
    exons_gr <- suppressWarnings(GenomicRanges::reduce(exons_gr))
    introns_gr <- if (length(exons_gr) > 0)
      suppressWarnings(GenomicRanges::setdiff(gene_gr, exons_gr))
    else
      gene_gr
    tx_names <- if (gid %in% names(tx_by_gene))
      tx_by_gene[[gid]]$tx_name
    else
      character(0)
    utr5_list <- if (length(tx_names) > 0)
      utr5_by_tx[names(utr5_by_tx) %in% tx_names]
    else
      GenomicRanges::GRangesList()
    utr3_list <- if (length(tx_names) > 0)
      utr3_by_tx[names(utr3_by_tx) %in% tx_names]
    else
      GenomicRanges::GRangesList()
    utr5_gr <- suppressWarnings(GenomicRanges::reduce(unlist(utr5_list)))
    utr3_gr <- suppressWarnings(GenomicRanges::reduce(unlist(utr3_list)))
    promoters_gr <- suppressWarnings(
      GenomicFeatures::promoters(
        gene_gr,
        upstream = upstream,
        downstream = downstream,
        use.names = TRUE
      )
    )
    df_g <- rbind(
      to_df(promoters_gr, "promoter", gid, gene_strand),
      to_df(utr5_gr, "utr5", gid, gene_strand),
      to_df(utr3_gr, "utr3", gid, gene_strand),
      to_df(exons_gr, "exon", gid, gene_strand),
      to_df(introns_gr, "intron", gid, gene_strand)
    )
    df_g <- df_g[stats::complete.cases(df_g), , drop = FALSE]
    df_all <- rbind(df_all, df_g)
  }
  gene_levels <- unique(gene_ids)
  df_all$gene_factor <- factor(df_all$gene_id, levels = gene_levels)
  df_all$y <- as.numeric(df_all$gene_factor)
  feature_order <- c("promoter", "utr5", "utr3", "exon", "intron")
  df_all$feature <- factor(df_all$feature, levels = feature_order)
  col_map <- c(
    promoter = promoter_color,
    utr5 = utr5_color,
    utr3 = utr3_color,
    exon = exon_color,
    intron = intron_color
  )
  df_intron_plot <- df_all[df_all$feature == "intron", , drop = FALSE]
  if (nrow(df_intron_plot) > 0) {
    df_intron_plot$x <- ifelse(df_intron_plot$strand == "-",
                               df_intron_plot$end,
                               df_intron_plot$start)
    df_intron_plot$xend <- ifelse(df_intron_plot$strand == "-",
                                  df_intron_plot$start,
                                  df_intron_plot$end)
  }

  # Arrows
  df_intron_arrows <- data.frame()
  if (nrow(df_intron_plot) > 0) {
    for (i in seq_len(nrow(df_intron_plot))) {
      s <- df_intron_plot$start[i]
      e <- df_intron_plot$end[i]
      y <- df_intron_plot$y[i]
      st <- df_intron_plot$strand[i]
      len <- e - s
      if (len <= 0)
        next
      cnt <- if (is.na(arrow_count)) {
        max(1, min(20, floor(len / 1000)))
      } else {
        max(1, as.integer(arrow_count))
      }
      seg_len <- max(1, len / (cnt * 6))
      if (st == "-") {
        pos <- seq(from = e,
                   to = s + seg_len,
                   length.out = cnt)
        x <- pos
        xend <- pos - seg_len
      } else {
        pos <- seq(from = s,
                   to = e - seg_len,
                   length.out = cnt)
        x <- pos
        xend <- pos + seg_len
      }
      df_i <- data.frame(
        x = x,
        xend = xend,
        y = y,
        yend = y,
        feature = "intron"
      )
      df_intron_arrows <- rbind(df_intron_arrows, df_i)
    }
  }
  xmin <- start
  xmax <- end
  xpad <- max(upstream, downstream)
  exon_count <- sum(df_all$feature == "exon")

  # Plot
  p <- ggplot2::ggplot(df_all) +
    ggplot2::geom_rect(
      data = df_all[df_all$feature %in% c("promoter", "utr5", "utr3", "exon"), , drop = FALSE],
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = y - 0.3,
        ymax = y + 0.3,
        fill = feature
      ),
      color = NA,
      alpha = feature_alpha
    ) +
    ggplot2::geom_segment(
      data = df_intron_plot,
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = y,
        color = feature
      ),
      linewidth = intron_width,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = df_intron_arrows,
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = y,
        color = feature
      ),
      linewidth = intron_width,
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(arrow_length, arrow_unit),
        type = "closed"
      )
    ) +
    ggplot2::scale_fill_manual(values = col_map, drop = FALSE) +
    ggplot2::scale_color_manual(values = col_map, drop = FALSE) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(gene_levels),
      labels = gene_levels,
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(
      x = "Genomic position",
      y = "Gene",
      title = paste0(
        "Chrom: ",
        chrom_id,
        "; Interval: ",
        start,
        "-",
        end,
        "; Genes: ",
        length(gene_levels),
        "; Exons: ",
        exon_count,
        "; Introns: ",
        sum(df_all$feature == "intron")
      )
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = x_breaks)) +
    ggplot2::coord_cartesian(xlim = c(xmin - xpad, xmax + xpad)) +
    my_theme()
  p
}
