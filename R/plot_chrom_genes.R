#' @title Plot chromosome structures and gene annotation
#' @description Plot \bold{\emph{chromosome structures}} and gene annotation.
#'
#' @return A ggplot object of chromosomes with gene labels.
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path.
#' @param gene_table A data.frame: first column gene ID, second column gene name.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param annotate Which to annotate on the right: (\bold{\emph{"id"}}, "name").
#' @param orientation Coordinate orientation: (\bold{\emph{"vertical"}}, "horizontal").
#' @param bar_width Chromosome bars relative width. (\bold{\emph{0.6}}).
#' @param chrom_alpha Chromosome bars alpha. (\bold{\emph{0.1}}).
#' @param gene_width Gene bar relative width. (\bold{\emph{0.5}}).
#' @param chrom_color Chromosome bar color. (\bold{\emph{"#008888"}}).
#' @param gene_color Gene rectangle color. (\bold{\emph{"#0088ff"}}).
#' @param telomere_color Telomere color. (\bold{\emph{"#ff0000"}}).
#' @param label_size Label text size. (\bold{\emph{3}}).
#' @param connector_dx1 First connector horizontal offset. (\bold{\emph{0.15}}).
#' @param connector_dx2 Second connector horizontal offset. (\bold{\emph{0.15}}).
#' @param min_gap_frac Minimum vertical gap between labels (fraction of chromosome length). (\bold{\emph{0.02}}).
#'
#' @export
#'
#' @examples
#' # Example GFF3 file in SlideAnno
#' gff_file <- system.file(
#'   "extdata",
#'   "example.gff",
#'   package = "SlideAnno")
#'
#' genes <- data.frame(
#'   gene_id = c("HdF029609", "HdF029610"),
#'   gene_name = c("GeneA", "GeneB"))
#'
#' # Vertical, annotate by name
#' plot_chrom_genes(
#'   gff_file = gff_file,
#'   gene_table = genes,
#'   annotate = "name",
#'   orientation = "vertical")
#'
#' # Horizontal, annotate by id
#' plot_chrom_genes(
#'   gff_file = gff_file,
#'   gene_table = genes,
#'   annotate = "id",
#'   orientation = "horizontal")
#'
plot_chrom_genes <- function(gff_file,
                             gene_table,
                             format = "auto",
                             annotate = "id",
                             orientation = "vertical",
                             bar_width = 0.6,
                             chrom_alpha = 0.1,
                             gene_width = 0.5,
                             chrom_color = "#008888",
                             gene_color = "#0088ff",
                             telomere_color = "#ff0000",
                             label_size = 3,
                             connector_dx1 = 0.2,
                             connector_dx2 = 0.2,
                             min_gap_frac = 0.02) {
  # Genes
  genes <- extract_genes(gff_file = gff_file,
                         format = format,
                         gene_info = "all")

  # Chrom
  chrom_raw <- as.character(GenomicRanges::seqnames(genes))
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  chrom_chr <- chrom_raw[keep_mask]
  genes <- genes[keep_mask]
  gene_ids <- genes$gene_id
  chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
  suppressWarnings({
    chrom_num <- as.numeric(chrom_num)
  })
  order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
  chrom_levels <- unique(chrom_chr[order_index])

  # Genes data.frame
  df_genes <- data.frame(
    chrom = chrom_chr,
    start = BiocGenerics::start(genes),
    end = BiocGenerics::end(genes),
    gene_id = gene_ids,
    stringsAsFactors = FALSE
  )
  df_genes <- df_genes[stats::complete.cases(df_genes), , drop = FALSE]
  df_genes$pos <- (df_genes$start + df_genes$end) / 2

  length_map <- stats::aggregate(end ~ chrom, data = df_genes, FUN = max)
  names(length_map) <- c("chrom", "length")

  x_idx <- seq_along(chrom_levels)
  pos_map <- data.frame(chrom = chrom_levels,
                        x = x_idx,
                        stringsAsFactors = FALSE)

  df_chrom <- merge(length_map, pos_map, by = "chrom", all.x = TRUE)
  df_chrom <- df_chrom[stats::complete.cases(df_chrom), , drop = FALSE]
  df_chrom$xmin <- df_chrom$x - bar_width / 2
  df_chrom$xmax <- df_chrom$x + bar_width / 2
  df_chrom$ymin <- 0
  df_chrom$ymax <- df_chrom$length

  # Chrom roundangle
  n_arc <- 50
  poly_list <- lapply(seq_len(nrow(df_chrom)), function(i) {
    ch <- df_chrom$chrom[i]
    x <- df_chrom$x[i]
    len <- df_chrom$length[i]
    xL <- x - bar_width / 2
    xR <- x + bar_width / 2
    yB <- 0
    yT <- len
    rx <- bar_width * 0.45
    ry <- min(max(1, 0.04 * len), len / 2)
    t1 <- seq(-pi / 2, 0, length.out = n_arc)
    t2 <- seq(0, pi / 2, length.out = n_arc)
    t3 <- seq(pi / 2, pi, length.out = n_arc)
    t4 <- seq(pi, 3 * pi / 2, length.out = n_arc)
    bx_x <- seq(xL + rx, xR - rx, length.out = n_arc)
    bx_y <- rep(yB, n_arc)
    br_x <- (xR - rx) + rx * cos(t1)
    br_y <- (yB + ry) + ry * sin(t1)
    re_x <- rep(xR, n_arc)
    re_y <- seq(yB + ry, yT - ry, length.out = n_arc)
    tr_x <- (xR - rx) + rx * cos(t2)
    tr_y <- (yT - ry) + ry * sin(t2)
    tx_x <- seq(xR - rx, xL + rx, length.out = n_arc)
    tx_y <- rep(yT, n_arc)
    tl_x <- (xL + rx) + rx * cos(t3)
    tl_y <- (yT - ry) + ry * sin(t3)
    le_x <- rep(xL, n_arc)
    le_y <- seq(yT - ry, yB + ry, length.out = n_arc)
    bl_x <- (xL + rx) + rx * cos(t4)
    bl_y <- (yB + ry) + ry * sin(t4)
    xs <- c(bx_x, br_x, re_x, tr_x, tx_x, tl_x, le_x, bl_x)
    ys <- c(bx_y, br_y, re_y, tr_y, tx_y, tl_y, le_y, bl_y)
    data.frame(chrom = ch, xp = xs, yp = ys)
  })
  df_chrom_poly <- do.call(rbind, poly_list)

  df_gene_rect <- merge(df_genes[, c("chrom", "start", "end", "pos")], pos_map, by = "chrom", all.x = TRUE)
  df_gene_rect <- df_gene_rect[stats::complete.cases(df_gene_rect), , drop = FALSE]
  df_gene_rect$xmin <- df_gene_rect$x - (bar_width * gene_width) / 2
  df_gene_rect$xmax <- df_gene_rect$x + (bar_width * gene_width) / 2
  df_gene_rect$ymin <- pmax(0, df_gene_rect$start)
  df_gene_rect$ymax <- pmax(df_gene_rect$ymin, df_gene_rect$end)

  df_tel <- df_chrom[, c("chrom", "length", "x")]
  df_tel$tw <- pmax(1, 0.02 * df_tel$length)
  df_tel_b <- data.frame(
    chrom = df_tel$chrom,
    x = df_tel$x,
    xmin = df_tel$x - (bar_width * 0.5) / 2,
    xmax = df_tel$x + (bar_width * 0.5) / 2,
    ymin = 0,
    ymax = df_tel$tw
  )
  df_tel_t <- data.frame(
    chrom = df_tel$chrom,
    x = df_tel$x,
    xmin = df_tel$x - (bar_width * 0.5) / 2,
    xmax = df_tel$x + (bar_width * 0.5) / 2,
    ymin = df_tel$length - df_tel$tw,
    ymax = df_tel$length
  )
  y_max <- max(df_chrom$length, na.rm = TRUE)

  # Gene id or name
  gene_table <- as.data.frame(gene_table, stringsAsFactors = FALSE)
  if (ncol(gene_table) < 2)
    stop("gene_table must have at least 2 columns: id and name")
  names(gene_table)[1:2] <- c("gene_id", "gene_name")
  df_points <- merge(df_genes[, c("chrom", "pos", "gene_id")], gene_table, by = "gene_id", all.x = FALSE)
  df_points <- merge(df_points, pos_map, by = "chrom", all.x = TRUE)
  df_points <- df_points[stats::complete.cases(df_points), , drop = FALSE]
  if (nrow(df_points) == 0)
    stop("No gene_id matched between GFF/GTF and gene_table first column.")
  df_points$label <- if (annotate == "id")
    df_points$gene_id
  else
    df_points$gene_name
  min_gap <- max(1, min_gap_frac * y_max)
  df_points <- df_points[order(df_points$x, df_points$pos), , drop = FALSE]
  adjusted_list <- split(df_points, df_points$chrom)
  adjusted_list <- lapply(adjusted_list, function(d) {
    if (nrow(d) == 0)
      return(d)
    ord <- order(d$pos)
    d <- d[ord, , drop = FALSE]
    y <- d$pos
    for (i in seq_along(y)) {
      if (i == 1)
        next
      if ((y[i] - y[i - 1]) < min_gap)
        y[i] <- y[i - 1] + min_gap
    }
    d$y_label <- y
    d
  })
  df_points2 <- do.call(rbind, adjusted_list)
  df_points2 <- merge(df_points2,
                      df_chrom[, c("chrom", "x", "xmin", "xmax")],
                      by = c("chrom", "x"),
                      all.x = TRUE)
  df_points2$xR <- df_points2$xmax
  df_points2$x1 <- df_points2$xR + connector_dx1
  df_points2$x2 <- df_points2$x1
  df_points2$x3 <- df_points2$x1 + connector_dx2
  df_points2$gid <- df_points2$gene_id
  df_poly <- do.call(rbind, lapply(seq_len(nrow(df_points2)), function(i) {
    r <- df_points2[i, ]
    data.frame(
      gid = r$gid,
      x = c(r$xR, r$x1, r$x2, r$x3),
      y = c(r$pos, r$pos, r$y_label, r$y_label)
    )
  }))

  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = df_chrom_poly,
      ggplot2::aes(x = xp, y = yp, group = chrom),
      fill = chrom_color,
      color = NA,
      alpha = chrom_alpha
    ) +
    ggplot2::geom_rect(
      data = df_gene_rect,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = gene_color,
      color = NA
    ) +
    ggplot2::geom_rect(
      data = df_tel_b,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = telomere_color,
      color = NA,
      alpha = 0.5
    ) +
    ggplot2::geom_rect(
      data = df_tel_t,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = telomere_color,
      color = NA,
      alpha = 0.5
    ) +
    ggplot2::geom_path(data = df_poly,
                       ggplot2::aes(x = x, y = y, group = gid),
                       color = "#333333") +
    ggplot2::geom_text(
      data = df_points2,
      ggplot2::aes(x = x3, y = y_label, label = label),
      size = label_size,
      color = "#000000",
      hjust = 0,
      vjust = 0.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_idx,
      labels = chrom_levels,
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_number(),
                                expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::labs(
      x = "Length (bp)",
      y = "Chromosome"
    ) +
    ggplot2::coord_cartesian(ylim = c(
      0,
      max(df_points2$y_label, df_chrom$length, na.rm = TRUE) * 1.05
    )) +
    my_theme()
  if (orientation == "horizontal")
    p <- p + ggplot2::coord_flip()
  p
}
