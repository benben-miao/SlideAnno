#' @title Plot DMGs manhattan plot across chromosomes
#' @description Plot genome-wide \bold{\emph{DMGs}} as a manhattan-style scatter plot, using chromosomes on x-axis and \code{meth.diff} on y-axis.
#' @author benben-miao
#'
#' @return A \bold{\emph{ggplot object}} of DMG manhattan plot.
#' @param dmr_file DEG table from \bold{\emph{MethylKit}} analysis.
#' @param y_transform Y-axis transformation for \code{meth.diff}. (\bold{\emph{"none"}}, "log2", "log10"). Log transform uses \code{sign(x) * log(1 + abs(x), base)}.
#' @param label_col Column name used for annotation labels. If \code{NULL}, uses \code{chr:start-end}.
#' @param chromosome_spacing Gap width (bp) inserted between chromosomes on x-axis. (\bold{\emph{1e6}}).
#' @param gff_file Genomic structural annotation \bold{\file{GFF3/GTF}} file path. If provided and \code{label_col} is \code{NULL}, will try to label by overlapped/nearest \code{gene_id}.
#' @param format Format of GFF3/GTF file. (\bold{\emph{"auto"}}, "gff3", "gtf").
#' @param label_type Label by gene name or gene id. (\bold{\emph{"name"}}, "id"). If \code{"name"} but no gene names can be inferred from \code{gff_file}, please provide \code{gene_table}.
#' @param gene_table Optional gene ID/name mapping table (first two columns: id, name). If provided, will use gene name for labels when possible.
#' @param label_top_n Number of top positive and top negative DMGs to label. (\bold{\emph{10}}).
#' @param hyper_color Color for hyper-methylated points. (\bold{\emph{"#ff0000"}}).
#' @param hypo_color Color for hypo-methylated points. (\bold{\emph{"#008800"}}).
#' @param point_size Point size. (\bold{\emph{1}}).
#' @param point_alpha Point alpha. (\bold{\emph{0.5}}).
#' @param label_size Text size for labels. (\bold{\emph{3}}).
#' @param connector_dx1 First connector horizontal offset (bp). Default adapts to genome size.
#' @param connector_dx2 Second connector horizontal offset (bp). Default adapts to genome size.
#' @param connector_elbow Scale factor applied to \code{connector_dx2}. (\bold{\emph{0.8}}).
#' @param connector_tilt_frac Tilt amplitude for the second connector segment as a fraction of \code{gap_frac}. (\bold{\emph{0.2}}).
#' @param gap_frac Minimum vertical gap between labels (fraction of y-range). (\bold{\emph{0.04}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' dmr_file <- system.file(
#'   "extdata",
#'   "example.dmr",
#'   package = "GAnnoViz")
#'
#' plot_dmg_manhattan(
#'   dmr_file = dmr_file,
#'   y_transform = "none",
#'   label_top_n = 10,
#'   point_size = 1,
#'   point_alpha = 0.5,
#'   hyper_color = "#ff0000",
#'   hypo_color = "#008800"
#' )
#'
plot_dmg_manhattan <- function(dmr_file,
                               y_transform = c("none", "log2", "log10"),
                               label_col = NULL,
                               label_top_n = 10,
                               hyper_color = "#ff0000",
                               hypo_color = "#008800",
                               point_size = 1,
                               point_alpha = 0.5,
                               label_size = 3,
                               connector_dx1 = NULL,
                               connector_dx2 = NULL,
                               gap_frac = 0.04,
                               chromosome_spacing = 1e6,
                               gff_file = NULL,
                               format = "auto",
                               gene_table = NULL,
                               label_type = c("name", "id"),
                               connector_elbow = 0.8,
                               connector_tilt_frac = 0.2) {
  y_transform <- match.arg(y_transform)
  label_type <- match.arg(label_type)
  dmr <- utils::read.table(
    dmr_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cols <- colnames(dmr)
  needed <- c("chr", "start", "meth.diff")
  missing <- setdiff(needed, cols)
  if (length(missing) > 0)
    stop(paste0("Missing required columns: ", paste(missing, collapse = ", ")))
  dmr <- dmr[stats::complete.cases(dmr$chr, dmr$start, dmr$meth.diff), , drop = FALSE]
  if (nrow(dmr) == 0)
    stop("No DMGs with valid positions found")

  chrom_raw <- as.character(dmr$chr)
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  dmr <- dmr[keep_mask, , drop = FALSE]
  if (nrow(dmr) == 0)
    stop("No DMGs found on selected chromosomes")

  dmr$end <- if ("end" %in% colnames(dmr)) dmr$end else dmr$start
  end_col <- dmr$end
  dmr$pos <- (dmr$start + dmr$end) / 2

  chrom_chr <- as.character(dmr$chr)
  chrom_num <- stringr::str_extract(chrom_chr, "\\d+")
  suppressWarnings({
    chrom_num <- as.numeric(chrom_num)
  })
  order_index <- order(is.na(chrom_num), chrom_num, chrom_chr)
  chrom_levels <- unique(chrom_chr[order_index])

  length_map <- stats::aggregate(end_col ~ chr,
                                 data = data.frame(chr = dmr$chr, end_col = end_col),
                                 FUN = max)
  length_map <- length_map[match(chrom_levels, length_map$chr), , drop = FALSE]
  chrom_gap <- suppressWarnings(as.numeric(chromosome_spacing))
  if (!is.finite(chrom_gap) || chrom_gap < 0)
    chrom_gap <- 0
  length_map$cum_start <- c(0, cumsum(utils::head(length_map$end_col + chrom_gap, -1)))
  length_map$mid <- length_map$cum_start + length_map$end_col / 2
  dmr$x <- dmr$pos + length_map$cum_start[match(dmr$chr, length_map$chr)]

  y_fun <- switch(
    y_transform,
    none = function(x) x,
    log2 = function(x) sign(x) * log(1 + abs(x), base = 2),
    log10 = function(x) sign(x) * log10(1 + abs(x))
  )
  dmr$y <- y_fun(dmr$meth.diff)
  dmr$state <- ifelse(dmr$meth.diff >= 0, "hyper", "hypo")
  col_map <- c(hyper = hyper_color, hypo = hypo_color)

  gene_id_mapped <- NULL
  label_vec <- NULL
  if (!is.null(label_col) &&
    is.character(label_col) &&
    length(label_col) == 1 &&
    label_col %in% colnames(dmr)) {
    label_vec <- as.character(dmr[[label_col]])
  } else {
    if (label_type == "name") {
      auto_cols <- c("gene_name", "GeneName", "SYMBOL", "symbol", "gene")
      hit <- auto_cols[auto_cols %in% colnames(dmr)][1]
      if (!is.null(hit) && !is.na(hit))
        label_vec <- as.character(dmr[[hit]])
      if (is.null(label_vec)) {
        auto_cols <- c("gene_id", "GeneID", "ID", "gene")
        hit <- auto_cols[auto_cols %in% colnames(dmr)][1]
        if (!is.null(hit) && !is.na(hit)) {
          gene_id_mapped <- as.character(dmr[[hit]])
          label_vec <- gene_id_mapped
        }
      }
    }
    if (is.null(label_vec) && label_type == "id") {
      auto_cols <- c("gene_id", "GeneID", "ID", "gene")
      hit <- auto_cols[auto_cols %in% colnames(dmr)][1]
      if (!is.null(hit) && !is.na(hit))
        label_vec <- as.character(dmr[[hit]])
    }
  }
  if (is.null(label_vec) && !is.null(gff_file)) {
    genes <- extract_genes(gff_file = gff_file, format = format, gene_info = "all")
    genes <- genes[as.character(GenomicRanges::seqnames(genes)) %in% chrom_levels]
    if (length(genes) > 0) {
      gr_dmr <- GenomicRanges::GRanges(
        seqnames = dmr$chr,
        ranges = IRanges::IRanges(
          start = as.integer(dmr$start),
          end = as.integer(dmr$end)
        )
      )
      hits <- GenomicRanges::findOverlaps(gr_dmr, genes, ignore.strand = TRUE)
      gene_id <- rep(NA_character_, length(gr_dmr))
      if (length(hits) > 0) {
        q <- S4Vectors::queryHits(hits)
        s <- S4Vectors::subjectHits(hits)
        keep_first <- !duplicated(q)
        q <- q[keep_first]
        s <- s[keep_first]
        gene_id[q] <- as.character(genes$gene_id[s])
      }
      na_idx <- which(is.na(gene_id) | gene_id == "")
      if (length(na_idx) > 0) {
        near <- GenomicRanges::nearest(gr_dmr[na_idx], genes, ignore.strand = TRUE)
        ok <- !is.na(near)
        if (any(ok))
          gene_id[na_idx[ok]] <- as.character(genes$gene_id[near[ok]])
      }
      gene_id_mapped <- gene_id
      label_vec <- gene_id
    }
  }
  if (label_type == "name" && !is.null(gene_id_mapped)) {
    gene_name <- rep(NA_character_, length(gene_id_mapped))
    if (!is.null(gene_table)) {
      gene_table <- as.data.frame(gene_table, stringsAsFactors = FALSE)
      if (ncol(gene_table) >= 2) {
        names(gene_table)[1:2] <- c("gene_id", "gene_name")
        idx <- match(gene_id_mapped, gene_table$gene_id)
        ok <- !is.na(idx)
        if (any(ok))
          gene_name[ok] <- as.character(gene_table$gene_name[idx[ok]])
      }
    }
    if (all(is.na(gene_name) | gene_name == "") && !is.null(gff_file)) {
      infer_gene_name_map <- function(gff_file, ids, fmt) {
        ids <- unique(as.character(ids))
        ids <- ids[!is.na(ids) & ids != ""]
        if (length(ids) == 0)
          return(stats::setNames(character(0), character(0)))
        ids_norm <- sub("^.+?:", "", ids)
        wanted <- unique(ids_norm)
        map <- stats::setNames(rep(NA_character_, length(wanted)), wanted)
        con <- if (grepl("\\.gz$", gff_file, ignore.case = TRUE))
          gzfile(gff_file, open = "rt")
        else
          file(gff_file, open = "rt")
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        repeat {
          lines <- readLines(con, n = 50000, warn = FALSE)
          if (length(lines) == 0)
            break
          lines <- lines[!startsWith(lines, "#")]
          if (length(lines) == 0)
            next
          parts <- strsplit(lines, "\t", fixed = TRUE)
          attrs <- vapply(parts, function(z) if (length(z) >= 9) z[[9]] else NA_character_, character(1))
          attrs <- attrs[!is.na(attrs)]
          if (length(attrs) == 0)
            next
          if (identical(fmt, "gtf")) {
            gid <- stringr::str_match(attrs, 'gene_id\\s+"([^"]+)"')[, 2]
            gnm <- stringr::str_match(attrs, 'gene_name\\s+"([^"]+)"')[, 2]
            ok <- !is.na(gid) & gid != "" & !is.na(gnm) & gnm != ""
            if (any(ok)) {
              gid_norm <- sub("^.+?:", "", gid[ok])
              hit <- gid_norm %in% wanted
              if (any(hit)) {
                idx <- match(gid_norm[hit], names(map))
                fill <- is.na(map[idx]) | map[idx] == ""
                if (any(fill))
                  map[idx[fill]] <- gnm[ok][hit][fill]
              }
            }
          } else {
            gid <- stringr::str_match(attrs, "(?:^|;)ID=([^;]+)")[, 2]
            gnm <- stringr::str_match(attrs, "(?:^|;)Name=([^;]+)")[, 2]
            if (all(is.na(gnm) | gnm == "")) {
              gnm <- stringr::str_match(attrs, "(?:^|;)gene_name=([^;]+)")[, 2]
            }
            ok <- !is.na(gid) & gid != "" & !is.na(gnm) & gnm != ""
            if (any(ok)) {
              gid_norm <- sub("^.+?:", "", gid[ok])
              hit <- gid_norm %in% wanted
              if (any(hit)) {
                idx <- match(gid_norm[hit], names(map))
                fill <- is.na(map[idx]) | map[idx] == ""
                if (any(fill))
                  map[idx[fill]] <- gnm[ok][hit][fill]
              }
            }
          }
          if (all(!is.na(map) & map != ""))
            break
        }
        map[!is.na(map) & map != ""]
      }

      fmt <- resolve_gff_format(gff_file, format)
      gmap <- infer_gene_name_map(gff_file, gene_id_mapped, fmt)
      if (length(gmap) > 0) {
        gid_norm <- sub("^.+?:", "", as.character(gene_id_mapped))
        idx <- match(gid_norm, names(gmap))
        ok <- !is.na(idx)
        if (any(ok))
          gene_name[ok] <- unname(gmap[idx[ok]])
      }
    }
    label_vec <- ifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id_mapped)
  }
  if (is.null(label_vec)) {
    end_str <- if ("end" %in% colnames(dmr)) as.character(dmr$end) else as.character(dmr$start)
    label_vec <- paste0(dmr$chr, ":", dmr$start, "-", end_str)
  } else {
    end_str <- if ("end" %in% colnames(dmr)) as.character(dmr$end) else as.character(dmr$start)
    interval_label <- paste0(dmr$chr, ":", dmr$start, "-", end_str)
    empty <- is.na(label_vec) | label_vec == ""
    label_vec[empty] <- interval_label[empty]
  }
  dmr$label <- label_vec

  n_label <- max(0, as.integer(label_top_n))
  df_pos <- dmr %>%
    dplyr::filter(is.finite(meth.diff) & meth.diff > 0) %>%
    dplyr::arrange(dplyr::desc(meth.diff)) %>%
    utils::head(n_label)
  df_neg <- dmr %>%
    dplyr::filter(is.finite(meth.diff) & meth.diff < 0) %>%
    dplyr::arrange(meth.diff) %>%
    utils::head(n_label)
  df_top <- dplyr::bind_rows(df_pos, df_neg)

  p_base <- ggplot2::ggplot(dmr, ggplot2::aes(x = x, y = y, color = state)) +
    ggplot2::geom_hline(yintercept = 0, color = "#333333", linewidth = 0.4) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = col_map, drop = FALSE, name = NULL, labels = c(hypo = "Hypo", hyper = "Hyper")) +
    ggplot2::scale_x_continuous(
      breaks = length_map$mid,
      labels = length_map$chr,
      expand = ggplot2::expansion(mult = c(0.01, 0.05))
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::labs(x = "Chromosome", y = "Meth Diff") +
    my_theme()

  chr_bound <- length_map$cum_start[-1]
  if (length(chr_bound) > 0) {
    p_base <- p_base + ggplot2::geom_vline(xintercept = chr_bound, color = "#999999", linewidth = 0.3)
  }

  if (nrow(df_top) == 0)
    return(p_base)

  y_min <- min(dmr$y, na.rm = TRUE)
  y_max <- max(dmr$y, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range == 0)
    y_range <- 1

  genome_len <- sum(length_map$end_col, na.rm = TRUE) + chrom_gap * max(0, nrow(length_map) - 1)
  dx1 <- if (is.null(connector_dx1)) genome_len * 0.005 else connector_dx1
  dx2 <- if (is.null(connector_dx2)) genome_len * 0.01 else connector_dx2
  elbow <- suppressWarnings(as.numeric(connector_elbow))
  if (!is.finite(elbow) || elbow <= 0)
    elbow <- 0.8
  elbow <- min(1, max(0.05, elbow))
  dx2 <- dx2 * elbow
  min_gap <- y_range * gap_frac
  tilt_frac <- suppressWarnings(as.numeric(connector_tilt_frac))
  if (!is.finite(tilt_frac) || tilt_frac < 0)
    tilt_frac <- 0
  tilt_amp <- min_gap * tilt_frac

  df_pos <- df_top %>% dplyr::filter(meth.diff > 0)
  df_neg <- df_top %>% dplyr::filter(meth.diff < 0)

  df_pos <- df_pos[order(df_pos$x), , drop = FALSE]
  if (nrow(df_pos) > 0) {
    base_y_pos <- y_max + y_range * 0.08
    cur <- base_y_pos
    df_pos$y_label <- NA_real_
    for (i in seq_len(nrow(df_pos))) {
      df_pos$y_label[i] <- cur
      cur <- cur + min_gap
    }
    idx <- seq_len(nrow(df_pos)) - 1
    df_pos$tilt <- ((idx %% 3L) - 1L) * tilt_amp
  }

  df_neg <- df_neg[order(df_neg$x), , drop = FALSE]
  if (nrow(df_neg) > 0) {
    base_y_neg <- y_min - y_range * 0.08
    cur <- base_y_neg
    df_neg$y_label <- NA_real_
    for (i in seq_len(nrow(df_neg))) {
      df_neg$y_label[i] <- cur
      cur <- cur - min_gap
    }
    idx <- seq_len(nrow(df_neg)) - 1
    df_neg$tilt <- ((idx %% 3L) - 1L) * tilt_amp
  }

  df_lab <- dplyr::bind_rows(df_pos, df_neg)
  df_lab$x1 <- df_lab$x
  df_lab$y1 <- df_lab$y
  df_lab$x2 <- df_lab$x1 + dx1
  df_lab$y2 <- df_lab$y_label
  df_lab$x3 <- df_lab$x2 + dx2
  df_lab$y3 <- df_lab$y_label + ifelse(is.na(df_lab$tilt), 0, df_lab$tilt)

  df_seg1 <- data.frame(x = df_lab$x1, y = df_lab$y1, xend = df_lab$x2, yend = df_lab$y2)
  df_seg2 <- data.frame(x = df_lab$x2, y = df_lab$y2, xend = df_lab$x3, yend = df_lab$y3)

  p_base +
    ggplot2::geom_segment(
      data = df_seg1,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "#333333",
      linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = df_seg2,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "#333333",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = df_lab,
      ggplot2::aes(x = x3, y = y3, label = label),
      size = label_size,
      color = "#000000",
      hjust = 0,
      vjust = 0.5
    )
}
