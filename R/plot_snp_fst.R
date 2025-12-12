#' @title Plot genomic weighted FST heatmap
#' @description Plot genomic \bold{\emph{weighted FST}} heatmap.
#' @author benben-miao
#'
#' @return A \bold{\emph{ggplot object}} showing genomic weighted FST heatmap.
#' @param fst_file FST sliding window results. (\bold{\emph{CHROM, BIN_START, BIN_END, WEIGHTED_FST, N_VARIANTS}}.
#' @param bin_size Bin size in base pairs. (\bold{\emph{1e6}}).
#' @param metric Aggregation metric for bin fill. (\bold{\emph{"fst_mean"}}, "variant_count").
#' @param orientation Coordinate orientation. (\bold{\emph{"horizontal"}}, "vertical").
#' @param palette Continuous color for weighted FST. (\bold{\emph{c("#ffffff", "#aa00aa")}}).
#' @param alpha Tile alpha for heatmap. (\bold{\emph{0.9}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example FST sliding window table in GAnnoViz
#' fst_file <- system.file(
#'   "extdata",
#'   "example.fst",
#'   package = "GAnnoViz")
#'
#' # Plot weighted FST
#' plot_snp_fst(
#'   fst_file = fst_file,
#'   bin_size = 1e6,
#'   metric = "fst_mean",
#'   orientation = "horizontal",
#'   palette = c("#ffffff", "#aa00aa"),
#'   alpha = 0.9
#' )
plot_snp_fst <- function(fst_file,
                         bin_size = 1e6,
                         metric = "fst_mean",
                         orientation = "horizontal",
                         palette = c("#ffffff", "#aa00aa"),
                         alpha = 0.9) {
  # Parameter check
  if (!(metric %in% c("fst_mean", "variant_count")))
    stop("metric must be one of: 'fst_mean','variant_count'")

  # FST results
  fst <- utils::read.table(
    fst_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    stringsAsFactors = FALSE
  )

  # Chromosomes
  chrom_raw <- as.character(fst$CHROM)
  chr_prefix <- any(stringr::str_detect(chrom_raw, "^(?i)chr"))
  keep_mask <- if (chr_prefix)
    stringr::str_detect(chrom_raw, "^(?i)chr")
  else
    rep(TRUE, length(chrom_raw))
  fst <- fst[keep_mask, , drop = FALSE]

  suppressWarnings({
    chrom_num <- as.numeric(stringr::str_extract(as.character(fst$CHROM), "\\d+"))
  })
  order_index <- order(is.na(chrom_num), chrom_num, fst$CHROM)
  chrom_levels <- unique(as.character(fst$CHROM)[order_index])

  # Bin analysis
  agg_one_chrom <- function(ch) {
    sub <- fst[fst$CHROM == ch, , drop = FALSE]
    if (nrow(sub) == 0)
      return(NULL)
    max_end <- max(sub$BIN_END, na.rm = TRUE)
    bin_starts <- seq(1, max_end, by = bin_size)
    bin_ends <- pmin(bin_starts + bin_size - 1, max_end)
    bins <- data.frame(
      chrom = ch,
      start = bin_starts,
      end = bin_ends,
      stringsAsFactors = FALSE
    )
    res <- lapply(seq_len(nrow(bins)), function(i) {
      bstart <- bins$start[i]
      bend <- bins$end[i]
      mask <- (sub$BIN_START <= bend) & (sub$BIN_END >= bstart)
      if (!any(mask))
        return(list(value = NA_real_, center = (bstart + bend) / 2))
      s <- sub[mask, , drop = FALSE]
      ov_len <- pmin(bend, s$BIN_END) - pmax(bstart, s$BIN_START) + 1
      ov_len[ov_len < 0] <- 0
      if (metric == "fst_mean") {
        w <- ov_len
        wsum <- sum(w, na.rm = TRUE)
        val <- if (wsum > 0)
          sum(s$WEIGHTED_FST * w, na.rm = TRUE) / wsum
        else
          NA_real_
      } else {
        val <- sum(s$N_VARIANTS[ov_len > 0], na.rm = TRUE)
      }
      list(value = val, center = (bstart + bend) / 2)
    })
    bins$value <- vapply(res, function(x)
      x$value, numeric(1))
    bins$center <- vapply(res, function(x)
      x$center, numeric(1))
    bins
  }

  df_bin <- do.call(rbind, lapply(chrom_levels, agg_one_chrom))
  if (is.null(df_bin) ||
      nrow(df_bin) == 0)
    stop("No bins computed from FST table")

  length_map <- stats::aggregate(end ~ chrom, data = df_bin, FUN = max)
  names(length_map) <- c("chrom", "length")
  x_max <- max(length_map$length, na.rm = TRUE)

  fill_label <- if (metric == "fst_mean")
    "Weighted FST"
  else
    "SNP count"

  # Plot
  p <- ggplot2::ggplot(df_bin, ggplot2::aes(x = center, y = factor(chrom, levels = chrom_levels))) +
    ggplot2::geom_tile(ggplot2::aes(
      width = end - start,
      height = 0.9,
      fill = value
    ),
    alpha = alpha) +
    ggplot2::scale_fill_gradientn(colors = palette, na.value = "#ffffff00") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::label_number(),
                                breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::labs(x = "Genomic position", y = "Chromosome", fill = fill_label) +
    ggplot2::coord_cartesian(xlim = c(0, x_max * 1.05)) +
    my_theme()
  if (orientation == "horizontal")
    p <- p + ggplot2::coord_flip()
  p
}