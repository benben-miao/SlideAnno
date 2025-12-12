#' @title Chromosome-level signal overlay (FST and DMR)
#' @description Overlay \bold{\emph{FST}} and \bold{\emph{DMR}} signals along a selected chromosome with optional smoothing.
#'
#' @return A \bold{\emph{ggplot object}} or a \bold{\emph{patchwork}} of two aligned panels.
#' @param chrom_id Chromosome identifier (e.g., \bold{\emph{"chr1"}}).
#' @param fst_file FST slide window table (\bold{\emph{tab-separated}}). (\bold{\emph{NULL}} for skip).
#' @param dmr_file DMR table from \bold{\emph{MethylKit}}. (\bold{\emph{NULL}} for skip).
#' @param smooth_span Span for local regression smoothing. (\bold{\emph{0.1}}).
#' @param fst_color Line/point color for FST. (\bold{\emph{"#0088ff"}}).
#' @param hyper_color Color for hyper-methylated DMRs. (\bold{\emph{"#ff000055"}}).
#' @param hypo_color Color for hypo-methylated DMRs. (\bold{\emph{"#00880055"}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example FST and DMR tables in GAnnoViz
#' fst_file <- system.file(
#'   "extdata",
#'   "example.fst",
#'   package = "GAnnoViz")
#'
#' dmr_file <- system.file(
#'   "extdata",
#'   "example.dmr",
#'   package = "GAnnoViz")
#'
#' # Overlay: FST on chr11
#' plot_chrom_signal_overlay(
#'   chrom_id = "chr11",
#'   fst_file = fst_file,
#'   dmr_file = NULL,
#'   smooth_span = 0.2,
#'   fst_color = "#0088ff"
#' )
#'
#' # Overlay: DMR on chr1
#' plot_chrom_signal_overlay(
#'   chrom_id = "chr1",
#'   fst_file = NULL,
#'   dmr_file = dmr_file,
#'   smooth_span = 0.2,
#'   hyper_color = "#ff000055",
#'   hypo_color = "#00880055"
#' )
plot_chrom_signal_overlay <- function(chrom_id,
                                      fst_file = NULL,
                                      dmr_file = NULL,
                                      smooth_span = 0.1,
                                      fst_color = "#0088ff",
                                      hyper_color = "#ff000055",
                                      hypo_color = "#00880055") {
  if (is.null(fst_file) && is.null(dmr_file))
    stop("Provide at least one of fst_file or dmr_file")

  p_list <- list()

  if (!is.null(fst_file)) {
    fst <- utils::read.table(fst_file, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
    fst <- fst %>% dplyr::filter(CHROM == chrom_id)
    if (nrow(fst) > 0) {
      p_fst <- ggplot2::ggplot(fst, ggplot2::aes(x = BIN_START, y = WEIGHTED_FST)) +
        ggplot2::geom_line(color = fst_color, linewidth = 0.6) +
        ggplot2::geom_point(color = fst_color, size = 0.6, alpha = 0.6) +
        ggplot2::geom_smooth(span = smooth_span, se = FALSE, color = fst_color) +
        ggplot2::labs(x = "Genomic position", y = "WEIGHTED FST") +
        my_theme()
      p_list[[length(p_list) + 1]] <- p_fst
    }
  }

  if (!is.null(dmr_file)) {
    dmr <- utils::read.table(dmr_file, header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE)
    dmr <- dmr %>% dplyr::filter(chr == chrom_id)
    if (nrow(dmr) > 0) {
      dmr$state <- ifelse(dmr$meth.diff >= 0, "hyper", "hypo")
      p_dmr <- ggplot2::ggplot(dmr, ggplot2::aes(x = start, y = meth.diff, color = state)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.6) +
        ggplot2::geom_smooth(span = smooth_span, se = FALSE, color = "#444444") +
        ggplot2::scale_color_manual(values = c(hyper = hyper_color, hypo = hypo_color)) +
        ggplot2::labs(x = "Genomic position", y = "Methylation diff (%)", color = NULL) +
        my_theme()
      p_list[[length(p_list) + 1]] <- p_dmr
    }
  }

  if (length(p_list) == 0)
    stop("No data on the selected chromosome")
  if (length(p_list) == 1) {
    p_list[[1]]
  } else {
    patchwork::wrap_plots(p_list, ncol = 1)
  }
}

