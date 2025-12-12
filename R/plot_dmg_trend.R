#' @title Plot chromosomal DMGs trend
#' @description Plot chromosomal DMGs \bold{\emph{trend}}.
#'
#' @return A \bold{\emph{ggplot object}} chromosomal DMGs trend.
#' @param chrom_id Chromosome identifier (e.g., \bold{\emph{"chr1"}}).
#' @param dmr_file DEG table from \bold{\emph{MethylKit}} analysis.
#' @param smooth_span Span for local regression smoothing. (\bold{\emph{0.1}}).
#' @param hyper_color Color for hyper-methylated DMRs. (\bold{\emph{"#ff0000"}}).
#' @param hypo_color Color for hypo-methylated DMRs. (\bold{\emph{"#008800"}}).
#' @param point_size Point size. (\bold{\emph{3}}).
#' @param point_alpha Point alpha. (\bold{\emph{0.5}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example DMR table in GAnnoViz
#' dmr_file <- system.file(
#'   "extdata",
#'   "example.dmr",
#'   package = "GAnnoViz")
#'
#' # Plot DMR trend
#' plot_dmg_trend(
#'   chrom_id = "chr1",
#'   dmr_file = dmr_file,
#'   smooth_span = 0.1,
#'   hyper_color = "#ff0000",
#'   hypo_color = "#008800",
#'   point_size = 3,
#'   point_alpha = 0.5
#' )
plot_dmg_trend <- function(chrom_id,
                           dmr_file,
                           smooth_span = 0.1,
                           hyper_color = "#ff000055",
                           hypo_color = "#00880055",
                           point_size = 3,
                           point_alpha = 0.5) {
  dmr <- utils::read.table(
    dmr_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    stringsAsFactors = FALSE
  )
  dmr <- dmr %>% dplyr::filter(chr == chrom_id)
  if (nrow(dmr) == 0)
    stop("No DMRs on the selected chromosome")

  dmr$state <- ifelse(dmr$meth.diff >= 0, "hyper", "hypo")
  ggplot2::ggplot(dmr, ggplot2::aes(x = start, y = meth.diff, color = state)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_smooth(span = smooth_span,
                         se = FALSE,
                         color = "#444444") +
    ggplot2::scale_color_manual(values = c(hyper = hyper_color, hypo = hypo_color)) +
    ggplot2::labs(x = "Genomic position", y = "Methylation diff (%)", color = NULL) +
    my_theme()
}
