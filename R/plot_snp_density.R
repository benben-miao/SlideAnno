#' @title Plot SNP density at chromosome level
#' @description Plot \bold{\emph{SNP density}} at chromosome level.
#'
#' @return Plot SNP density at chromosome level.
#' @param fst_file FST slide window file (\bold{\emph{tab-separated}}).
#' @param LOG10 Whether to transform values by log10 before plotting. (\bold{\emph{FALSE}}).
#' @param bin_size Window size (bp) for density calculation. (\bold{\emph{1e6}}).
#' @param density_color Color palette for chromosome density. (\bold{\emph{c("#0088ff", "#ff8800", "#ff0000")}}).
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example data in SlideAnno
#' fst_file <- system.file(
#'   "extdata",
#'   "example.fst",
#'   package = "SlideAnno")
#'
#' # Plot SNP density
#' plot_snp_density(
#'   fst_file = fst_file,
#'   LOG10 = FALSE,
#'   bin_size = 1e6,
#'   density_color = c("#0088ff", "#ff8800", "#ff0000")
#' )
#'
plot_snp_density <- function(fst_file,
                             LOG10 = FALSE,
                             bin_size = 1e6,
                             density_color = c("#0088ff", "#ff8800", "#ff0000")) {
  # FST
  fst <- utils::read.table(
    fst_file,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    stringsAsFactors = FALSE
  )

  # Top5 and top1
  top5_top1 <- stats::quantile(fst$WEIGHTED_FST,
                               probs = c(0.95, 0.99),
                               na.rm = TRUE)
  top5_value <- top5_top1[1]
  top1_value <- top5_top1[2]

  fst_top5 <- fst %>%
    dplyr::arrange(dplyr::desc(WEIGHTED_FST)) %>%
    dplyr::filter(WEIGHTED_FST >= top5_value) %>%
    dplyr::mutate(id = paste0(CHROM, "_", BIN_START))

  fst_top1 <- fst %>%
    dplyr::arrange(dplyr::desc(WEIGHTED_FST)) %>%
    dplyr::filter(WEIGHTED_FST >= top1_value) %>%
    dplyr::mutate(id = paste0(CHROM, "_", BIN_START))

  data <- fst_top5 %>%
    dplyr::mutate(CHROM = as.numeric(gsub("chr", "", CHROM))) %>%
    dplyr::select(id, CHROM, BIN_START, WEIGHTED_FST)

  threshold_vals <- c(top5_value, top1_value)

  CMplot::CMplot(
    data,
    LOG10 = LOG10,
    ### chromosome
    bin.size = bin_size,
    band = 1,
    chr.labels = NULL,
    chr.border = FALSE,
    chr.pos.max = FALSE,
    chr.labels.angle = 0,
    chr.den.col = density_color,
    cir.band = 0.1,
    cir.chr = TRUE,
    cir.chr.h = 0.5,
    cir.axis = TRUE,
    cir.axis.col = "black",
    cir.axis.grid = TRUE,
    ### multitrack
    multraits = FALSE,
    multracks = FALSE,
    multracks.xaxis = FALSE,
    H = 1,
    # ylim = c(0.1, 0.8),
    ### points
    pch = 16,
    type = "p",
    # "p" (point), "l" (cross line), "h" (vertical lines)
    col = rep(c("#aaaaaa55", "#33333355"), 9),
    # points.alpha = 100L,
    ### circle
    plot.type = "d",
    # "m","c","q","d"
    cex = c(0.7, 1, 1),
    # circle, manhattan, qqplot
    r = 1,
    # circle
    outward = TRUE,
    # circle
    axis.cex = 1,
    # circle
    axis.lwd = 1,
    # circle
    lab.cex = 1.5,
    # circle
    lab.font = 1,
    # circle
    ### manhattan
    # ylab = expression(-log[10](italic(p))), # manhattan and qqplot
    ylab = "WEIGHTED FST (Top 5%)",
    ylab.pos = 3,
    # manhattan
    xticks.pos = 1,
    # manhattan
    mar = c(0, 0, 0, 0),
    # manhattan, bottom, left, up, and right
    mar.between = 1,
    # manhattan
    ### threshold
    threshold = threshold_vals,
    threshold.col = c("#0000ff", "#ff0000"),
    threshold.lty = c(2, 2),
    threshold.lwd = 2,
    ### signal
    amplify = FALSE,
    signal.cex = 2,
    signal.pch = 16,
    signal.line = 2,
    signal.col = NULL,
    ### highlight
    highlight = NULL,
    highlight.text = NULL,
    highlight.cex = 1,
    highlight.pch = 16,
    highlight.type = "p",
    highlight.col = "#ff0000",
    highlight.text.font = 1,
    # 1: normal, 2: bold, 3: italic
    highlight.text.cex = 0.7,
    highlight.text.col = "#000000",
    ### save
    conf.int = TRUE,
    conf.int.col = "white",
    file.output = FALSE,
    file.name = "file_name",
    file = "pdf",
    # "jpg", "pdf", "tiff", "png"
    box = FALSE,
    main = "",
    main.cex = 1,
    main.font = 1,
    legend.ncol = NULL,
    legend.cex = 1,
    legend.pos = "right",
    # "left","middle","right"
    dpi = 300,
    width = NULL,
    height = NULL,
    verbose = FALSE
  )
}
